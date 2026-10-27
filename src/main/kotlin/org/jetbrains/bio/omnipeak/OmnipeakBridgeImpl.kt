package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.omnipeak.fit.*
import org.jetbrains.bio.omnipeak.modes.KdeModeFinder
import org.jetbrains.bio.omnipeak.peaks.AutoCorrelations
import org.jetbrains.bio.omnipeak.peaks.MultipleTesting
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks
import org.jetbrains.bio.omnipeak.peaks.OmnipeakResult
import org.jetbrains.bio.omnipeak.peaks.Peak
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap
import org.jetbrains.bio.omnipeak.peaks.savePeaks
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.stem
import org.jetbrains.bio.viktor.F64Array
import java.nio.file.Path

/**
 * Default implementation of [OmnipeakBridge] residing inside omnipeak.jar.
 */
class OmnipeakBridgeImpl : OmnipeakBridge {

    override fun isLowMemory(): Boolean = OmnipeakCLA.isLowMemory()

    override fun recommendedMemoryGb(): Int = OmnipeakCLA.RECOMMENDED_MEMORY_GB

    override fun loadFitResults(genomeQuery: GenomeQuery, modelPath: Path): OmnipeakFitResults =
        OmnipeakModelFitExperiment.loadResults(genomeQuery, modelPath)

    override fun modelInformation(fitResults: OmnipeakFitResults, modelPath: Path): List<TrackAboutMetricValue<*>> {
        val model = fitResults.model
        return when (model) {
            is NB2ZHMM -> {
                val outOfSnrHitDown = model.outOfSignalToNoiseRatioRangeDown
                val outOfLowerNoise = model.outOfLowerNoise
                val statesSwitched = model.statesSwitched
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                listOf(
                    OmnipeakFitResults.CT_MODEL_FILE to modelPath,
                    OmnipeakFitResults.CT_MODEL_TYPE to "Negative binomial HMM 2 states with zero inflation",
                    OmnipeakFitResults.CT_SIGNAL_MEAN to signalMean,
                    OmnipeakFitResults.CT_NOISE_MEAN to noiseMean,
                    OmnipeakFitResults.CT_SIGNAL_TO_NOISE to ((signalMean + 1e-10) / (noiseMean + 1e-10)),
                    OmnipeakFitResults.CT_OUT_OF_SNR_DOWN to outOfSnrHitDown,
                    OmnipeakFitResults.CT_STATES_SWITCHED to statesSwitched,
                    OmnipeakFitResults.CT_OUT_OF_NOISE_DOWN to outOfLowerNoise
                )
            }
            else -> fitResults.modelInformation(modelPath)
        }
    }

    override fun fitInformationDataPaths(fitInfo: OmnipeakFitInformation): List<Pair<Path, Path?>> =
        when (fitInfo) {
            is OmnipeakAnalyzeFitInformation -> fitInfo.paths.map { it.treatment to it.control }
            is OmnipeakCompareFitInformation -> (fitInfo.data1 + fitInfo.data2).map { it.treatment to it.control }
            else -> emptyList()
        }

    override fun getModelMeanAndFailures(fitResults: OmnipeakFitResults, state: Int): Pair<Double, Double>? {
        val model = fitResults.model
        return if (model is NB2ZHMM) {
            model.means[state] to model.failures[state]
        } else {
            null
        }
    }

    override fun getModelPriorsAndTransitions(fitResults: OmnipeakFitResults): Pair<F64Array, F64Array>? {
        val model = fitResults.model
        return if (model is NB2ZHMM) {
            model.priorProbabilities to model.transitionProbabilities
        } else {
            null
        }
    }

    override fun getLogNulls(fitResults: OmnipeakFitResults, chromosome: Chromosome): F64Array =
        OmnipeakModelToPeaks.getLogNulls(fitResults, chromosome)

    override fun estimateSensitivity(
        genomeQuery: GenomeQuery,
        fitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        parallel: Boolean,
        name: String?,
        cancellableState: CancellableState
    ): Double? = SensitivityGap.estimateSensitivity(
        genomeQuery,
        fitResults,
        logNullMembershipsMap,
        bitList2reuseMap,
        parallel,
        name,
        cancellableState
    ).minAdditionalSensitivity

    override fun estimateCandidateGapNs(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        maxGap: Int,
        parallel: Boolean
    ): IntArray = SensitivityGap.estimateCandidateGapNs(
        genomeQuery,
        fitInfo,
        logNullMembershipsMap,
        bitList2reuseMap,
        sensitivity,
        maxGap,
        parallel
    )

    override fun estimateGapFromCandidateNs(candidatesNs: IntArray): Int =
        SensitivityGap.estimateGap(candidatesNs)

    override fun estimateGap(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        maxGap: Int,
        parallel: Boolean
    ): Pair<Int, IntArray> = SensitivityGap.estimateGap(
        genomeQuery,
        fitInfo,
        logNullMembershipsMap,
        bitList2reuseMap,
        sensitivity,
        maxGap,
        parallel
    )

    override fun estimateSensitivityAndGap(
        genomeQuery: GenomeQuery,
        fitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        fdr: Double,
        parallel: Boolean,
        name: String?,
        cancellableState: CancellableState
    ): Pair<Double, Int> {
        val res = SensitivityGap.estimateSensitivityAndGap(
            genomeQuery, fitResults, logNullMembershipsMap, bitList2reuseMap,
            fdr, null, null, false, parallel, name, cancellableState
        )
        return res.sensitivity to res.gap
    }

    override fun getChromosomeCandidates(
        chromosome: Chromosome,
        logNullMemberships: F64Array,
        bitList2reuse: BitList,
        sensitivity: Double,
        gap: Int
    ): List<Range> = OmnipeakModelToPeaks.getChromosomeCandidates(
        chromosome,
        logNullMemberships,
        bitList2reuse,
        sensitivity,
        gap
    )

    override fun getChromosomeCandidates(
        chromosome: Chromosome,
        logNullMemberships: F64Array,
        bitList2reuse: BitList,
        sensitivities: DoubleArray,
        gap: Int
    ): List<Range> = OmnipeakModelToPeaks.getChromosomeCandidates(
        chromosome,
        logNullMemberships,
        bitList2reuse,
        sensitivities,
        gap
    )

    override fun candidatesToPeaks(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        candidatesMap: GenomeMap<List<Range>>,
        sensitivity: Double,
        avgAutocorrelation: Double,
        summits: Boolean,
        singleModeLength: Int,
        clip: Double,
        fdr: Double,
        failures: Double?,
        multipleTesting: MultipleTesting,
        blacklistPath: Path?,
        name: String?,
        parallel: Boolean,
        cancellableState: CancellableState
    ): List<Peak> = OmnipeakModelToPeaks.candidatesToPeaks(
        genomeQuery, fitInfo, logNullMembershipsMap, avgAutocorrelation,
        failures = failures, sensitivity, candidatesMap, summits, singleModeLength,
        fdr, clip, testing = multipleTesting, blackListPath = blacklistPath, name = name, parallel = parallel,
        cancellableState = cancellableState
    )

    override fun getPeaks(
        fitResults: OmnipeakFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        multipleTesting: MultipleTesting,
        sensitivity: Double?,
        gap: Int?,
        summits: Boolean,
        clip: Double,
        blackListPath: Path?,
        name: String?,
        parallel: Boolean,
        cancellableState: CancellableState
    ): OmnipeakResult = OmnipeakModelToPeaks.getPeaks(
        fitResults, genomeQuery, fdr, multipleTesting,
        sensitivity, gap, summits, clip, blackListPath, parallel, name,
        cancellableState
    )

    override fun callPeaks(
        genomeQuery: GenomeQuery,
        paths: List<Pair<Path, Path?>>,
        bin: Int,
        fdr: Double,
        modelPath: Path,
        peaksPath: Path?,
        cancellableState: CancellableState
    ): List<Peak> {
        val experiment = OmnipeakPeakCallingExperiment.getExperiment(
            genomeQuery,
            paths.map { OmnipeakDataPaths(it.first, it.second) },
            null,
            Fragment.fromString("auto"),
            true,
            bin,
            true,
            OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR,
            OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD,
            modelPath,
            OmnipeakConstants.OMNIPEAK_DEFAULT_FIT_THRESHOLD,
            OmnipeakConstants.OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo = false,
            keepCacheFiles = false
        )
        cancellableState.checkCanceled()
        val results = experiment.results
        cancellableState.checkCanceled()
        val peaks = OmnipeakModelToPeaks.getPeaks(
            results,
            genomeQuery,
            fdr,
            OmnipeakConstants.OMNIPEAK_DEFAULT_MULTIPLE_TEST_CORRECTION,
            null,
            null,
            false,
            clip = OmnipeakConstants.OMNIPEAK_DEFAULT_CLIP_MAX_SIGNAL,
            blackListPath = null,
            name = modelPath.fileName.stem,
            cancellableState = cancellableState
        )
        cancellableState.checkCanceled()
        val peaksList = OmnipeakCLAAnalyze.processBlackList(genomeQuery, peaks, null)
        if (peaksPath != null) {
            Peak.savePeaks(peaksList, peaksPath)
        }
        return peaksList
    }

    override fun savePeaks(peaks: List<Peak>, path: Path, peakName: String) {
        Peak.savePeaks(peaks, path, peakName)
    }

    override fun processBlackList(genomeQuery: GenomeQuery, peaks: OmnipeakResult, blacklistPath: Path?): List<Peak> =
        OmnipeakCLAAnalyze.processBlackList(genomeQuery, peaks, blacklistPath)

    override fun estimateSingleModeLength(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        candidatesMap: GenomeMap<List<Range>>,
        parallel: Boolean
    ): Int = OmnipeakModelToPeaks.estimateSingleModeLength(
        genomeQuery,
        fitInfo,
        candidatesMap,
        parallel
    )

    override fun getChromosomeSummits(
        chromosome: Chromosome,
        fitInfo: OmnipeakFitInformation,
        candidates: List<Range>,
        singleModeLength: Int
    ): List<Range> = OmnipeakModelToPeaks.getChromosomeSummits(
        chromosome, fitInfo, candidates, singleModeLength
    )

    override fun kdeSmooth(signal: DoubleArray, bandwidth: Int): DoubleArray =
        KdeModeFinder.kdeSmooth(signal, bandwidth)

    override fun avgAutocorrelation(fitResults: OmnipeakFitResults, fitInfo: OmnipeakFitInformation): Double =
        AutoCorrelations.avgAutocorrelation(fitResults, fitInfo)

}
