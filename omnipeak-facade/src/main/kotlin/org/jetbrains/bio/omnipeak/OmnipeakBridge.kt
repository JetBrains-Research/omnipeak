package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.omnipeak.peaks.MultipleTesting
import org.jetbrains.bio.omnipeak.peaks.OmnipeakResult
import org.jetbrains.bio.omnipeak.peaks.Peak
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import java.nio.file.Path

/**
 * Architectural facade for the Omnipeak runtime.
 *
 * Defines the public ABI contract between JBR and omnipeak.jar,
 * eliminating direct coupling to internal static helper functions and
 * enabling automated signature verification.
 */
interface OmnipeakBridge {

    fun isLowMemory(): Boolean

    fun recommendedMemoryGb(): Int

    fun loadFitResults(genomeQuery: GenomeQuery, modelPath: Path): OmnipeakFitResults

    fun modelInformation(fitResults: OmnipeakFitResults, modelPath: Path): List<TrackAboutMetricValue<*>>

    fun fitInformationDataPaths(fitInfo: OmnipeakFitInformation): List<Pair<Path, Path?>>

    fun getModelMeanAndFailures(fitResults: OmnipeakFitResults, state: Int): Pair<Double, Double>?

    fun getModelPriorsAndTransitions(fitResults: OmnipeakFitResults): Pair<F64Array, F64Array>?

    fun getLogNulls(fitResults: OmnipeakFitResults, chromosome: Chromosome): F64Array

    fun estimateSensitivity(
        genomeQuery: GenomeQuery,
        fitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        parallel: Boolean = true,
        name: String? = null,
        cancellableState: CancellableState = CancellableState.current()
    ): Double?

    fun estimateCandidateGapNs(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        maxGap: Int,
        parallel: Boolean = true
    ): IntArray

    fun estimateGapFromCandidateNs(candidatesNs: IntArray): Int

    fun estimateGap(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        maxGap: Int,
        parallel: Boolean = true
    ): Pair<Int, IntArray>

    fun estimateSensitivityAndGap(
        genomeQuery: GenomeQuery,
        fitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        fdr: Double,
        parallel: Boolean = true,
        name: String? = null,
        cancellableState: CancellableState = CancellableState.current()
    ): Pair<Double, Int>

    fun getChromosomeCandidates(
        chromosome: Chromosome,
        logNullMemberships: F64Array,
        bitList2reuse: BitList,
        sensitivity: Double,
        gap: Int
    ): List<Range>

    fun getChromosomeCandidates(
        chromosome: Chromosome,
        logNullMemberships: F64Array,
        bitList2reuse: BitList,
        sensitivities: DoubleArray,
        gap: Int
    ): List<Range>

    fun candidatesToPeaks(
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
        failures: Double? = null,
        multipleTesting: MultipleTesting,
        blacklistPath: Path?,
        name: String?,
        parallel: Boolean,
        cancellableState: CancellableState
    ): List<Peak>

    fun getPeaks(
        fitResults: OmnipeakFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        multipleTesting: MultipleTesting = MultipleTesting.BH,
        sensitivity: Double? = null,
        gap: Int? = null,
        summits: Boolean = false,
        clip: Double = 0.4,
        blackListPath: Path? = null,
        name: String? = null,
        parallel: Boolean = true,
        cancellableState: CancellableState = CancellableState.current()
    ): OmnipeakResult

    fun callPeaks(
        genomeQuery: GenomeQuery,
        paths: List<Pair<Path, Path?>>,
        bin: Int,
        fdr: Double,
        modelPath: Path,
        peaksPath: Path?,
        cancellableState: CancellableState = CancellableState.current()
    ): List<Peak>

    fun savePeaks(peaks: List<Peak>, path: Path, peakName: String = "peak")

    fun processBlackList(genomeQuery: GenomeQuery, peaks: OmnipeakResult, blacklistPath: Path?): List<Peak>

    fun avgAutocorrelation(fitResults: OmnipeakFitResults, fitInfo: OmnipeakFitInformation): Double

    fun estimateSingleModeLength(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        candidatesMap: GenomeMap<List<Range>>,
        parallel: Boolean = false
    ): Int

    fun getChromosomeSummits(
        chromosome: Chromosome,
        fitInfo: OmnipeakFitInformation,
        candidates: List<Range>,
        singleModeLength: Int
    ): List<Range>

    fun kdeSmooth(signal: DoubleArray, bandwidth: Int): DoubleArray


    companion object {
        const val BRIDGE_CLASS_NAME = "org.jetbrains.bio.omnipeak.OmnipeakBridgeImpl"
    }
}
