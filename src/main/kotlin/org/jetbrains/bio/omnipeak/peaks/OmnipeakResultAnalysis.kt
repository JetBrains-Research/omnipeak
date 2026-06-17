package org.jetbrains.bio.omnipeak.peaks

import com.google.common.collect.MinMaxPriorityQueue
import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.FragmentSize
import org.jetbrains.bio.genome.coverage.PairedEndCoverage
import org.jetbrains.bio.genome.coverage.SingleEndCoverage
import org.jetbrains.bio.omnipeak.OmnipeakCLA
import org.jetbrains.bio.omnipeak.coverage.NormalizedBinnedCoverageQuery
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_FRAGMENTATION_MAX_GAP_BP
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.omnipeak.peaks.AutoCorrelations.computeAutoCorrelations
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.estimateSingleModeLength
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.getChromosomeCandidates
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.getLogNullPvals
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.getLogNulls
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.analyzeAdditiveCandidates
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.detectSensitivityTriangle
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.estimateCandidatesNumberLens
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.estimateGap
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.getSensitivitiesAndCandidatesCharacteristics
import org.jetbrains.bio.omnipeak.peaks.Signal.computeSignalToControlAverage
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.util.deleteIfExists
import org.jetbrains.bio.util.stem
import org.jetbrains.bio.util.toPath
import org.jetbrains.bio.viktor.F64Array
import org.knowm.xchart.BitmapEncoder
import org.knowm.xchart.XYChartBuilder
import org.knowm.xchart.XYSeries
import org.knowm.xchart.style.Styler
import org.knowm.xchart.style.markers.SeriesMarkers
import org.knowm.xchart.style.theme.GGPlot2Theme
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.awt.Color
import java.io.BufferedWriter
import java.io.FileWriter
import java.nio.file.Path
import kotlin.math.*

object OmnipeakResultAnalysis {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun doDeepAnalysis(
        actualModelPath: Path,
        omnipeakFitResults: OmnipeakFitResults,
        fitInfo: OmnipeakAnalyzeFitInformation,
        genomeQuery: GenomeQuery,
        fdr: Double,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        blackListPath: Path?,
        peaksList: List<Peak>,
        peaksPath: Path?
    ) {
        val name = actualModelPath.fileName.stem
        check(fitInfo.binnedCoverageQueries != null) {
            "$name Please use prepareData before!"
        }
        check(fitInfo.binnedCoverageQueries!!.all { it.areCachesPresent() }) {
            "$name Coverage information is not available"
        }
        val infoFile = if (peaksPath != null) "$peaksPath.txt" else null
        infoFile?.toPath()?.deleteIfExists()
        val infoWriter = if (infoFile != null) BufferedWriter(FileWriter(infoFile)) else null

        // Save basic stats to infoFile
        LOG.info("$name Processing basic info")
        if (infoWriter != null) {
            val aboutModel = omnipeakFitResults.modelInformation(actualModelPath)
            val aboutPeaks = PeaksInfo.compute(
                genomeQuery,
                peaksList.map { it.location }.stream(),
                peaksPath!!.toUri(),
                fitInfo.paths.map { it.treatment }
            )
            infoWriter.write((aboutModel + aboutPeaks).joinToString("\n") { (k, v) ->
                "${k.name}: ${k.render(v)}"
            } + "\n")
        }
        val modelLowSignalToNoise = omnipeakFitResults.model is NB2ZHMM &&
                omnipeakFitResults.model.outOfSignalToNoiseRatioRangeDown
        logInfo("Model low signal to noise: $modelLowSignalToNoise", infoWriter)

        val blackList = if (blackListPath != null) {
            LOG.info("$name Loading blacklist regions: $blackListPath")
            LocationsMergingList.load(genomeQuery, blackListPath)
        } else null


        LOG.info("$name Analysing auto correlations...")
        var coverage: DoubleArray
        var coverageCorrelations: DoubleArray = DoubleArray(0)
        var logNullPValsCorrelations: DoubleArray = DoubleArray(0)
        var fragmentSizeCorrelations: DoubleArray? = null
        var detectedFragment = -1
        val binnedQuery = fitInfo.binnedCoverageQueries!!.first()
        if (binnedQuery is NormalizedBinnedCoverageQuery) {
            val ncq = binnedQuery.ncq
            val (controlScale, beta, minCorrelation) = ncq.coveragesNormalizedInfo
            val treatmentCoverage = ncq.treatmentReads.coverage()
            val treatmentTotal = genomeQuery.get().sumOf {
                treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange).toLong()
            }
            logInfo("Treatment coverage: $treatmentTotal", infoWriter)
            val controlCoverage = ncq.controlReads?.coverage()
            if (controlCoverage != null) {
                val controlTotal = genomeQuery.get().sumOf {
                    controlCoverage.getBothStrandsCoverage(it.chromosomeRange).toLong()
                }
                logInfo("Control coverage: $controlTotal", infoWriter)
                logInfo("Control scale: ${"%.3f".format(controlScale)}", infoWriter)
                logInfo("Beta: ${"%.3f".format(beta)}", infoWriter)
                logInfo("Min control correlation: ${"%.3f".format(minCorrelation)}", infoWriter)
            }
            LOG.info("$name Analysing coverage distribution...")
            coverage = computeCoverageScores(
                genomeQuery,
                treatmentCoverage, controlCoverage, controlScale,
                beta, fitInfo.binSize, blackList
            )
            LOG.info("$name Analysing tracks variance...")
            val normVariance = computeAverageVariance(
                genomeQuery, treatmentCoverage, controlCoverage, controlScale, beta, blackList
            )
            logInfo("Track normalized variance: ${"%.3f".format(normVariance)}", infoWriter)

            LOG.debug("$name Analysing coverage autocorrelation...")
            coverageCorrelations = computeAutoCorrelations(coverage)

            if (treatmentCoverage is SingleEndCoverage) {
                LOG.debug("$name Analysing fragment size cross-correlation...")
                fragmentSizeCorrelations = FragmentSize.computePearsonCorrelationTransform(
                    (0..FragmentSize.MAX_FRAGMENT_SIZE).toList(),
                    treatmentCoverage.data
                ).map { it.pearsonTransform }.toDoubleArray()
                detectedFragment = treatmentCoverage.detectedFragment
                logInfo("Detected fragment size: $detectedFragment", infoWriter)
            } else if (treatmentCoverage is PairedEndCoverage) {
                detectedFragment = treatmentCoverage.averageFragmentSize
                logInfo("Average fragment size: $detectedFragment", infoWriter)
            }
        } else {
            coverage = DoubleArray(genomeQuery.get().size) { 0.0 }
        }
        val positiveCoverage = coverage.filter { it > 0 }.toDoubleArray()
        val positivePercentage = if (coverage.isNotEmpty()) (100.0 * positiveCoverage.size / coverage.size).toInt() else 0
        val max = positiveCoverage.maxOrNull() ?: 0.0
        val mean = if (positiveCoverage.isNotEmpty()) positiveCoverage.average() else 0.0
        val median = if (positiveCoverage.isNotEmpty()) StatUtils.percentile(positiveCoverage, 50.0) else 0.0
        val std = if (positiveCoverage.isNotEmpty()) positiveCoverage.standardDeviation() else 0.0
        logInfo(
            "Positive coverage distribution: coverage >0 %: $positivePercentage, " +
                    "max: ${"%.3f".format(max)}, mean: ${"%.3f".format(mean)}, " +
                    "median: ${"%.3f".format(median)}, std: ${"%.3f".format(std)}",
            infoWriter
        )
        coverage = positiveCoverage

        LOG.info("$name Analysing log null pvalues distribution...")
        val logNullPvals = getLogNullPvals(genomeQuery, omnipeakFitResults, blackList)
        logInfo("LogNullPVals mean: ${"%.3f".format(logNullPvals.average())}", infoWriter)
        logInfo("LogNullPVals std: ${"%.3f".format(logNullPvals.standardDeviation())}", infoWriter)

        LOG.debug("$name Analysing log null pvalues autocorrelation...")
        logNullPValsCorrelations = computeAutoCorrelations(logNullPvals)
        val avgAutoCorrelation = logNullPValsCorrelations.average()
        logInfo("Average autocorrelation score: ${"%.3f".format(avgAutoCorrelation)}", infoWriter)

        // Collect candidates from model
        val logNullMembershipsMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                // HACK, it shouldn't be used
                // Cannot be null because of GenomeMap and empty F64Array is not supported
                return@genomeMap F64Array.full(1, 0.0)
            }
            getLogNulls(omnipeakFitResults, chromosome)
        }
        val bitList2reuseMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
            BitList(logNullMembershipsMap[chromosome].length)
        }

        LOG.debug("$name Analysing autocorrelation...")
        // Already computed above

        LOG.info("$name Analysing sensitivity...")
        val (sensitivities, candidatesNs, candidatesALs) =
            getSensitivitiesAndCandidatesCharacteristics(
                genomeQuery,
                omnipeakFitResults,
                logNullMembershipsMap,
                bitList2reuseMap
            )
        val st = detectSensitivityTriangle(sensitivities, candidatesNs, candidatesALs)

        LOG.info("$name Analysing additive candidates...")
        val (totals, news) = analyzeAdditiveCandidates(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
            sensitivities, true
        )

        val sensitivity2use: Double
        when {
            sensitivityCmdArg != null ->
                sensitivity2use = sensitivityCmdArg

            st != null -> {
                val (beforeMerge, stable, beforeNoise) = st
                logInfo("Sensitivity beforeMerge: ${"%.3f".format(sensitivities[beforeMerge])}", infoWriter)
                logInfo("Sensitivity beforeMerge index: $beforeMerge", infoWriter)
                logInfo("Sensitivity stable: ${"%.3f".format(sensitivities[stable])}", infoWriter)
                logInfo("Sensitivity stable index: $stable", infoWriter)
                logInfo("Sensitivity beforeNoise: ${"%.3f".format(sensitivities[beforeNoise])}", infoWriter)
                logInfo("Sensitivity beforeNoise index: $beforeNoise", infoWriter)

                val minAdditionalIdx = (st.beforeMerge until st.stable)
                    .minByOrNull { if (totals[it] == 0) 0.0 else news[it].toDouble() / totals[it].toDouble() }!!
                val minAdditionalSensitivity = sensitivities[minAdditionalIdx]
                logInfo("Minimal additional: ${"%.3f".format(minAdditionalSensitivity)}", infoWriter)
                logInfo("Minimal additional index: $minAdditionalIdx", infoWriter)
                sensitivity2use = minAdditionalSensitivity
            }

            else -> {
                LOG.error("$name Failed to automatically estimate sensitivity")
                sensitivity2use = ln(fdr)
            }
        }
        logInfo("Sensitivity2use: ${"%.3f".format(sensitivity2use)}", infoWriter)

        LOG.debug("$name Analysing gap...")
        val candidateGapNs = IntArray(OMNIPEAK_FRAGMENTATION_MAX_GAP_BP / fitInfo.binSize) {
            estimateCandidatesNumberLens(
                genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                sensitivity2use, it
            ).n
        }
        val gap2use = if (gapCmdArg != null) {
            gapCmdArg
        } else {
            estimateGap(candidateGapNs)
        }

        logInfo("Gap2use: $gap2use", infoWriter)

        val candidatesMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(chromosome, logNullMemberships, bitList2reuse, sensitivity2use, gap2use)
        }

        val avgModeLen = estimateSingleModeLength(genomeQuery, fitInfo, candidatesMap, true)
        logInfo("Single model length: ${"%.3f".format(avgModeLen.toDouble())}", infoWriter)

        val peakScorer = PeakScorer.create(fitInfo, logNullMembershipsMap)

        // Estimate signal and noise average signal by candidates
        val (avgSignalDensity, avgNoiseDensity) =
            peakScorer.analyzeSignalAndNoise(genomeQuery, fitInfo, candidatesMap, true)
        logInfo("Candidates signal density: ${"%.3f".format(avgSignalDensity ?: 0.0)}", infoWriter)
        logInfo("Candidates noise density: ${"%.3f".format(avgNoiseDensity ?: 0.0)}", infoWriter)
        val signalToNoise = if (avgSignalDensity != null && avgNoiseDensity != null)
            avgSignalDensity / avgNoiseDensity else 0.0
        logInfo("Coverage signal to noise: ${"%.3f".format(signalToNoise)}", infoWriter)

        if (peakScorer.coverageControlComputable != null) {
            val signalToControl = computeSignalToControlAverage(
                genomeQuery, fitInfo, candidatesMap, peakScorer.coverageControlComputable, true
            )
            logInfo("Coverage signal to control: ${"%.3f".format(signalToControl)}", infoWriter)
        }
        infoWriter?.close()

        prepareCoveragePercentilesTsvFile(coverage, peaksPath)

        prepareLogNullPercentilesTsvFile(logNullPvals, peaksPath)

        prepareAutocorrelationTsvFile(coverageCorrelations, ".ac.coverage.tsv", peaksPath)

        prepareAutocorrelationTsvFile(logNullPValsCorrelations, ".ac.pvals.tsv", peaksPath)

        prepareGapsFile(peaksPath, candidateGapNs)

        prepareSensitivitiesTsvFile(
            genomeQuery, omnipeakFitResults, logNullMembershipsMap, bitList2reuseMap,
            peaksPath, sensitivities
        )

        prepareDeepAnalysisPlots(
            peaksPath, sensitivities, candidatesNs, candidatesALs, totals, news, st, sensitivity2use,
            coverageCorrelations, logNullPValsCorrelations, fragmentSizeCorrelations, detectedFragment,
            candidateGapNs, gap2use, coverage, positivePercentage
        )

        prepareSegmentsTsvFile(sensitivities, totals, news, peaksPath)
        OmnipeakModelPlots.prepareModelPlots(peaksPath, omnipeakFitResults)
        LOG.info("$name Done analysis")
    }

    private fun prepareAutocorrelationTsvFile(correlations: DoubleArray, suffix: String, peaksPath: Path?) {
        val autocorrelationFile = if (peaksPath != null) "$peaksPath$suffix" else null
        autocorrelationFile?.toPath()?.deleteIfExists()
        if (autocorrelationFile != null) {
            LOG.info("See $autocorrelationFile")
        }

        val autocorrelationWriter = if (autocorrelationFile != null)
            BufferedWriter(FileWriter(autocorrelationFile))
        else
            null
        logInfo("D\tCorrelation", autocorrelationWriter, false)
        correlations.forEachIndexed { i, d ->
            logInfo("${i + 1}\t$d", autocorrelationWriter, false)
        }
        autocorrelationWriter?.close()
    }

    private fun prepareCoveragePercentilesTsvFile(
        coverage: DoubleArray,
        peaksPath: Path?
    ) {
        LOG.info("Analysing coverage percentiles...")
        val coveragePercFile = if (peaksPath != null) "$peaksPath.coverage.tsv" else null
        coveragePercFile?.toPath()?.deleteIfExists()
        if (coveragePercFile != null) {
            LOG.info("See $coveragePercFile")
        }

        val coveragePercWriter = if (coveragePercFile != null)
            BufferedWriter(FileWriter(coveragePercFile))
        else
            null
        logInfo("Q\tCoverage", coveragePercWriter, false)
        for (i in 0 until 100) {
            logInfo("${i + 1}\t${StatUtils.percentile(coverage, i + 1.0)}", coveragePercWriter, false)
        }
        coveragePercWriter?.close()
    }

    private fun prepareLogNullPercentilesTsvFile(
        logNullPvals: DoubleArray,
        peaksPath: Path?,
        maxQ: Double = 0.01,
        step: Double = 1e-5,
    ) {
        LOG.info("Analysing log nulls percentiles...")
        logNullPvals.sort()
        val logNullPsFile = if (peaksPath != null) "$peaksPath.logps.tsv" else null
        logNullPsFile?.toPath()?.deleteIfExists()
        if (logNullPsFile != null) {
            LOG.info("See $logNullPsFile")
        }

        val logNullPsWriter = if (logNullPsFile != null)
            BufferedWriter(FileWriter(logNullPsFile))
        else
            null
        logInfo("Q\tLogNullP", logNullPsWriter, false)
        var realStep = step
        while (1 / step > logNullPvals.size) {
            realStep *= 10
        }
        var q = 0.0
        while (q < maxQ) {
            logInfo("$q\t${logNullPvals[(logNullPvals.size * q).toInt()]}", logNullPsWriter, false)
            q += realStep
        }
        logNullPsWriter?.close()
    }

    private fun prepareGapsFile(
        peaksPath: Path?,
        candidateGapNs: IntArray,
    ) {
        LOG.info("Analysing candidates wrt gap...")
        val gapsDetailsFile = if (peaksPath != null) "$peaksPath.gaps.tsv" else null
        if (gapsDetailsFile != null) {
            LOG.info("See $gapsDetailsFile")
        }
        gapsDetailsFile?.toPath()?.deleteIfExists()
        val gapsDetailsWriter = if (gapsDetailsFile != null)
            BufferedWriter(FileWriter(gapsDetailsFile))
        else
            null
        logInfo(
            "Gap\tCandidatesN",
            gapsDetailsWriter, false
        )
        candidateGapNs.forEachIndexed { g, n ->
            logInfo("$g\t$n", gapsDetailsWriter, false)
        }
        gapsDetailsWriter?.close()
    }

    private fun prepareSensitivitiesTsvFile(
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        peaksPath: Path?,
        sensitivities: DoubleArray,
    ) {
        LOG.info("Analysing candidates characteristics wrt sensitivity and gap...")
        val sensDetailsFile = if (peaksPath != null) "$peaksPath.sensitivity.tsv" else null
        if (sensDetailsFile != null) {
            LOG.info("See $sensDetailsFile")
        }
        sensDetailsFile?.toPath()?.deleteIfExists()
        val sensDetailsWriter = if (sensDetailsFile != null)
            BufferedWriter(FileWriter(sensDetailsFile))
        else
            null
        logInfo(
            "Sensitivity\tCandidatesN\tCandidatesAL\tCandidatesML",
            sensDetailsWriter, false
        )
        for (s in sensitivities) {
            val candidatesMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
                if (!omnipeakFitResults.fitInfo.containsChromosomeInfo(chromosome)) {
                    return@genomeMap emptyList<Range>()
                }
                val logNullMemberships = logNullMembershipsMap[chromosome]
                val bitList2reuse = bitList2reuseMap[chromosome]
                getChromosomeCandidates(chromosome, logNullMemberships, bitList2reuse, s, 0)
            }
            val candidatesList = genomeQuery.get().flatMap { chromosome ->
                candidatesMap[chromosome].map {
                    Location(it.startOffset, it.endOffset, chromosome)
                }
            }
            val total = candidatesList.size
            val lengths = DoubleArray(total) { candidatesList[it].length().toDouble() }
            val avgL = lengths.average()
            val medianL = StatUtils.percentile(lengths, 50.0)
            logInfo("$s\t$total\t$avgL\t$medianL", sensDetailsWriter, false)
        }
        sensDetailsWriter?.close()
    }

    private fun prepareSegmentsTsvFile(
        sensitivities: DoubleArray,
        totals: IntArray,
        news: IntArray,
        peaksPath: Path?,
    ) {
        LOG.info("Analysing segments...")
        val segmentsFile = if (peaksPath != null) "$peaksPath.segments.tsv" else null
        if (segmentsFile != null) {
            LOG.info("See $segmentsFile")
        }
        segmentsFile?.toPath()?.deleteIfExists()
        val segmentsWriter = if (segmentsFile != null)
            BufferedWriter(FileWriter(segmentsFile))
        else
            null
        logInfo(
            "Sensitivity\tCandidatesN\tNew\tOld",
            segmentsWriter, false
        )
        sensitivities.forEachIndexed { i, s ->
            val total = totals[i]
            val new = news[i]
            val old = total - new
            logInfo(
                "$s\t$total\t$new\t$old",
                segmentsWriter, false
            )
        }
        segmentsWriter?.close()
    }

    private fun prepareDeepAnalysisPlots(
        peaksPath: Path?,
        sensitivities: DoubleArray,
        candidatesNs: IntArray,
        candidatesALs: DoubleArray,
        totals: IntArray,
        news: IntArray,
        st: SensitivityGap.PepInfo?,
        sensitivity2use: Double,
        coverageCorrelations: DoubleArray,
        logNullPValsCorrelations: DoubleArray,
        fragmentSizeCorrelations: DoubleArray?,
        detectedFragment: Int,
        candidateGapNs: IntArray,
        gap2use: Int,
        coverage: DoubleArray,
        positivePercentage: Int
    ) {
        if (peaksPath == null) return
        LOG.info("Plotting deep analysis plots...")

        val pep = exp(sensitivity2use)
        if (coverage.isNotEmpty()) {
            val max = coverage.maxOrNull() ?: 0.0
            val mean = coverage.average()
            val median = StatUtils.percentile(coverage, 50.0)
            val std = coverage.standardDeviation()

            val title = "Positive coverage distribution: coverage >0 %: $positivePercentage, " +
                    "max: ${"%.3f".format(max)}, mean: ${"%.3f".format(mean)}, " +
                    "median: ${"%.3f".format(median)}, std: ${"%.3f".format(std)}"

            val xData = DoubleArray(100) { (it + 1).toDouble() }
            val yData = DoubleArray(100) { StatUtils.percentile(coverage, it + 1.0) }

            savePlot(
                peaksPath, "coverage_distribution", title,
                "Percentile", "Coverage",
                xData, yData, null, -1
            )
        }

        val logNs = DoubleArray(candidatesNs.size) { ln1p(candidatesNs[it].toDouble()) }
        val logALs = DoubleArray(candidatesALs.size) { ln1p(candidatesALs[it]) }
        val newPercentages = DoubleArray(sensitivities.size) {
            if (totals[it] == 0) 0.0 else news[it].toDouble() / totals[it] * 100.0
        }

        val pepRanks = DoubleArray(sensitivities.size) { (it + 1).toDouble() }
        val estimatedIdx = sensitivities.indices.minByOrNull { abs(sensitivities[it] - sensitivity2use) } ?: -1
        val estimatedPEPLabel = "Estimated PEP: ${"%.3f".format(pep)}"

        // 1) log Number of candidates vs log Average lengths
        savePlot(
            peaksPath, "logN_vs_logAL", "log Number of candidates vs log Average lengths",
            "log(Number of candidates + 1)", "log(Average lengths + 1)",
            logNs, logALs, st, estimatedIdx,
            estimatedLabel = estimatedPEPLabel
        )

        // 2) PEP vs log number of candidates
        savePlot(
            peaksPath, "pep_vs_logN", "PEP rank vs log number of candidates",
            "PEP rank", "log(Number of candidates + 1)",
            pepRanks, logNs, st, estimatedIdx,
            estimatedLabel = estimatedPEPLabel
        )

        // 3) PEP vs log average lengths
        savePlot(
            peaksPath, "pep_vs_logAL", "PEP rank vs log average lengths",
            "PEP rank", "log(Average lengths + 1)",
            pepRanks, logALs, st, estimatedIdx,
            estimatedLabel = estimatedPEPLabel
        )

        // 4) PEP vs new candidates%
        savePlot(
            peaksPath, "pep_vs_new", "PEP rank vs new candidates%",
            "PEP rank", "New candidates%",
            pepRanks, newPercentages, st, estimatedIdx,
            estimatedLabel = estimatedPEPLabel
        )

        // Also update the original sensitivity vs candidates plot with markers
        val distances = DoubleArray(coverageCorrelations.size) { (it + 1).toDouble() }
        savePlot(
            peaksPath, "ac_coverage", "Autocorrelation: Coverage",
            "Distance", "Correlation",
            distances, coverageCorrelations, null, -1
        )
        savePlot(
            peaksPath, "ac_pvals", "Autocorrelation: LogNullPVals",
            "Distance", "Correlation",
            distances, logNullPValsCorrelations, null, -1
        )

        if (fragmentSizeCorrelations != null) {
            val fragments = DoubleArray(fragmentSizeCorrelations.size) { it.toDouble() }
            savePlot(
                peaksPath, "fragment_size", "Fragment Size Estimation",
                "Fragment Size", "Cross-correlation",
                fragments, fragmentSizeCorrelations, null, detectedFragment,
                estimatedLabel = "Estimated Fragment: $detectedFragment"
            )
        }

        val gapNs = DoubleArray(candidateGapNs.size) { candidateGapNs[it].toDouble() }
        val gaps = DoubleArray(candidateGapNs.size) { it.toDouble() }
        savePlot(
            peaksPath, "gap_candidates", "Gap Candidates Estimation",
            "Gap", "Number of Candidates",
            gaps, gapNs, null, gap2use,
            estimatedLabel = "Estimated Gap: $gap2use"
        )
    }

    private fun savePlot(
        peaksPath: Path,
        name: String,
        title: String,
        xTitle: String,
        yTitle: String,
        xData: DoubleArray,
        yData: DoubleArray,
        st: SensitivityGap.PepInfo?,
        estimatedIdx: Int,
        estimatedLabel: String = "Estimated PEP"
    ) {
        val plotFile = "$peaksPath.$name.png"
        val chart = XYChartBuilder()
            .width(800).height(600)
            .title(title)
            .xAxisTitle(xTitle)
            .yAxisTitle(yTitle)
            .build()
        chart.styler.theme = GGPlot2Theme()
        chart.styler.plotBackgroundColor = Color.WHITE
        chart.styler.chartBackgroundColor = Color.WHITE
        chart.styler.plotGridLinesColor = Color.LIGHT_GRAY
        chart.styler.legendPosition = Styler.LegendPosition.InsideNE
        chart.styler.defaultSeriesRenderStyle = XYSeries.XYSeriesRenderStyle.Line
        chart.styler.xAxisDecimalPattern = if (xTitle in listOf("PEP rank", "Counts", "Gap", "Fragment Size", "Number of Candidates", "Percentile", "Distance")) "0" else "0.000"
        chart.styler.yAxisDecimalPattern = if (yTitle in listOf("PEP rank", "Counts", "Gap", "Fragment Size", "Number of Candidates", "Percentile", "Distance")) "0" else "0.000"

        val series = chart.addSeries("Data", xData, yData)
        series.lineColor = Color.BLACK
        series.markerColor = Color.BLACK
        series.marker = SeriesMarkers.NONE

        if (st != null) {
            val triangleIndices = listOf(st.beforeMerge, st.stable, st.beforeNoise)
            val tx = triangleIndices.map { xData[it] }.toDoubleArray()
            val ty = triangleIndices.map { yData[it] }.toDoubleArray()
            val series = chart.addSeries("Triangle PEPs", tx, ty)
            series.xySeriesRenderStyle = XYSeries.XYSeriesRenderStyle.Scatter
            series.marker = SeriesMarkers.TRIANGLE_UP
            series.markerColor = Color.RED
        }

        if (estimatedIdx != -1 && estimatedIdx < xData.size) {
            val ex = doubleArrayOf(xData[estimatedIdx])
            val ey = doubleArrayOf(yData[estimatedIdx])
            val series = chart.addSeries(estimatedLabel, ex, ey)
            series.xySeriesRenderStyle = XYSeries.XYSeriesRenderStyle.Scatter
            series.marker = SeriesMarkers.DIAMOND
            series.markerColor = Color.BLUE
        }

        try {
            BitmapEncoder.saveBitmapWithDPI(chart, plotFile, BitmapEncoder.BitmapFormat.PNG, 300)
            LOG.info("See $plotFile")
        } catch (e: Exception) {
            LOG.error("Failed to save plot $plotFile", e)
        }
    }

    private fun logInfo(msg: String, infoWriter: BufferedWriter?, useLog: Boolean = true) {
        if (useLog)
            OmnipeakCLA.LOG.info(msg)
        if (infoWriter != null) {
            infoWriter.write(msg)
            infoWriter.newLine()
        }
    }

    const val REGION_LEN = 10_000
    const val TOP_REGIONS = 10_000
    const val WORK_REGIONS = 200
    const val RESOLUTION = 100

    /**
     * Detects normalized variance as average std / mean for candidates
     */

    private fun computeAverageVariance(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        blackList: LocationsMergingList? = null,
        regionLen: Int = REGION_LEN,
        topRegions: Int = TOP_REGIONS,
        workRegions: Int = WORK_REGIONS,
        resolution: Int = RESOLUTION
    ): Double {
        LOG.debug("Compute coverage in regions")
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange) > 0 }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val genomeRegions = limitedQuery.get().sumOf { floor(it.length.toDouble() / regionLen).toLong() }
        check(topRegions < genomeRegions) {
            "Too many top regions $topRegions > $genomeRegions"
        }

        val comparator = Comparator<Triple<Chromosome, Int, Double>> { o1, o2 -> o2.third.compareTo(o1.third) }
        val regionCoverages: MinMaxPriorityQueue<Triple<Chromosome, Int, Double>> =
            MinMaxPriorityQueue
                .orderedBy(comparator)
                .maximumSize(topRegions)
                .create()
        var regions = 0
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            for (i in 0 until floor(chr.length.toDouble() / regionLen).toInt()) {
                regions += 1
                val start = regionLen * i
                val end = regionLen * (i + 1)
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                val c = coverage(
                    chr, start, end,
                    treatmentCoverage, controlCoverage, controlScale, beta
                )
                regionCoverages.add(Triple(chr, i, c))
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, regions)
        }

        val regionCoveragesArray = regionCoverages.toTypedArray()
        regionCoveragesArray.sortWith(comparator)

        val step = if (regionCoveragesArray.size > workRegions) {
            LOG.debug("Pick $workRegions / $topRegions uniform regions for computation speedup")
            ceil(regionCoveragesArray.size.toDouble() / workRegions).toInt()
        } else
            1

        val stdMeans = DoubleArray(workRegions)
        for (i in regionCoveragesArray.indices) {
            if (i % step != 0) {
                continue
            }
            val (chr, start, _) = regionCoveragesArray[i]
            val stats = DoubleArray(regionLen / resolution) {
                coverage(
                    chr, start + it * resolution, start + (it + 1) * resolution,
                    treatmentCoverage, controlCoverage, controlScale, beta
                )
            }
            val mean = stats.average()
            val std = stats.standardDeviation()
            stdMeans[i / step] = if (mean > 0) std / mean else 0.0
        }
        return stdMeans.average()
    }


    private fun computeCoverageScores(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        bin: Int,
        blackList: LocationsMergingList?
    ): DoubleArray {
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange) > 0 }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val totalBins = limitedQuery.get().sumOf { floor(it.length.toDouble() / bin).toInt() }
        val coverage = DoubleArray(totalBins) { 0.0 }
        var i = 0
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            for (j in 0 until floor(chr.length.toDouble() / bin).toInt()) {
                val start = j * bin
                val end = (j + 1) * bin
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                coverage[i++] = coverage(
                    chr, start, end,
                    treatmentCoverage, controlCoverage, controlScale, beta
                )
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, totalBins)
        }
        return coverage
    }


    private fun coverage(
        chromosome: Chromosome,
        start: Int, end: Int,
        treatmentCoverage: Coverage,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double
    ): Double {
        val chromosomeRange = ChromosomeRange(start, end, chromosome)
        val tc = treatmentCoverage.getBothStrandsCoverage(chromosomeRange).toDouble()
        return if (controlCoverage != null && controlScale != null) {
            val cc = controlCoverage.getBothStrandsCoverage(chromosomeRange) * controlScale
            max(0.0, tc - beta * cc)
        } else {
            tc
        }
    }

    fun DoubleArray.standardDeviation(): Double {
        var sum = 0.0
        var sumSq = 0.0
        for (value in this) {
            sum += value
            sumSq += value * value
        }
        return sqrt((sumSq - sum * sum / size) / size)
    }

}

