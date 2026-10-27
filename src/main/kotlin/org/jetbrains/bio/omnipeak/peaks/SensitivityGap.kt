package org.jetbrains.bio.omnipeak.peaks

import gnu.trove.list.array.TDoubleArrayList
import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.containers.toRangeMergingList
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FDR
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_GAP
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_FRAGMENTATION_MAX_GAP_BP
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_MIN_SENSITIVITY
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_SENSITIVITY_N
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.getChromosomeCandidates
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.await
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import java.util.concurrent.Callable
import kotlin.math.*

object SensitivityGap {
    private val LOG = LoggerFactory.getLogger(SensitivityGap::class.java)

    data class PepInfo(
        val beforeMerge: Int, val stable: Int, val beforeNoise: Int,
    )

    data class CandidatesInfo(
        val n: Int, val averageLen: Double, val medianLen: Double, val maxLen: Int
    )

    @Suppress("ArrayInDataClass")
    data class SensitivityInfo(
        val sensitivities: DoubleArray,
        val candidatesNs: IntArray,
        val candidatesALs: DoubleArray,
        val pepInfo: PepInfo? = null,
        val totals: IntArray? = null,
        val news: IntArray? = null,
        val minAdditionalIdx: Int? = null,
        val minAdditionalSensitivity: Double? = null,
    )

    /**
     * Result of sensitivity and gap estimation.
     *
     * @property sensitivity sensitivity for candidate selection.
     * @property gap gap size used to merge adjacent candidates.
     * @property avgAutocorrelation broad index (autocorrelation of posterior error probabilities).
     * @property sensitivityInfo detailed information collected during sensitivity estimation.
     * @property candidateGapNs candidate counts evaluated across tested gap sizes.
     */
    @Suppress("ArrayInDataClass")
    data class SensitivityGapResult(
        val sensitivity: Double,
        val gap: Int,
        val avgAutocorrelation: Double = 0.0,
        val sensitivityInfo: SensitivityInfo = SensitivityInfo(doubleArrayOf(), intArrayOf(), doubleArrayOf()),
        val candidateGapNs: IntArray = intArrayOf()
    )

    fun estimateSensitivityAndGap(
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        fdr: Double = OMNIPEAK_DEFAULT_FDR,
        sensitivityCmdArg: Double? = null,
        gapCmdArg: Int? = null,
        summits: Boolean = false,
        parallel: Boolean = true,
        name: String? = null,
        cancellableState: CancellableState? = null
    ): SensitivityGapResult {
        val fitInfo = omnipeakFitResults.fitInfo
        val avgAutocorrelation = AutoCorrelations.avgAutocorrelation(omnipeakFitResults, fitInfo, genomeQuery)
        LOG.info("${name?.let { "$it " } ?: ""}Broad index (average PEP autocorrelation): ${
            "%.3f".format(
                avgAutocorrelation
            )
        }")
        // Estimate sensitivity
        val sensitivityInfo = estimateSensitivity(
            genomeQuery, omnipeakFitResults, logNullMembershipsMap, bitList2reuseMap,
            parallel, name, cancellableState
        )
        val beforeMergeIdx = if (sensitivityCmdArg == null) sensitivityInfo.pepInfo?.beforeMerge else null
        val mergeSensitivity = if (beforeMergeIdx != null && beforeMergeIdx in sensitivityInfo.sensitivities.indices)
            sensitivityInfo.sensitivities[beforeMergeIdx] else null
        if (mergeSensitivity != null) {
            LOG.debug("${name?.let { "$it " } ?: ""}Estimated merge sensitivity: ${"%.3f".format(mergeSensitivity)}")
        }
        val minAdditionalSensitivity =
            if (sensitivityCmdArg == null) sensitivityInfo.minAdditionalSensitivity else null
        if (minAdditionalSensitivity != null) {
            LOG.debug("${name?.let { "$it " } ?: ""}Estimated stable sensitivity ${
                "%.3f".format(minAdditionalSensitivity)}")
        }
        val logFdr = ln(fdr)
        val sensitivity = when {
            sensitivityCmdArg != null -> {
                LOG.info("${name?.let { "$it " } ?: ""}Using command line sensitivity")
                sensitivityCmdArg
            }

            minAdditionalSensitivity != null -> {
                val mergedSensitivity = computeMergedSensitivity(
                    avgAutocorrelation,
                    minAdditionalSensitivity,
                    logFdr
                )
                if (mergedSensitivity != minAdditionalSensitivity) {
                    LOG.info(
                        "${name?.let { "$it " } ?: ""}Merged sensitivity (${"%.3f".format(minAdditionalSensitivity)} + ${"%.3f".format(logFdr)}): ${"%.3f".format(mergedSensitivity)}"
                    )
                } else {
                    LOG.info("${name?.let { "$it " } ?: ""}Using estimated sensitivity")
                }
                mergedSensitivity
            }

            else -> {
                LOG.info("${name?.let { "$it " } ?: ""}Using FDR sensitivity ${"%.3f".format(logFdr)}")
                logFdr
            }
        }

        // Estimate gap
        val (estimatedGap, candidateGapNs) = if (summits && gapCmdArg == null) {
            0 to intArrayOf()
        } else {
            estimateGap(
                genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                if (minAdditionalSensitivity != null) minAdditionalSensitivity else sensitivity,
                parallel = parallel,
                avgAutocorrelation = avgAutocorrelation
            )
        }
        val gap = when {
            gapCmdArg != null -> {
                LOG.info("${name?.let { "$it " } ?: ""}Using command line gap")
                gapCmdArg
            }

            summits -> {
                LOG.info("${name?.let { "$it " } ?: ""}Using no gap for summits")
                0
            }

            else -> {
                LOG.info("${name?.let { "$it " } ?: ""}Using estimated gap")
                estimatedGap
            }
        }

        LOG.info("${name?.let { "$it " } ?: ""}Sensitivity: ${"%.3f".format(sensitivity)}")
        LOG.info("${name?.let { "$it " } ?: ""}Gap: $gap")
        return SensitivityGapResult(
            sensitivity = sensitivity,
            gap = gap,
            avgAutocorrelation = avgAutocorrelation,
            sensitivityInfo = sensitivityInfo,
            candidateGapNs = candidateGapNs
        )
    }

    fun estimateSensitivity(
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        parallel: Boolean,
        name: String?,
        cancellableState: CancellableState?
    ): SensitivityInfo {
        cancellableState?.checkCanceled()

        LOG.debug("${name ?: ""} Analyzing candidates by sensitivities...")
        val (sensitivities, candidatesNs, candidatesALs) = getSensitivitiesAndCandidatesCharacteristics(
            genomeQuery,
            omnipeakFitResults,
            logNullMembershipsMap,
            bitList2reuseMap
        )
        val st = detectSensitivityTriangle(sensitivities, candidatesNs, candidatesALs)
        cancellableState?.checkCanceled()
        if (st != null) {
            LOG.debug("${name ?: ""} Analyzing candidates additive numbers...")
            val sensitivitiesLimited =
                sensitivities.slice(st.beforeMerge until st.stable).toDoubleArray()
            val (totals, news) = analyzeAdditiveCandidates(
                genomeQuery,
                omnipeakFitResults.fitInfo,
                logNullMembershipsMap,
                bitList2reuseMap,
                sensitivitiesLimited,
                parallel
            )
            val newCandidatesList = sensitivitiesLimited.indices
                .map { news[it].toDouble() / totals[it].toDouble() }
            val minAdditionalIdx = newCandidatesList.indices.minByOrNull { newCandidatesList[it] }!!
            val minAdditionalSensitivity = sensitivitiesLimited[minAdditionalIdx]
            LOG.debug(
                "${name ?: ""} " +
                        "Minimal additional ${st.beforeMerge + minAdditionalIdx}: $minAdditionalSensitivity"
            )
            return SensitivityInfo(
                sensitivities, candidatesNs, candidatesALs,
                st, totals, news, minAdditionalIdx, minAdditionalSensitivity
            )
        } else {
            LOG.error("${name ?: ""} Failed to estimate sensitivity, using defaults.")
            return SensitivityInfo(sensitivities, candidatesNs, candidatesALs)
        }
    }


    /**
     * Interpolates sensitivity between [stricterSensitivity] and [stableSensitivity] using [avgAutocorrelation]
     * (broadness index).
     *
     * The track's broadness index ([avgAutocorrelation]) governs the interpolation:
     * - Narrow marks (low autocorrelation, ~0.0) shift sensitivity toward [stricterSensitivity] (e.g. logFdr).
     * - Broad marks (high autocorrelation, ~1.0) keep more of original [stableSensitivity] (minAdditionalSensitivity).
     */
    fun interpolateSensitivity(
        avgAutocorrelation: Double,
        stricterSensitivity: Double,
        stableSensitivity: Double
    ): Double {
        val broadness = avgAutocorrelation.coerceIn(0.0, 1.0)
        return stricterSensitivity + (stableSensitivity - stricterSensitivity) * broadness
    }

    /**
     * Computes merged main sensitivity combining stable sensitivity ([minAdditionalSensitivity]) and [logFdr]
     * using [avgAutocorrelation] interpolation.
     *
     * When [logFdr] is stricter than stable sensitivity, narrow marks shift toward [logFdr],
     * whereas broad domains remain close to stable sensitivity.
     */
    fun computeMergedSensitivity(
        avgAutocorrelation: Double,
        minAdditionalSensitivity: Double?,
        logFdr: Double
    ): Double {
        if (minAdditionalSensitivity == null) {
            return logFdr
        }
        if (logFdr >= minAdditionalSensitivity) {
            return minAdditionalSensitivity
        }
        return interpolateSensitivity(avgAutocorrelation, logFdr, minAdditionalSensitivity)
    }

    fun linSpace(min: Double, max: Double, n: Int): DoubleArray {
        // If we use linear space here, we can't see the plots with merge
        // return DoubleArray(n) {
        //     min + (max - min) * it.toDouble() / (n - 1)
        // }
        require(min * max >= 0) { "Both min and max should have same sign, got $min, $max" }
        val sign = if (min + max < 0) -1 else 1
        val maxLog = log10(max * sign)
        val minLog = log10(min * sign)
        return DoubleArray(n) {
            sign * exp((minLog + (maxLog - minLog) * it.toDouble() / (n - 1)) * ln(10.0))
        }
    }

    private fun triangleSignedSquare(
        x1: Double, y1: Double,
        x2: Double, y2: Double,
        x3: Double, y3: Double
    ) =
        x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3

    /**
     * Detects major semantic changes along sensitivity values.
     * 1) Before merge adjacent candidates
     * 2) Merge - stable
     * 3) Before noise calling
     * See [PepInfo] for details
     */
    fun detectSensitivityTriangle(
        sensitivities: DoubleArray,
        candidatesNs: IntArray,
        candidatesALs: DoubleArray
    ): PepInfo? {
        LOG.debug("Compute sensitivity triangle...")
        require(sensitivities.size == candidatesNs.size)
        require(sensitivities.size == candidatesALs.size)
        val logNs = DoubleArray(candidatesNs.size) { ln1p(candidatesNs[it].toDouble()) }
        val logALs = DoubleArray(candidatesALs.size) { ln1p(candidatesALs[it]) }
        val n = sensitivities.size
        var maxArea = 0.0
        var i1 = -1
        var i2 = -1
        var i3 = -1
        val im1 = (n * 0.2).toInt()
        val im2 = sensitivities.indices.maxBy { logNs[it] } - 1
        for (i in im1..im2) {
            val i1mab = findSensitivityTriangleMaxAreaBetween(
                logNs, logALs, 0, i, -1
            )
            val i3mab = findSensitivityTriangleMaxAreaBetween(
                logNs, logALs, i, n - 1, -1
            )
            if (i1mab.first == -1 || i3mab.first == -1) {
                continue
            }
            // We want both parts to be balanced so geometric mean optimization is better here
            val area = sqrt(i1mab.second * i3mab.second)
            if (area > maxArea) {
                maxArea = area
                i1 = i1mab.first
                i2 = i
                i3 = i3mab.first
            }
        }
        if (i1 == -1 || i2 == -1 || i3 == -1) {
            return null
        }
        // Update i3, i1 points to be closer to i2 for more accurate pivot estimations
        val i3mab = findSensitivityTriangleMaxAreaBetween(logNs, logALs, i3, i2, -1)
        if (i3mab.first != -1) {
            i3 = i3mab.first
        }
        val i1mab = findSensitivityTriangleMaxAreaBetween(logNs, logALs, i2, i1, -1)
        if (i1mab.first != -1) {
            i1 = i1mab.first
        }
        val result = PepInfo(i1, i2, i3)
        LOG.debug(
            "Result beforeMerge: {}: {}, stable: {}: {}, beforeNoise: {}: {}",
            i1, sensitivities[i1], i2, sensitivities[i2], i3, sensitivities[i3],
        )
        return result
    }

    private fun findSensitivityTriangleMaxAreaBetween(
        candidatesNs: DoubleArray,
        candidatesALs: DoubleArray,
        start: Int, end: Int, sign: Int = 1
    ): Pair<Int, Double> {
        if (start > end) {
            return findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, end, start, sign)
        }
        var maxI = -1
        var maxArea = 0.0

        val startN = candidatesNs[start]
        val startAL = candidatesALs[start]
        val endN = candidatesNs[end]
        val endAL = candidatesALs[end]

        for (i in start + 1 until end) {
            val n = candidatesNs[i]
            val al = candidatesALs[i]
            var area = triangleSignedSquare(startN, startAL, n, al, endN, endAL)
            if (area * sign > 0)
                continue
            area = abs(area)
            if (area > maxArea) {
                maxI = i
                maxArea = area
            }
        }
        return maxI to maxArea
    }


    /**
     * Estimates the gap to merge fragmented candidates for a given sensitivity.
     *
     * Gap merging compensates for local fragmentation in broad chromatin domains, where dips in
     * coverage split a single continuous region into multiple pieces.
     *
     * As gap increases from 0, the number of candidate intervals N(gap) declines because adjacent
     * candidates join. The *shape* of that decline in log space tells narrow from broad marks:
     *
     * - Narrow marks (e.g. H3K27ac) merge at a roughly constant relative rate, so `ln(N)` falls
     *   almost linearly - a straight line on the log scale. There is little fragmentation to undo,
     *   so the gap stays close to [OMNIPEAK_DEFAULT_GAP].
     * - Broad marks (e.g. H3K36me3) fragment heavily: many candidates merge already at small gaps
     *   and the merge rate decelerates, so `ln(N)` bows below the straight chord connecting its
     *   ends. The more pronounced that convexity, the larger the gap.
     *
     * Combines empirical fragmentation curvature (the area between the chord and the `ln(N)` curve
     * normalized by the total log decline) with the spatial autocorrelation correlation length scale:
     *   correlationLength = -1 / ln(broadness)
     *   gap = empiricalGap * broadness + correlationLength * broadness^2
     *
     * @param candidatesNs Numbers of candidates by gap to use for the gap estimation.
     * @param avgAutocorrelation Average autocorrelation of posterior error probabilities in [0.0, 1.0].
     * @return An integer representing the gap.
     */
    fun estimateGap(
        candidatesNs: IntArray,
        avgAutocorrelation: Double = 0.0
    ): Int {
        val n = candidatesNs.size
        if (n < 2) return OMNIPEAK_DEFAULT_GAP
        // Laplace-smoothed (+1) log counts, robust for small candidate numbers
        val logNs = DoubleArray(n) { ln1p(candidatesNs[it].toDouble()) }
        val totalDecline = logNs.first() - logNs.last()
        // Curve doesn't decline (no merging at all) - nothing to compensate.
        if (totalDecline <= 0.0) return OMNIPEAK_DEFAULT_GAP
        val chordSum = (logNs.first() + logNs.last()) * n / 2.0
        val curveSum = logNs.sum()
        val area = chordSum - curveSum
        if (area <= 0.0) return OMNIPEAK_DEFAULT_GAP
        val broadness = avgAutocorrelation.coerceIn(0.0, 1.0)
        val correlationLength = if (broadness > 0.0 && broadness < 1.0) {
            min(n.toDouble(), -1.0 / ln(broadness))
        } else if (broadness >= 1.0) {
            n.toDouble()
        } else {
            0.0
        }
        val empiricalGap = area / totalDecline
        val gap = (empiricalGap * broadness + correlationLength * broadness.pow(2)).roundToInt()
        return gap.coerceIn(OMNIPEAK_DEFAULT_GAP, n - 1)
    }

    /**
     * Estimates candidate numbers for gaps 0..[maxGap] in a single pass over chromosomes at gap=0.
     *
     * Adjacent candidate intervals [s_k, e_k] and [s_{k+1}, e_{k+1}] merge under aggregate(gap)
     * iff (s_{k+1} - e_k) <= gap. Each merged adjacent pair reduces the total candidate count by 1.
     * Therefore, collecting the histogram of adjacent gap distances allows computing N(g) for all
     * gap values directly in O(total candidates) instead of scanning the full genome for every gap.
     */
    fun estimateCandidateGapNs(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        maxGap: Int = OMNIPEAK_FRAGMENTATION_MAX_GAP_BP / fitInfo.binSize,
        parallel: Boolean = true
    ): IntArray {
        val totalDistanceCounts = IntArray(maxGap + 1)
        val totalCandidates0 = java.util.concurrent.atomic.AtomicInteger(0)

        genomeQuery.get().map { chromosome ->
            Callable {
                if (!fitInfo.containsChromosomeInfo(chromosome)) {
                    return@Callable
                }
                val logNullMemberships = logNullMembershipsMap[chromosome]
                val bitList2reuse = bitList2reuseMap[chromosome]
                val candidates = getChromosomeCandidates(
                    chromosome, logNullMemberships, bitList2reuse, sensitivity, 0
                )
                if (candidates.isEmpty()) {
                    return@Callable
                }
                totalCandidates0.addAndGet(candidates.size)
                val localDistanceCounts = IntArray(maxGap + 1)
                for (i in 0 until candidates.size - 1) {
                    val dist = candidates[i + 1].startOffset - candidates[i].endOffset
                    if (dist in 1..maxGap) {
                        localDistanceCounts[dist]++
                    }
                }
                synchronized(totalDistanceCounts) {
                    for (i in 1..maxGap) {
                        totalDistanceCounts[i] += localDistanceCounts[i]
                    }
                }
            }
        }.await(parallel)

        val candidateGapNs = IntArray(maxGap + 1)
        candidateGapNs[0] = totalCandidates0.get()
        for (g in 1..maxGap) {
            candidateGapNs[g] = candidateGapNs[g - 1] - totalDistanceCounts[g]
        }
        return candidateGapNs
    }

    /**
     * Estimates the gap to merge fragmented candidates for a given sensitivity,
     * returning both the optimal gap and the candidate counts array by gap.
     */
    fun estimateGap(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        maxGap: Int = OMNIPEAK_FRAGMENTATION_MAX_GAP_BP / fitInfo.binSize,
        parallel: Boolean = true,
        avgAutocorrelation: Double = 0.0
    ): Pair<Int, IntArray> {
        val candidateGapNs = estimateCandidateGapNs(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap, sensitivity, maxGap, parallel
        )
        val gap = estimateGap(candidateGapNs, avgAutocorrelation)
        return Pair(gap, candidateGapNs)
    }


    internal fun getSensitivitiesAndCandidatesCharacteristics(
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>
    ): Triple<DoubleArray, IntArray, DoubleArray> {
        val minLogNull = genomeQuery.get().minOf { logNullMembershipsMap[it].min() }
        // Limit value due to floating point errors
        val maxLogNull = min(OMNIPEAK_MIN_SENSITIVITY, genomeQuery.get().maxOf { logNullMembershipsMap[it].max() })
        return getSensitivitiesAndCandidatesCharacteristics(
            minLogNull,
            maxLogNull,
            genomeQuery,
            omnipeakFitResults,
            logNullMembershipsMap,
            bitList2reuseMap
        )
    }

    private fun getSensitivitiesAndCandidatesCharacteristics(
        minLogNull: Double,
        maxLogNull: Double,
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>
    ): Triple<DoubleArray, IntArray, DoubleArray> {
        val sensitivities = linSpace(minLogNull, maxLogNull, OMNIPEAK_SENSITIVITY_N)
        // Compute candidates characteristics
        val (candidatesNs, candidatesALs) = candidatesNumbersLengths(
            sensitivities,
            genomeQuery,
            omnipeakFitResults,
            logNullMembershipsMap,
            bitList2reuseMap
        )
        val equalTail = candidatesNs.indices.reversed().takeWhile {
            candidatesNs[it] == candidatesNs.last()
        }.count()
        if (equalTail <= 5) {
            if (LOG.isDebugEnabled) {
                LOG.debug("Sensitivity table")
                println("Sensitivity\tCandidatesN\tCandidatesAL")
                for (i in sensitivities.indices) {
                    println("${sensitivities[i]}\t${candidatesNs[i]}\t${candidatesALs[i]}")
                }
            }
            return Triple(sensitivities, candidatesNs, candidatesALs)
        }
        return getSensitivitiesAndCandidatesCharacteristics(
            minLogNull,
            sensitivities[sensitivities.size - equalTail + 1],
            genomeQuery,
            omnipeakFitResults,
            logNullMembershipsMap,
            bitList2reuseMap
        )
    }

    private fun candidatesNumbersLengths(
        sensitivities: DoubleArray,
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>
    ): Pair<IntArray, DoubleArray> {
        val n = sensitivities.size
        val candidatesNs = IntArray(n)
        val candidatesALs = DoubleArray(n)
        for ((i, s) in sensitivities.withIndex()) {
            val ci = estimateCandidatesNumberLens(
                genomeQuery, omnipeakFitResults.fitInfo, logNullMembershipsMap, bitList2reuseMap,
                s, 0
            )
            candidatesNs[i] = ci.n
            candidatesALs[i] = ci.averageLen
        }
        return Pair(candidatesNs, candidatesALs)
    }

    fun estimateCandidatesNumberLens(
        genomeQuery: GenomeQuery,
        omnipeakFitInformation: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        gap: Int
    ): CandidatesInfo {
        val lens = TDoubleArrayList()
        genomeQuery.get().map { chromosome ->
            Callable {
                if (!omnipeakFitInformation.containsChromosomeInfo(chromosome)) {
                    return@Callable
                }
                val logNullMemberships = logNullMembershipsMap[chromosome]
                val bitList2reuse = bitList2reuseMap[chromosome]
                val candidates = getChromosomeCandidates(
                    chromosome, logNullMemberships, bitList2reuse, sensitivity, gap,
                )
                if (candidates.isNotEmpty()) {
                    synchronized(lens) {
                        lens.add(DoubleArray(candidates.size) { candidates[it].length().toDouble() })
                    }
                }
            }
        }.await(true)
        return CandidatesInfo(
            lens.size(),
            if (lens.size() == 0) 0.0 else lens.sum() / lens.size(),
            if (lens.size() == 0) 0.0 else StatUtils.percentile(lens.toArray(), 50.0),
            if (lens.size() == 0) 0 else lens.max().toInt(),
        )
    }


    fun analyzeAdditiveCandidates(
        genomeQuery: GenomeQuery,
        omnipeakFitInformation: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivities: DoubleArray,
        parallel: Boolean
    ): Pair<IntArray, IntArray> {
        // Collect candidates from model for most strict sensitivity
        // We assume that sensitivities are sorted ascending
        var candidatesPrev = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            if (!omnipeakFitInformation.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>().toRangeMergingList()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(
                chromosome, logNullMemberships, bitList2reuse, sensitivities[0], 0
            ).toRangeMergingList()
        }

        val totals = IntArray(sensitivities.size)
        val news = IntArray(sensitivities.size)
        val candidatesPrevSize = genomeQuery.get().sumOf { candidatesPrev[it].size }
        totals[0] = candidatesPrevSize
        news[0] = candidatesPrevSize
        for (i in 1 until sensitivities.size) {
            val s = sensitivities[i]
            val candidates = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
                if (!omnipeakFitInformation.containsChromosomeInfo(chromosome)) {
                    return@genomeMap emptyList<Range>().toRangeMergingList()
                }
                val logNullMemberships = logNullMembershipsMap[chromosome]
                val bitList2reuse = bitList2reuseMap[chromosome]
                getChromosomeCandidates(
                    chromosome, logNullMemberships, bitList2reuse, s, 0,
                ).toRangeMergingList()
            }
            val total = genomeQuery.get().sumOf { candidates[it].size }
            val old = genomeQuery.get().sumOf { chromosome ->
                val chrCandidatesPrev = candidatesPrev[chromosome]
                candidates[chromosome].count { chrCandidatesPrev.intersectionLength(it) > 0 }
            }
            val new = total - old
            totals[i] = total
            news[i] = new
            candidatesPrev = candidates
        }
        return totals to news
    }

}
