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
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FRAGMENTATION_THRESHOLD_BP
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

    fun estimateSensitivity(
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        parallel: Boolean,
        name: String?,
        cancellableState: CancellableState?
    ): Double {
        cancellableState?.checkCanceled()

        OmnipeakModelToPeaks.LOG.info("${name ?: ""} Analyzing candidates by sensitivity...")
        val (sensitivities, candidatesNs, candidatesALs) = getSensitivitiesAndCandidatesCharacteristics(
            genomeQuery,
            omnipeakFitResults,
            logNullMembershipsMap,
            bitList2reuseMap
        )
        val st = detectSensitivityTriangle(sensitivities, candidatesNs, candidatesALs)
        cancellableState?.checkCanceled()
        if (st != null) {
            OmnipeakModelToPeaks.LOG.debug("${name ?: ""} Analyzing candidates additive numbers...")
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
            OmnipeakModelToPeaks.LOG.debug(
                "${name ?: ""} " +
                        "Minimal additional ${st.beforeMerge + minAdditionalIdx}: $minAdditionalSensitivity"
            )
            return minAdditionalSensitivity
        } else {
            OmnipeakModelToPeaks.LOG.error("${name ?: ""} Failed to estimate sensitivity, using defaults.")
            return ln(OMNIPEAK_DEFAULT_FDR)
        }
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
     * Estimates the fragmentation gap based on candidate data and a fragmentation threshold.
     * Fragmentation score is computed as the area above the curve of the plot of
     * candidates number with gap to the candidates number without gap.
     *
     * @param candidatesNs Numbers of candidates by gap to use for the fragmentation score estimation.
     * @param name An optional name or identifier for logging purposes.
     * @param binSize The model bin size to use for the fragmentation score estimation.
     * @param fragmentationThresholdBp The fragmentation threshold value which determines if further processing is needed.
     * @return An integer representing the fragmentation compensation gap.
     * Returns 0 if no fragmentation is detected.
     */
    fun estimateGap(
        candidatesNs: IntArray,
        name: String?,
        binSize: Int,
        fragmentationThresholdBp: Int = OMNIPEAK_DEFAULT_FRAGMENTATION_THRESHOLD_BP,
    ): Int {
        val fragmentations = DoubleArray(candidatesNs.size) {
            candidatesNs[it].toDouble() / candidatesNs[0]
        }
        // Area above the curve
        val fragmentationScore = (fragmentations.size - fragmentations.sum())
        OmnipeakModelToPeaks.LOG.debug(
            "${name ?: ""} Fragmentation $fragmentationScore"
        )
        if (fragmentationScore < fragmentationThresholdBp / binSize) {
            OmnipeakModelToPeaks.LOG.info("${name ?: ""} No fragmentation detected!")
            return 0
        }
        val compensationGap = (fragmentationScore - fragmentationThresholdBp / binSize).toInt()
        OmnipeakModelToPeaks.LOG.info("${name ?: ""} Fragmentation compensation gap: $compensationGap")
        return compensationGap
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
            if (OmnipeakModelToPeaks.LOG.isDebugEnabled) {
                OmnipeakModelToPeaks.LOG.debug("Sensitivity table")
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