package org.jetbrains.bio.omnipeak.peaks

import com.google.common.util.concurrent.AtomicDouble
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_CLIP_MAX_LENGTH
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_CLIP_STEPS
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.LOG
import org.jetbrains.bio.util.await
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong
import kotlin.math.min

object Signal {

    fun estimateGenomeSignalNoiseAverage(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        candidates: GenomeMap<List<Range>>,
        coverageComputable: (ChromosomeRange) -> Double,
        parallel: Boolean
    ): Pair<Double, Double> {
        val sumSignalScoreA = AtomicDouble()
        val sumSignalLengthA = AtomicLong()
        val sumNoiseScoreA = AtomicDouble()
        val sumNoiseLengthA = AtomicLong()
        genomeQuery.get().map { chromosome ->
            Callable {
                if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                    return@Callable
                }
                val chrCandidates = candidates[chromosome]
                val offsets = fitInfo.offsets(chromosome)
                var prevNoiseStart = 0
                chrCandidates.forEach { (from, to) ->
                    val start = offsets[from]
                    val end = if (to < offsets.size) offsets[to] else chromosome.length
                    val range = ChromosomeRange(start, end, chromosome)
                    val score = coverageComputable(range)
                    sumSignalScoreA.addAndGet(score)
                    sumSignalLengthA.addAndGet(end.toLong() - start)
                    val rangeNoise = ChromosomeRange(prevNoiseStart, start, chromosome)
                    val scoreNoise = coverageComputable(rangeNoise)
                    sumNoiseScoreA.addAndGet(scoreNoise)
                    sumNoiseLengthA.addAndGet(start.toLong() - prevNoiseStart)
                    prevNoiseStart = end
                }
                if (prevNoiseStart < chromosome.length) {
                    val rangeNoise = ChromosomeRange(prevNoiseStart, chromosome.length, chromosome)
                    val scoreNoise = coverageComputable(rangeNoise)
                    sumNoiseScoreA.addAndGet(scoreNoise)
                    sumNoiseLengthA.addAndGet(chromosome.length.toLong() - prevNoiseStart)
                }
            }
        }.await(parallel)
        val sumSignalScore = sumSignalScoreA.get()
        val sumSignalLength = sumSignalLengthA.get()
        val sumNoiseScore = sumNoiseScoreA.get()
        val sumNoiseLength = sumNoiseLengthA.get()
        val avgSignalDensity = if (sumSignalLength > 0) sumSignalScore / sumSignalLength else 0.0
        val avgNoiseDensity = if (sumNoiseLength > 0) sumNoiseScore / sumNoiseLength else 0.0
        if (sumSignalLength != 0L && sumNoiseLength != 0L && avgSignalDensity <= avgNoiseDensity) {
            LOG.warn("Average signal density $avgSignalDensity <= average noise density $avgNoiseDensity")
        }
        return avgSignalDensity to avgNoiseDensity
    }

    fun computeSignalToControlAverage(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        candidates: GenomeMap<List<Range>>,
        coverageControlComputable: (ChromosomeRange) -> Pair<Double, Double>,
        parallel: Boolean
    ): Double {
        val signalToControls = AtomicDouble()
        val signalToControlsN = AtomicInteger()
        genomeQuery.get().map { chromosome ->
            Callable {
                if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                    return@Callable
                }
                val offsets = fitInfo.offsets(chromosome)
                candidates[chromosome].forEach { (from, to) ->
                    val start = offsets[from]
                    val end = if (to < offsets.size) offsets[to] else chromosome.length
                    val chromosomeRange = ChromosomeRange(start, end, chromosome)
                    val (score, controlScore) = coverageControlComputable(chromosomeRange)
                    val ratio = if (controlScore != 0.0) score / controlScore else 0.0
                    signalToControls.addAndGet(ratio)
                    signalToControlsN.addAndGet(1)
                }
            }
        }.await(parallel)
        return if (signalToControls.get() > 0) signalToControls.get() / signalToControls.get() else 0.0
    }

    fun clipPeakBySignal(
        chromosome: Chromosome,
        start: Int,
        end: Int,
        fitInfo: OmnipeakFitInformation,
        avgSignalDensity: Double,
        avgNoiseDensity: Double,
        clipSignal: Double = OMNIPEAK_DEFAULT_CLIP_MAX_SIGNAL,
        clipLength: Double = OMNIPEAK_CLIP_MAX_LENGTH,
        clipSteps: DoubleArray = OMNIPEAK_CLIP_STEPS,
    ): Pair<Int, Int> {
        if (avgSignalDensity <= avgNoiseDensity) {
            return start to end
        }
        // Additionally, clip peaks by local coverage signal
        val maxClippedDensity = avgNoiseDensity + clipSignal * (avgSignalDensity - avgNoiseDensity)
        val maxClippedSideLength = (end - start) * clipLength / 2
        val bin = fitInfo.binSize

        // Try to change the left boundary
        val maxStart = start + maxClippedSideLength
        var clippedStart = start
        var step = clipSteps.size - 1
        while (step >= 0 && clippedStart <= maxStart) {
            val newStart = clippedStart + (clipSteps[step] * bin).toInt()
            if (newStart > maxStart) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clippedDensity = fitInfo.score(ChromosomeRange(start, newStart, chromosome)) / (newStart - start)
            if (clippedDensity < maxClippedDensity) {
                clippedStart = newStart
                step = min(step + 1, clipSteps.size - 1)
            } else {
                step -= 1
            }
        }
        // Try to change the right boundary
        val minEnd = end - maxClippedSideLength
        var clippedEnd = end
        step = clipSteps.size - 1
        while (step >= 0 && clippedEnd >= minEnd) {
            val newEnd = clippedEnd - (clipSteps[step] * bin).toInt()
            if (newEnd < minEnd) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clippedDensity = fitInfo.score(ChromosomeRange(newEnd, end, chromosome)) / (end - newEnd)
            if (clippedDensity < maxClippedDensity) {
                clippedEnd = newEnd
                step = min(step + 1, clipSteps.size - 1)
            } else {
                step -= 1
            }
        }
        return clippedStart to clippedEnd
    }

}