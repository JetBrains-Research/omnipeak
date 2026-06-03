package org.jetbrains.bio.omnipeak.peaks

import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.peaks.Signal.estimateGenomeSignalNoiseAverage
import org.jetbrains.bio.omnipeak.statistics.util.PoissonUtil
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.pow

class PeakScorer(
    val modelScoreComputable: (ChromosomeRange) -> Double,
    val coverageControlComputable: ((ChromosomeRange) -> Pair<Double, Double>)?,
    val coverageComputable: ((ChromosomeRange) -> Double)?,
    val binSize: Int,
    var avgNoiseDensity: Double? = null
) {
    /**
     * Whether [logP] yields p-values rather than model posterior error probabilities (PEPs).
     */
    val producesPValues: Boolean
        get() = (coverageControlComputable != null || avgNoiseDensity != null)

    /**
     * @return logPValue or logPEP
     */
    fun logP(cr: ChromosomeRange): Double {
        val logPs = ArrayList<Double>()
        if (coverageControlComputable != null) {
            val (score, controlScore) = coverageControlComputable(cr)
            logPs.add(
                // Evidence vs. control
                PoissonUtil.logPoissonCdf(ceil(score).toInt() + 1, controlScore + 1)
            )
        }

        if (avgNoiseDensity != null) {
            checkNotNull(coverageComputable)
            val score = coverageComputable(cr)
            logPs.add(
                PoissonUtil.logPoissonCdf(
                    // Evidence vs. background noise
                    ceil(score).toInt() + 1, avgNoiseDensity!! * (cr.endOffset - cr.startOffset) + 1
                )
            )
        }

        // Model PEP over the block's bins
        logPs.add(modelScoreComputable(cr))

        // Cannot apply FisherCombine here because of dependency between p-values
        // FisherCombine inflates two or more significant pvalues
        return -logPs.fold(1.0) { acc, v -> acc * abs(v) }.pow(1.0 / logPs.size.toDouble())
    }

    fun valueScore(cr: ChromosomeRange): Double {
        when {
            coverageControlComputable != null -> {
                val (score, controlScore) = coverageControlComputable(cr)
                return (ceil(score) + 1) / (controlScore + 1)
            }

            avgNoiseDensity != null -> {
                checkNotNull(coverageComputable)
                val score = coverageComputable(cr)
                return (ceil(score).toInt() + 1) / (avgNoiseDensity!! * (cr.endOffset - cr.startOffset) + 1)
            }

            else -> {
                // Fallback: model posterior log error probability for block
                return -modelScoreComputable(cr)
            }
        }
    }

    fun analyzeSignalAndNoise(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        candidatesMap: GenomeMap<List<Range>>,
        parallel: Boolean
    ): Pair<Double?, Double?> {
        val (s, n) = if (coverageComputable != null)
            estimateGenomeSignalNoiseAverage(
                genomeQuery, fitInfo, candidatesMap, coverageComputable, parallel
            )
        else
            null to null
        avgNoiseDensity = n
        return s to n
    }

    companion object {
        fun create(
            fitInfo: OmnipeakFitInformation,
            logNullMembershipsMap: GenomeMap<F64Array>
        ): PeakScorer {
            fitInfo.prepareData()
            // TODO[oleg] support OmnipeakCompareFitInformation
            val coverageControlComputable =
                if (fitInfo is OmnipeakAnalyzeFitInformation && fitInfo.isControlAvailable())
                    fitInfo.getTreatmentControlComputable()
                else
                    null
            val coverageComputable = if (fitInfo is OmnipeakAnalyzeFitInformation)
                fitInfo.getTreatmentComputable()
            else
                null
            val modelScoreComputable = { cr: ChromosomeRange ->
                val logNullMemberships = logNullMembershipsMap[cr.chromosome]
                val startBin = cr.startOffset / fitInfo.binSize
                val endBin = cr.endOffset / fitInfo.binSize
                // Mean log PEP over the block's bins (length-normalized so that it is
                // comparable across blocks of different sizes); 0.0 for sub-bin blocks.
                if (endBin > startBin)
                    (startBin until endBin).sumOf { logNullMemberships[it] } / (endBin - startBin)
                else
                    0.0
            }
            return PeakScorer(
                modelScoreComputable,
                coverageControlComputable,
                coverageComputable,
                fitInfo.binSize,
            )
        }
    }

}
