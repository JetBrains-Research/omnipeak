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
import kotlin.math.ceil

class PeakScorer(
    val modelScoreComputable: (ChromosomeRange) -> Double,
    val coverageControlComputable: ((ChromosomeRange) -> Pair<Double, Double>)?,
    val coverageComputable: ((ChromosomeRange) -> Double)?,
    val binSize: Int,
    var avgNoiseDensity: Double? = null
) {
    /**
     * Whether [logP] yields genuine p-values rather than model posterior
     * error probabilities (PEPs).
     *
     * A Poisson enrichment term is available (against control or estimated noise
     * background) exactly when a coverage source is present.
     * Without it, the only evidence is the model PEP, which is not a p-value and
     * must be handled by a PEP-based FDR procedure
     * (see [org.jetbrains.bio.statistics.hypothesis.Fdr.qvalidatePEPs]).
     */
    val producesPValues: Boolean
        get() = coverageControlComputable != null || avgNoiseDensity != null

    /**
     * @return logPValue or logPEP
     */
    fun logP(cr: ChromosomeRange): Double {
        return when {
            coverageControlComputable != null -> {
                val (score, controlScore) = coverageControlComputable(cr)
                // Evidence vs. control
                PoissonUtil.logPoissonCdf(ceil(score).toInt() + 1, controlScore + 1)
            }
            avgNoiseDensity != null -> {
                checkNotNull(coverageComputable)
                val score = coverageComputable(cr)
                PoissonUtil.logPoissonCdf(
                    // Evidence vs. background noise
                    ceil(score).toInt() + 1, avgNoiseDensity!! * (cr.endOffset - cr.startOffset) + 1
                )
            }
            else ->
                // Model evidence: mean log PEP over the block's bins
                if (cr.endOffset / binSize > cr.startOffset / binSize)
                    modelScoreComputable(cr)
                else
                   // No evidence available -> p-value 1 (log 0).
                    1.0
        }
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
