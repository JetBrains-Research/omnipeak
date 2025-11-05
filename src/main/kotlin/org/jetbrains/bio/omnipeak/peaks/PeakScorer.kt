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
import kotlin.math.sqrt

class PeakScorer(
    val modelScoreComputable: (ChromosomeRange) -> Double,
    val coverageControlComputable: ((ChromosomeRange) -> Pair<Double, Double>)?,
    val coverageComputable: ((ChromosomeRange) -> Double)?,
    var avgNoiseDensity: Double? = null
) {
    fun logPValue(cr: ChromosomeRange): Double {
        // Model posterior log error probability for block
        val modelLogPs = modelScoreComputable(cr)
        when {
            coverageControlComputable != null -> {
                val (score, controlScore) = coverageControlComputable(cr)
                // Combine both model and signal estimations
                return -sqrt(
                    modelLogPs * PoissonUtil.logPoissonCdf(
                        ceil(score).toInt() + 1, controlScore + 1
                    )
                )
            }

            avgNoiseDensity != null -> {
                checkNotNull(coverageComputable)
                val score = coverageComputable(cr)
                // Combine both model and signal estimations
                return -sqrt(
                    modelLogPs * PoissonUtil.logPoissonCdf(
                        ceil(score).toInt() + 1, avgNoiseDensity!! * (cr.endOffset - cr.startOffset) + 1
                    )
                )
            }

            else ->
                // Fallback
                return modelLogPs

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
                (cr.startOffset / fitInfo.binSize until cr.endOffset / fitInfo.binSize).sumOf {
                    logNullMemberships[it]
                }
            }
            return PeakScorer(
                modelScoreComputable,
                coverageControlComputable,
                coverageComputable,
            )
        }
    }

}
