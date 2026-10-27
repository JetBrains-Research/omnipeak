package org.jetbrains.bio.omnipeak.peaks

import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_LOCAL_WINDOW
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_MIN_EXPECTED_NOISE
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.peaks.Signal.estimateGenomeSignalNoiseAverage
import org.jetbrains.bio.omnipeak.statistics.util.NegBinUtil
import org.jetbrains.bio.omnipeak.statistics.util.PoissonUtil
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import kotlin.math.ceil
import kotlin.math.max
import kotlin.math.min


/**
 * Scores candidate peak regions using the Hybrid Model: candidate spatial segmentation
 * from HMM combined with empirical read count Negative Binomial / Poisson p-values evaluated
 * against multi-scale local and global background noise. Local background window is dynamically
 * scaled using the track's lag-1 autocorrelation ([avgAutocorrelation]), e.g., 10k * avgAutocorrelation.
 * For datasets without coverage tracks, falls back to model log posterior error
 * probabilities (PEPs).
 */
class CandidatesScores(
    val coverageControlComputable: ((ChromosomeRange) -> Pair<Double, Double>)?,
    val coverageComputable: ((ChromosomeRange) -> Double)?,
    val binSize: Int,
    var avgNoiseDensity: Double? = null,
    val avgAutocorrelation: Double = 0.0,
    val failures: Double? = null,
    val modelScoreComputable: (ChromosomeRange) -> Double
) {
    val localWindow: Int = max(binSize, (OMNIPEAK_DEFAULT_LOCAL_WINDOW * avgAutocorrelation).toInt())
    val halfWindow: Int = localWindow / 2

    /**
     * Whether [logP] yields empirical p-values rather than model posterior error probabilities (PEPs).
     */
    val producesPValues: Boolean
        get() = (coverageControlComputable != null || coverageComputable != null || avgNoiseDensity != null)

    init {
        val scoreType = when {
            !producesPValues -> "posterior error probabilities (PEPs)"
            failures != null && failures > 0.0 && !failures.isInfinite() -> "Negative Binomial p-values"
            else -> "Poisson p-values"
        }
        LOG.info("Scoring candidates using $scoreType...")
    }

    /**
     * @return logPValue or logPEP
     */
    fun logP(cr: ChromosomeRange): Double {
        if (!producesPValues) {
            return modelScoreComputable(cr)
        }

        val length = cr.endOffset - cr.startOffset
        if (length <= 0) return 0.0

        val score: Double
        var muLocal = 0.0
        var muControl = 0.0
        var muControlLocal = 0.0

        if (coverageControlComputable != null) {
            val (treatInCr, controlInCr) = coverageControlComputable.invoke(cr)
            score = treatInCr
            val k = ceil(score).toInt()
            if (k <= 0) return 0.0
            muControl = controlInCr

            val wStart = max(0, cr.startOffset - halfWindow)
            val wEnd = min(cr.chromosome.length, cr.endOffset + halfWindow)
            val wLength = wEnd - wStart
            val flankLength = wLength - length
            if (flankLength > 0) {
                val windowRange = ChromosomeRange(wStart, wEnd, cr.chromosome)
                val (_, wCtrl) = coverageControlComputable.invoke(windowRange)
                val flankCtrl = max(0.0, wCtrl - controlInCr)
                muControlLocal = (flankCtrl / flankLength) * length
            }
        } else if (coverageComputable != null) {
            score = coverageComputable.invoke(cr)
            val k = ceil(score).toInt()
            if (k <= 0) return 0.0

            val wStart = max(0, cr.startOffset - halfWindow)
            val wEnd = min(cr.chromosome.length, cr.endOffset + halfWindow)
            val wLength = wEnd - wStart
            val flankLength = wLength - length
            if (flankLength > 0) {
                val windowRange = ChromosomeRange(wStart, wEnd, cr.chromosome)
                val wTreat = coverageComputable.invoke(windowRange)
                val flankTreat = max(0.0, wTreat - score)
                muLocal = (flankTreat / flankLength) * length
            }
        } else {
            return 0.0
        }

        val k = ceil(score).toInt()
        val muGlobal = if (avgNoiseDensity != null) avgNoiseDensity!! * length else 0.0

        val expectedNoise = if (coverageControlComputable != null) {
            // When control is available, background is determined by control and global noise.
            // Treatment flank coverage is excluded to avoid shadowing broad peaks near sharp spikes.
            maxOf(muGlobal, muControl, muControlLocal)
        } else {
            // When only treatment coverage is available, bound local flank noise relative to global noise
            // to prevent neighboring point-wise peaks from inflating expected noise for broad candidates.
            if (avgNoiseDensity != null && muGlobal > 0.0) {
                maxOf(muGlobal, min(muLocal, muGlobal * 2.0))
            } else {
                maxOf(muGlobal, muLocal)
            }
        }

        val mu = max(OMNIPEAK_MIN_EXPECTED_NOISE, expectedNoise)

        return if (failures != null && failures > 0.0 && !failures.isInfinite()) {
            val numBins = max(1, (length + binSize - 1) / binSize)
            val effBins = effectiveBins(numBins, avgAutocorrelation)
            val scaledFailures = failures * (numBins * numBins) / effBins
            NegBinUtil.logNegBinCdf(k, mu, scaledFailures)
        } else {
            PoissonUtil.logPoissonCdf(k - 1, mu)
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
        private val LOG: Logger = LoggerFactory.getLogger(CandidatesScores::class.java)

        fun create(
            fitInfo: OmnipeakFitInformation,
            logNullMembershipsMap: GenomeMap<F64Array>,
            avgAutocorrelation: Double = 0.0,
            failures: Double? = null
        ): CandidatesScores {
            fitInfo.prepareData()
            // TODO[oleg] support OmnipeakCompareFitInformation
            val coverageControlComputable =
                if (fitInfo is OmnipeakAnalyzeFitInformation && fitInfo.isControlAvailable())
                    fitInfo.getTreatmentControlComputable()
                else
                    null
            val coverageComputable = if (fitInfo is OmnipeakAnalyzeFitInformation && fitInfo.paths.isNotEmpty())
                fitInfo.getTreatmentComputable()
            else
                null
            val modelScoreComputable = { cr: ChromosomeRange ->
                val logNullMemberships = logNullMembershipsMap[cr.chromosome]
                val startBin = cr.startOffset / fitInfo.binSize
                val endBin = cr.endOffset / fitInfo.binSize
                if (endBin > startBin) {
                    val length = endBin - startBin
                    val meanLogPep = (startBin until endBin).sumOf { logNullMemberships[it] } / length
                    val effectiveBins = effectiveBins(length, avgAutocorrelation)
                    meanLogPep * effectiveBins
                } else
                    0.0
            }
            return CandidatesScores(
                coverageControlComputable,
                coverageComputable,
                fitInfo.binSize,
                avgAutocorrelation = avgAutocorrelation,
                failures = failures,
                modelScoreComputable = modelScoreComputable
            )
        }

        /**
         * Computes the effective number of independent informative bins N_eff = 1 + (length - 1) * (1 - rho)
         * for a candidate segment of [length] bins given lag-1 autocorrelation [avgAutocorrelation].
         * For length = 1, returns 1.0 (single locus evidence).
         * For broad correlated domains (rho -> 1.0), returns 1.0 (domain evaluated as a single coherent locus).
         * For uncorrelated tracks (rho -> 0.0), returns length (linear accumulation of independent evidence).
         */
        fun effectiveBins(length: Int, avgAutocorrelation: Double): Double {
            if (length <= 0) return 0.0
            if (length == 1) return 1.0
            val rho = avgAutocorrelation.coerceIn(0.0, 1.0)
            return 1.0 + (length - 1) * (1.0 - rho)
        }

    }

}
