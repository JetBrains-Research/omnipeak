package org.jetbrains.bio.omnipeak.statistics.util

import org.apache.commons.math3.special.Beta
import org.apache.commons.math3.special.Gamma
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_HMM_ESTIMATE_LOW
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_HMM_MAX_MEAN_TO_STD
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_HMM_NB_VAR_MEAN_MULTIPLIER
import org.jetbrains.bio.statistics.average
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution.Companion.estimateFailuresUsingMoments
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.viktor.logAddExp
import org.slf4j.LoggerFactory
import kotlin.math.abs
import kotlin.math.ln
import kotlin.math.max
import kotlin.math.pow

object NegBinUtil {

    private val LOG = LoggerFactory.getLogger(NegBinUtil::class.java)

    @Suppress("ArrayInDataClass")
    data class Guess(
        val means: DoubleArray,
        val failures: DoubleArray,
        val lowMin: Double,
        val signalToNoise: Double
    )

    /**
     * Guess initial approximation for mixture / HMM given data for several Negative Binomial distributions.
     */
    fun guessByData(
        emissions: IntArray,
        n: Int,
        estimateSNRFraction: Double = OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR,
        estimateLowFraction: Double = OMNIPEAK_HMM_ESTIMATE_LOW,
        estimateLowMinThreshold: Double = OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD,
        maxMeanToStd: Double = OMNIPEAK_HMM_MAX_MEAN_TO_STD
    ): Guess {
        require(estimateSNRFraction + estimateLowFraction < 1.0) {
            "Estimate SNR fraction is too high $estimateSNRFraction"
        }
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        val vars = sd * sd
        LOG.debug("All emissions mean $mean\t std $sd")
        emissions.sortDescending()

        val fraction = if (estimateSNRFraction > 0) estimateSNRFraction else OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR
        val highCount = (emissions.size * fraction).toInt()
        var meanH = emissions.average(0, highCount)
        val sdH = emissions.standardDeviation(0, highCount)
        LOG.debug("High $fraction emissions mean $meanH\t std $sdH")
        if (meanH > maxMeanToStd * (sdH + 1e-10)) {
            LOG.warn("High mean / std > $maxMeanToStd, adjusting...")
            meanH = (sdH + 1e-10) * maxMeanToStd
            LOG.debug("Adjusted high $fraction emissions mean $meanH\t std $sdH")
        }

        // Most likely close to 1, but in case when all the coverage is extremely high this helps
        val lowCount = (emissions.size * estimateLowFraction).toInt()
        val lowOffset = emissions.size - lowCount
        val meanL = emissions.average(lowOffset, lowCount)
        val sdL = emissions.standardDeviation(lowOffset, lowCount)
        LOG.debug("Low $estimateLowFraction emissions mean $meanL\t std $sdL")

        val signalToNoiseRatio = meanH / meanL

        // Make means proportional to snr, minimal should be equal to the noise state
        val means = DoubleArray(n) { meanL * signalToNoiseRatio.pow(it.toDouble() / (n - 1)) }
        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val failures = DoubleArray(n) {
            estimateFailuresUsingMoments(means[it], max(means[it] * OMNIPEAK_HMM_NB_VAR_MEAN_MULTIPLIER, vars))
        }
        LOG.debug(
            "Guess NegBinEmissionScheme means ${means.joinToString(",")}, " +
                    "failures ${failures.joinToString(",")}"
        )
        val snrMin = if (estimateSNRFraction > 0) signalToNoiseRatio else 0.0
        LOG.debug("Minimal signal to noise ratio: $snrMin")
        val lowMin = meanL * estimateLowMinThreshold
        LOG.debug("Minimal low state: $lowMin")

        return Guess(
            means = means,
            failures = failures,
            lowMin = lowMin,
            signalToNoise = snrMin
        )
    }

    /**
     * Negative Binomial CDF evaluator for upper tail P(X >= k) in log space.
     * @param k observation count
     * @param mu expected mean background
     * @param r number of failures (overdispersion parameter)
     * @return log(p-value) = ln P(X >= k)
     */
    fun logNegBinCdf(
        k: Int,
        mu: Double,
        r: Double,
        maxM: Int = 10_000,
        epsilon: Double = 1e-5
    ): Double {
        if (k <= 0) return 0.0
        if (mu <= 0.0) return Double.NEGATIVE_INFINITY
        if (r.isInfinite() || r <= 0.0 || r > 1e6) {
            return PoissonUtil.logPoissonCdf(k - 1, mu, maxM, epsilon)
        }
        val p = mu / (mu + r)
        if (p <= 0.0) return Double.NEGATIVE_INFINITY
        if (p >= 1.0) return 0.0

        val betaProb = try {
            Beta.regularizedBeta(p, k.toDouble(), r)
        } catch (_: Exception) {
            0.0
        }
        if (betaProb > 1e-15) {
            return ln(betaProb)
        }

        val lnP = ln(p)
        val ln1MinusP = ln(r / (mu + r))

        // Initial log probability for X = k
        var logX = Gamma.logGamma(k + r) - PoissonUtil.logFactorial(k) - Gamma.logGamma(r) +
                r * ln1MinusP + k * lnP
        var residue = logX

        for (i in k + 1..k + maxM) {
            val logY = logX + ln(i - 1 + r) - PoissonUtil.lnTab(i) + lnP
            val preResidue = residue
            residue = preResidue logAddExp logY
            if (abs(preResidue - residue) < epsilon) {
                break
            }
            logX = logY
        }
        return residue
    }
}