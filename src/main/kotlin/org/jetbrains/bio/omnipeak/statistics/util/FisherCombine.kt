package org.jetbrains.bio.omnipeak.statistics.util

import org.jetbrains.bio.omnipeak.statistics.util.PoissonUtil.PREFIX_LOGS_SUMS
import org.jetbrains.bio.omnipeak.statistics.util.PoissonUtil.lnTab
import org.jetbrains.bio.viktor.logAddExp
import kotlin.collections.sumOf
import kotlin.math.ln

object FisherCombine {
    private fun lnFactorial(i: Int): Double =
        if (i < PREFIX_LOGS_SUMS.size)
            PREFIX_LOGS_SUMS[i]
        else
            PREFIX_LOGS_SUMS.last() + (PREFIX_LOGS_SUMS.size..i).sumOf { lnTab(it) }

    /**
     * Natural log of the regularized upper incomplete gamma function `Q(k, x)`
     * for a positive integer shape `k`, evaluated in log space for stability.
     *
     * For integer `k`, `Q(k, x) = e^{-x} * sum_{i=0}^{k-1} x^i / i!`.
     */
    private fun logRegularizedUpperGammaInt(k: Int, x: Double): Double {
        require(k >= 1) { "Shape should be >= 1, got $k" }
        if (x <= 0.0) return 0.0                          // Q(k, 0) = 1 -> log(1) = 0
        if (x.isInfinite()) return Double.NEGATIVE_INFINITY
        val lnX = ln(x)
        var logSum = Double.NEGATIVE_INFINITY
        for (i in 0 until k) {
            // ln(x^i / i!) = i * ln(x) - ln(i!)
            logSum = logSum logAddExp (i * lnX - lnFactorial(i))
        }
        // Cap at log(1) = 0 to guard against floating point overshoot.
        return minOf(0.0, logSum - x)
    }

    /**
     * Combines independent p-values using Fisher's method, in log space.
     *
     * Fisher's statistic `S = -2 * sum(ln p_i)` follows a chi-squared distribution
     * with `2k` degrees of freedom under the null, so the combined p-value is
     * `P(chi^2_{2k} >= S) = Q(k, S/2)`, where `Q` is the regularized upper
     * incomplete gamma function. As `k` (the number of combined p-values) is an
     * integer, `Q(k, .)` has a closed form that is evaluated entirely in log space
     * to remain accurate for the extremely small p-values produced at genome scale.
     *
     * Combining a single p-value returns it unchanged.
     *
     * NOTE: Fisher's method assumes the combined p-values are independent.
     *
     * @param logPs natural-log p-values, each expected to be <= 0.
     * @return natural log of the combined p-value.
     */
    fun logFisherCombinedP(logPs: DoubleArray): Double {
        require(logPs.isNotEmpty()) { "No p-values to combine" }
        // x = S / 2 = -sum(ln p_i); clamp per-term at log(1) = 0 against fp noise.
        var x = 0.0
        for (logP in logPs) {
            x -= minOf(0.0, logP)
        }
        return logRegularizedUpperGammaInt(logPs.size, x)
    }
}