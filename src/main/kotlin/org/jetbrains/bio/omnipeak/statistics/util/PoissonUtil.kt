package org.jetbrains.bio.omnipeak.statistics.util

import org.jetbrains.bio.viktor.logAddExp
import kotlin.math.PI
import kotlin.math.abs
import kotlin.math.ln

object PoissonUtil {

    // Tabulating for ln(x) and prefix sums of ln(x)
    private const val MAX_TAB = 10001

    private val LOGS = DoubleArray(MAX_TAB) { ln(it.toDouble()) }

    fun lnTab(i: Int): Double = if (i < LOGS.size) LOGS[i] else ln(i.toDouble())

    internal val PREFIX_LOGS_SUMS = DoubleArray(MAX_TAB).apply {
        for (i in 1 until size) {
            this[i] = this[i - 1] + LOGS[i]
        }
    }

    private val LN_SQRT_PI = 0.5 * ln(PI)

    /**
     * Approximates ln(m!) for large m using Ramanujan's formula for Stirling's approximation:
     *   ln(m!) ≈ m * ln(m) - m + (1/6) * ln(m * (1 + 4m * (1 + 2m))) + (1/2) * ln(π)
     *
     * Reference:
     *   Srinivasa Ramanujan, "The Lost Notebook and Other Unpublished Papers" (1988).
     *   See also: Stirling's approximation (Ramanujan's formula).
     */
    fun logFactorial(m: Int): Double {
        if (m < PREFIX_LOGS_SUMS.size) {
            return PREFIX_LOGS_SUMS[m]
        }
        val md = m.toDouble()
        return md * ln(md) - md + (ln(md * (1.0 + 4.0 * md * (1.0 + 2.0 * md))) / 6.0) + LN_SQRT_PI
    }

    /**
     * Poisson CDF evaluater for upper tail which allow calculation in log space.
     * @param k observation
     * @param lbd: Lambda
     * @return log(pvalue)
     *
     * See MACS2 sources Prob.pyx for original source code (log10PoissonCdfQLargeLambda):
     * ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
     * \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
     * Calculate \ln( sum{exp{N} ) by logspace_add function
     */
    fun logPoissonCdf(k: Int, lbd: Double, maxM: Int = 10_000, epsilon: Double = 1e-5): Double {
        require(lbd > 0) {
            "Lambda should be > 0, got $lbd"
        }
        var residue: Double
        var logX: Double
        val lnLbd = ln(lbd)
        // first residue
        val m = k + 1
        val sumLns = logFactorial(m)
        logX = m * lnLbd - sumLns
        residue = logX
        var logy: Double
        for (i in m + 1..maxM) { // Limit
            logy = logX + lnLbd - lnTab(i)
            val preResidue = residue
            residue = preResidue logAddExp logy
            if (abs(preResidue - residue) < epsilon)
                break
            logX = logy
        }
        return residue - lbd
    }
}