package org.jetbrains.bio.omnipeak.statistics.util

import org.jetbrains.bio.viktor.logAddExp
import kotlin.math.abs
import kotlin.math.ln

object PoissonUtil {
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

    // Tabulating for ln(x) and prefix sums of ln(x)
    private const val MAX_TAB = 10001

    private val LOGS = DoubleArray(MAX_TAB) { ln(it.toDouble()) }

    private fun lnTab(i: Int): Double = if (i < LOGS.size) LOGS[i] else ln(i.toDouble())

    private val PREFIX_SUMS = DoubleArray(MAX_TAB).apply {
        for (i in 1 until size) {
            this[i] = this[i - 1] + LOGS[i]
        }
    }

    fun logPoissonCdf(k: Int, lbd: Double, maxM: Int = 10_000, epsilon: Double = 1e-5): Double {
        require(lbd > 0) {
            "Lambda should be > 0, got $lbd"
        }
        var residue: Double
        var logX: Double
        val lnLbd = ln(lbd)
        // first residue
        val m = k + 1
        val sumLns = if (m < PREFIX_SUMS.size)
            PREFIX_SUMS[m]
        else
            PREFIX_SUMS.last() + (PREFIX_SUMS.size .. m).sumOf { lnTab(it) }
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