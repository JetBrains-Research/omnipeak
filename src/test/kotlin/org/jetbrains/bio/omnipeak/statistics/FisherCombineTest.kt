package org.jetbrains.bio.omnipeak.statistics

import org.apache.commons.math3.distribution.ChiSquaredDistribution
import org.jetbrains.bio.omnipeak.statistics.util.FisherCombine
import org.junit.Assert.assertEquals
import org.junit.Test
import kotlin.math.exp
import kotlin.math.ln

class FisherCombineTest {

    /** Reference Fisher combination via the chi-squared survival function. */
    private fun fisherReference(ps: DoubleArray): Double {
        val stat = -2.0 * ps.sumOf { ln(it) }
        return 1.0 - ChiSquaredDistribution((2 * ps.size).toDouble()).cumulativeProbability(stat)
    }

    @Test
    fun singlePValueIsUnchanged() {
        for (p in doubleArrayOf(0.5, 0.1, 0.01, 1e-8)) {
            assertEquals(ln(p), FisherCombine.logFisherCombinedP(doubleArrayOf(ln(p))), 1e-9)
        }
    }

    @Test
    fun allOnesCombineToOne() {
        assertEquals(0.0, FisherCombine.logFisherCombinedP(doubleArrayOf(0.0, 0.0, 0.0)), 1e-12)
    }

    @Test
    fun matchesChiSquaredReference() {
        val cases = listOf(
            doubleArrayOf(0.1, 0.1),
            doubleArrayOf(0.05, 0.2, 0.5),
            doubleArrayOf(0.01, 0.4, 0.7, 0.9),
            doubleArrayOf(0.2, 0.3)
        )
        for (ps in cases) {
            val expected = fisherReference(ps)
            val actual = exp(FisherCombine.logFisherCombinedP(DoubleArray(ps.size) { ln(ps[it]) }))
            assertEquals(expected, actual, 1e-6)
        }
    }

    @Test
    fun stableForTinyPValues() {
        // Genome-scale p-values underflow a naive (non-log) Fisher implementation;
        // the log-space combination must stay finite and below log(1) = 0.
        val logPs = doubleArrayOf(-2000.0, -1500.0, -3000.0)
        val combined = FisherCombine.logFisherCombinedP(logPs)
        assert(combined.isFinite()) { "Combined log p-value should be finite, got $combined" }
        assert(combined <= 0.0) { "Combined log p-value should be <= 0, got $combined" }
        // Fisher with k=3: log Q(3, x) ~ -x + 2 ln x for large x; dominated by -x.
        val x = 2000.0 + 1500.0 + 3000.0
        assertEquals(-x + ln(1.0 + x + x * x / 2.0), combined, 1e-6)
    }
}
