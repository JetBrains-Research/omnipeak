package org.jetbrains.bio.omnipeak.statistics

import org.apache.commons.math3.special.Beta
import org.jetbrains.bio.omnipeak.statistics.util.NegBinUtil
import org.jetbrains.bio.omnipeak.statistics.util.PoissonUtil
import org.junit.Test
import kotlin.math.exp
import kotlin.math.ln
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class NegBinUtilTest {

    @Test
    fun testLogNegBinCdfEdgeCases() {
        // k <= 0 -> probability 1.0 (log = 0.0)
        assertEquals(0.0, NegBinUtil.logNegBinCdf(0, 5.0, 2.0))
        assertEquals(0.0, NegBinUtil.logNegBinCdf(-5, 5.0, 2.0))

        // mu <= 0 -> log = -Infinity
        assertEquals(Double.NEGATIVE_INFINITY, NegBinUtil.logNegBinCdf(5, 0.0, 2.0))
        assertEquals(Double.NEGATIVE_INFINITY, NegBinUtil.logNegBinCdf(5, -1.0, 2.0))

        // r <= 0 or infinite -> delegates to Poisson
        val poissonLogP = PoissonUtil.logPoissonCdf(4, 5.0)
        assertEquals(poissonLogP, NegBinUtil.logNegBinCdf(5, 5.0, Double.POSITIVE_INFINITY), 1e-6)
        assertEquals(poissonLogP, NegBinUtil.logNegBinCdf(5, 5.0, -1.0), 1e-6)
    }

    @Test
    fun testLogNegBinCdfConcordanceWithRegularizedBeta() {
        val testCases = listOf(
            Triple(5, 2.0, 3.0),
            Triple(10, 3.0, 5.0),
            Triple(15, 5.0, 4.0),
            Triple(8, 4.0, 10.0),
            Triple(12, 2.5, 2.0)
        )

        for ((k, mu, r) in testCases) {
            val p = mu / (mu + r)
            val expectedProb = Beta.regularizedBeta(p, k.toDouble(), r)
            val expectedLogProb = ln(expectedProb)
            val actualLogProb = NegBinUtil.logNegBinCdf(k, mu, r)

            assertEquals(expectedLogProb, actualLogProb, 1e-4, "Mismatch for k=$k, mu=$mu, r=$r")
            assertEquals(expectedProb, exp(actualLogProb), 1e-4, "Prob mismatch for k=$k, mu=$mu, r=$r")
        }
    }

    @Test
    fun testLogNegBinCdfMonotonicityInK() {
        val mu = 3.0
        val r = 4.0
        var prevLogP = 0.0
        for (k in 1..25) {
            val logP = NegBinUtil.logNegBinCdf(k, mu, r)
            assertTrue(logP < prevLogP, "Log p-value should strictly decrease with k: k=$k, logP=$logP, prev=$prevLogP")
            prevLogP = logP
        }
    }

    @Test
    fun testLogNegBinCdfMonotonicityInMu() {
        val k = 10
        val r = 4.0
        var prevLogP = Double.NEGATIVE_INFINITY
        for (i in 1..10) {
            val mu = i * 0.5
            val logP = NegBinUtil.logNegBinCdf(k, mu, r)
            assertTrue(logP > prevLogP, "Log p-value should strictly increase with higher background mu: mu=$mu, logP=$logP")
            prevLogP = logP
        }
    }

    @Test
    fun testLogNegBinCdfLargeKNoUnderflow() {
        val mu = 2.0
        val r = 5.0
        val logP = NegBinUtil.logNegBinCdf(100, mu, r)
        assertTrue(logP < -100.0, "Log p-value for large count should be very negative: $logP")
        assertTrue(!logP.isNaN() && !logP.isInfinite(), "Log p-value should not be NaN or Infinite: $logP")
    }

    @Test
    fun testLogNegBinCdfPoissonLimit() {
        val k = 8
        val mu = 3.0
        val rLarge = 1e8
        val actual = NegBinUtil.logNegBinCdf(k, mu, rLarge)
        val expected = PoissonUtil.logPoissonCdf(k - 1, mu)
        assertEquals(expected, actual, 1e-4, "Large r should converge to Poisson")
    }
}
