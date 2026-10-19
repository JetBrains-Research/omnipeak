package org.jetbrains.bio.omnipeak.statistics

import org.apache.commons.math3.special.Gamma
import org.jetbrains.bio.omnipeak.statistics.util.PoissonUtil
import org.junit.Test
import kotlin.math.abs
import kotlin.math.ln
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PoissonUtilTest {

    @Test
    fun testLogFactorialSmallNumbers() {
        assertEquals(0.0, PoissonUtil.logFactorial(0), 1e-12)
        assertEquals(0.0, PoissonUtil.logFactorial(1), 1e-12)
        assertEquals(ln(2.0), PoissonUtil.logFactorial(2), 1e-12)
        assertEquals(ln(6.0), PoissonUtil.logFactorial(3), 1e-12)
        assertEquals(ln(24.0), PoissonUtil.logFactorial(4), 1e-12)
        assertEquals(ln(120.0), PoissonUtil.logFactorial(5), 1e-12)
        assertEquals(ln(3628800.0), PoissonUtil.logFactorial(10), 1e-12)
    }

    @Test
    fun testLogFactorialTableRange() {
        val testValues = intArrayOf(1, 2, 10, 50, 100, 500, 1000, 5000, 10000)
        for (m in testValues) {
            val expected = Gamma.logGamma(m + 1.0)
            val actual = PoissonUtil.logFactorial(m)
            assertEquals(expected, actual, 1e-9, "Mismatch for m = $m")
        }
    }

    @Test
    fun testLogFactorialLargeNumbersApproximation() {
        val testValues = intArrayOf(10001, 10002, 10050, 15000, 20000, 50000, 100000)
        for (m in testValues) {
            val expected = Gamma.logGamma(m + 1.0)
            val actual = PoissonUtil.logFactorial(m)
            // Ramanujan's Stirling approximation has extremely low relative and absolute error
            assertEquals(expected, actual, 1e-6, "Mismatch for large m = $m")
            val relativeError = abs((actual - expected) / expected)
            assertTrue(relativeError < 1e-9, "Relative error too high for m = $m: $relativeError")
        }
    }

    @Test
    fun testLogFactorialBoundaryContinuity() {
        val v10000 = PoissonUtil.logFactorial(10000)
        val v10001 = PoissonUtil.logFactorial(10001)
        // ln(10001!) - ln(10000!) = ln(10001)
        val diff = v10001 - v10000
        val expectedStep = ln(10001.0)
        assertEquals(expectedStep, diff, 1e-6, "Step at MAX_TAB boundary should equal ln(10001)")
    }

    @Test
    fun testLogPoissonCdf() {
        // Upper tail p-value for Poisson(lambda = 10), observation k = 15
        val logP = PoissonUtil.logPoissonCdf(15, 10.0)
        assertTrue(logP < 0.0, "Log p-value should be negative")
        // Expected tail probability P(X > 15) for Poisson(10) is ~0.04874 -> ln(0.04874) ~ -3.021
        assertEquals(-3.02, logP, 0.05)

        // Large observation and lambda
        val logPLarge = PoissonUtil.logPoissonCdf(10050, 10000.0)
        assertTrue(logPLarge < 0.0, "Log p-value for large lambda should be negative")
    }
}
