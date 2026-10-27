package org.jetbrains.bio.omnipeak.peaks

import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class AutoCorrelationsTest {

    @Test
    fun testSyntheticNarrowVsBroadPeaks() {
        // 1. Narrow point spikes every 100 bins
        val narrowSignal = DoubleArray(1000)
        for (i in 0 until 1000 step 100) {
            narrowSignal[i] = 10.0
        }
        val narrowCorr = AutoCorrelations.computeAutoCorrelations(narrowSignal, maxCorrelationDistance = 20)
        val narrowAvg = narrowCorr.average()

        // 2. Broad domain blocks of 20 consecutive bins
        val broadSignal = DoubleArray(1000)
        for (i in 0 until 1000 step 100) {
            for (j in 0 until 20) {
                if (i + j < 1000) broadSignal[i + j] = 10.0
            }
        }
        val broadCorr = AutoCorrelations.computeAutoCorrelations(broadSignal, maxCorrelationDistance = 20)
        val broadAvg = broadCorr.average()

        assertTrue(broadAvg > narrowAvg, "Broad signal must have higher average autocorrelation than narrow signal: broad=$broadAvg, narrow=$narrowAvg")
    }

    @Test
    fun testConstantSignal() {
        val constant = DoubleArray(20) { 1.0 }
        val corr = AutoCorrelations.computeAutoCorrelations(constant, maxCorrelationDistance = 5)
        assertEquals(1.0, corr[0], 1e-6)
    }
}
