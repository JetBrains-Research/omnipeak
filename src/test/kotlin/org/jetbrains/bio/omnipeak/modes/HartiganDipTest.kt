package org.jetbrains.bio.omnipeak.modes

import org.jetbrains.bio.omnipeak.modes.HartiganDip.dipTest
import org.junit.Test
import java.util.*
import kotlin.test.assertTrue

class HartiganDipTest {

    private val R = Random(42)

    @Test
    fun testBimodal() {
        val x = DoubleArray(400) { if (it < 200) R.nextGaussian() - 2 else R.nextGaussian() + 2 }
        val r = dipTest(x, 1000)
        assertTrue(r.pValue < 0.05)
    }

    @Test
    fun testTrimodal() {
        val x = DoubleArray(600) {
            if (it < 200) R.nextGaussian() - 2 else
                if (it < 400) R.nextGaussian() else R.nextGaussian() + 2
        }
        val r = dipTest(x, 1000)
        assertTrue(r.pValue < 0.05)
    }

    @Test
    fun testUnimodal() {
        val x = DoubleArray(400) { R.nextDouble() }
        val r = dipTest(x, 1000)
        assertTrue(r.pValue > 0.05)
    }
}