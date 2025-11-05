package org.jetbrains.bio.omnipeak.modes

import org.jetbrains.bio.omnipeak.modes.KdeModeFinder.detectModes
import org.junit.Test
import kotlin.test.assertEquals

class KdeModeFinderTest {

    val signal = doubleArrayOf(
        15.0,
        12.0,
        19.0,
        10.0,
        21.0,
        20.0,
        54.0,
        30.0,
        16.0,
        21.0,
        29.0,
        12.0,
        13.0,
        13.0,
        19.0,
        20.0,
        21.0,
        16.0,
        10.0,
        15.0,
        16.0,
        21.0,
        20.0,
        16.0,
        10.0,
        13.0,
        21.0,
        12.0,
        25.0,
        20.0,
        25.0,
        16.0,
        19.0,
        14.0,
        13.0,
        17.0,
        31.0,
        63.0,
        104.0,
        138.0,
        157.0,
        198.0,
        200.0,
        198.0,
        174.0,
        125.0,
        81.0,
        74.0,
        75.0,
        70.0,
        85.0,
        87.0,
        76.0,
        97.0,
        123.0,
        162.0,
        177.0,
        191.0,
        149.0,
        103.0,
        55.0,
        57.0,
        54.0,
        69.0,
        124.0,
        151.0,
        124.0,
        94.0,
        43.0,
        32.0,
        16.0,
        19.0,
        35.0,
        26.0,
        41.0,
        34.0,
        65.0,
        112.0,
        141.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        150.0,
        140.0,
        130.0,
        120.0,
        92.0,
        84.0,
        113.0,
        145.0,
        159.0,
        117.0,
        44.0
    )

    @Test
    fun testSignal() {
        assertEquals(8, detectModes(signal, 6).size)
        assertEquals(3, detectModes(signal, 20).size)
        assertEquals(1, detectModes(signal, 100).size)
    }

    @Test
    fun testMinLength() {
        assertEquals(8, detectModes(signal, 6, minLx = 0.5).size)
        assertEquals(6, detectModes(signal, 6, minLx = 0.8).size)
        assertEquals(6, detectModes(signal, 6, minLx = 1.0).size)
        assertEquals(0, detectModes(signal, 6, minLx = 10.0).size)
    }


    @Test
    fun testMinDistance() {
        assertEquals(8, detectModes(signal, 6, minDx = 0.5).size)
        assertEquals(6, detectModes(signal, 6, minDx = 1.0).size)
        assertEquals(1, detectModes(signal, 6, minDx = 50.0).size)
    }

}