package org.jetbrains.bio.omnipeak.peaks

import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.detectSensitivityTriangle
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.withResource
import org.junit.Test
import org.junit.runner.RunWith
import org.junit.runners.Parameterized
import org.junit.runners.Parameterized.Parameters
import kotlin.test.assertEquals


@RunWith(Parameterized::class)
class SensitivityGapTest(private val name: String, private val triangle: SensitivityGap.PepInfo) {
    @Test
    fun testDetectTriangle() {
        withResource(
            SensitivityGapTest::class.java,
            name
        ) { tablePath ->
            val sensitivities = arrayListOf<Double>()
            val candidatesNs = arrayListOf<Int>()
            val candidatesALs = arrayListOf<Double>()
            CSVFormat.TDF.parse(tablePath.bufferedReader()).use { parser ->
                parser.records.drop(1).forEach { record ->
                    sensitivities.add(record[0].toDouble())
                    candidatesNs.add(record[1].toInt())
                    candidatesALs.add(record[2].toDouble())
                }
            }
            val detected = detectSensitivityTriangle(
                sensitivities.toDoubleArray(), candidatesNs.toIntArray(), candidatesALs.toDoubleArray()
            )
            assertEquals(triangle, detected, "Incorrect triangle detected for $name")
        }

    }

    companion object {
        @JvmStatic
        @Parameters
        fun `data`(): Collection<Array<out Any>>  = listOf(
            arrayOf("cd14_chr15_peps.tsv",
                SensitivityGap.PepInfo(beforeMerge = 13, stable = 29, beforeNoise = 48)
            ),
            arrayOf("chips_H3K27ac_1_1.0.tsv",
                SensitivityGap.PepInfo(beforeMerge = 1, stable = 27, beforeNoise = 49)
            ),

        )
    }
}