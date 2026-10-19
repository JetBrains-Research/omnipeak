package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.coverage.SingleEndCoverage
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.withTempDirectory
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals

class NormalizedCoverageQueryTest {
    private var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
    private var chromosome1: Chromosome = genomeQuery.get()[0]

    @Test
    fun testScoresWithControl() {
        withTempDirectory("custom_genome") { dir ->
            val chromSizes = dir / "to100.chrom.sizes"
            chromSizes.bufferedWriter().use {
                it.write("chr1\t100\n")
            }
            val genome100 = Genome[chromSizes, "to100"]
            val genomeQuery100 = GenomeQuery(genome100, "chr1")
            val chr1 = genomeQuery100.get()[0]

            withTempFile("track", ".bed.gz") { path ->
                val cond = SingleEndCoverage.builder(genomeQuery100).apply {
                    data[chr1, Strand.PLUS].addAll(intArrayOf(1, 2, 3, 4, 5, 10, 11, 15))
                }.build(unique = false)

                val control = SingleEndCoverage.builder(genomeQuery100).apply {
                    data[chr1, Strand.PLUS].addAll(intArrayOf(0, 2, 4, 6, 10, 12, 14, 20, 21, 22, 25))
                }.build(unique = false)

                val (scaleControl, beta, _) =
                    NormalizedCoverageQuery.analyzeCoverage(genomeQuery100, path, cond, path, control, 1, true)
                assertEquals(8.0 / 11.0, scaleControl, 0.01)
                assertEquals(0.30, beta, 0.01)
            }
        }
    }

    @Test
    fun testScoresWithControlNoRegression() {
        withTempDirectory("custom_genome") { dir ->
            val chromSizes = dir / "to100.chrom.sizes"
            chromSizes.bufferedWriter().use {
                it.write("chr1\t100\n")
            }
            val genome100 = Genome[chromSizes, "to100_no_reg"]
            val genomeQuery100 = GenomeQuery(genome100, "chr1")
            val chr1 = genomeQuery100.get()[0]

            withTempFile("track", ".bed.gz") { path ->
                val cond = SingleEndCoverage.builder(genomeQuery100).apply {
                    data[chr1, Strand.PLUS].addAll(intArrayOf(1, 2, 3, 4, 5, 10, 11, 15))
                }.build(unique = false)

                val control = SingleEndCoverage.builder(genomeQuery100).apply {
                    data[chr1, Strand.PLUS].addAll(intArrayOf(0, 2, 4, 6, 10, 12, 14, 20, 21, 22, 25))
                }.build(unique = false)

                val (scaleControl, beta, _) =
                    NormalizedCoverageQuery.analyzeCoverage(genomeQuery100, path, cond, path, control, 1, false)
                assertEquals(8.0 / 11.0, scaleControl, 0.01)
                assertEquals(0.0, beta, 0.01)
            }
        }
    }

    @Test
    fun testBetaEstimationWithCorrelatedBackground() {
        withTempFile("track", ".bed.gz") { path ->
            // Bin 0 (0..200): Treatment has 10 tags, Control has 20 tags
            // Bin 1 (200..400): Treatment has 5 tags, Control has 10 tags
            // Bin 2 (400..600): Treatment has 50 tags (true peak), Control has 0 tags
            val treatmentTags = mutableListOf<Int>()
            repeat(10) { treatmentTags.add(50) } // bin 0
            repeat(5) { treatmentTags.add(250) }  // bin 1
            repeat(50) { treatmentTags.add(450) } // bin 2 (peak)

            val controlTags = mutableListOf<Int>()
            repeat(20) { controlTags.add(50) }  // bin 0
            repeat(10) { controlTags.add(250) } // bin 1
            repeat(0) { controlTags.add(450) }  // bin 2

            val cond = SingleEndCoverage.builder(genomeQuery).apply {
                data[chromosome1, Strand.PLUS].addAll(treatmentTags.toIntArray())
            }.build(unique = false)

            val control = SingleEndCoverage.builder(genomeQuery).apply {
                data[chromosome1, Strand.PLUS].addAll(controlTags.toIntArray())
            }.build(unique = false)

            val (scaleControl, beta, _) =
                NormalizedCoverageQuery.analyzeCoverage(genomeQuery, path, cond, path, control, 200, true)

            // Total treatment = 65, total control = 30 -> scaleControl = 65 / 30 = 2.1667
            // In bins 0 and 1, treatment has 0.5 * control.
            // When control is scaled by 2.1667, treatment background is (0.5 / 2.1667) * scaledControl ≈ 0.23 * scaledControl.
            // Background subtraction should find beta > 0.
            assertEquals(65.0 / 30.0, scaleControl, 0.01)
            assertEquals(0.23, beta, 0.05)
        }
    }

}


