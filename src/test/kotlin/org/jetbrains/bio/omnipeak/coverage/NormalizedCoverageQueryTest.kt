package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.coverage.SingleEndCoverage
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals

class NormalizedCoverageQueryTest {
    private var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
    private var chromosome1: Chromosome = genomeQuery.get()[0]

    @Test
    fun testScoresWithControl() {
        withTempFile("track", ".bed.gz") { path ->
            val cond = SingleEndCoverage.builder(genomeQuery).apply {
                data[chromosome1, Strand.PLUS].addAll(intArrayOf(1, 2, 3, 4, 5, 10, 11, 15))
            }.build(unique = false)

            val control = SingleEndCoverage.builder(genomeQuery).apply {
                data[chromosome1, Strand.PLUS].addAll(intArrayOf(0, 2, 4, 6, 10, 12, 14, 20, 21, 22, 25))
            }.build(unique = false)

            val (scaleControl, beta, _) =
                NormalizedCoverageQuery.analyzeCoverage(genomeQuery, path, cond, path, control, 200, true)
            assertEquals(0.72, scaleControl, 0.01)
            assertEquals(0.27, beta, 0.01)
        }
    }

    @Test
    fun testScoresWithControlNoRegression() {
        withTempFile("track", ".bed.gz") { path ->
            val cond = SingleEndCoverage.builder(genomeQuery).apply {
                data[chromosome1, Strand.PLUS].addAll(intArrayOf(1, 2, 3, 4, 5, 10, 11, 15))
            }.build(unique = false)

            val control = SingleEndCoverage.builder(genomeQuery).apply {
                data[chromosome1, Strand.PLUS].addAll(intArrayOf(0, 2, 4, 6, 10, 12, 14, 20, 21, 22, 25))
            }.build(unique = false)

            val (scaleControl, beta, _) =
                NormalizedCoverageQuery.analyzeCoverage(genomeQuery, path, cond, path, control, 200, false)
            assertEquals(0.72, scaleControl, 0.01)
            assertEquals(0.0, beta, 0.01)
        }
    }

}


