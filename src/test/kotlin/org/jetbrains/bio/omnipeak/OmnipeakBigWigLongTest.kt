package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.omnipeak.coverage.BigWigCoverageWriter
import org.jetbrains.bio.omnipeak.coverage.CoverageSampler.sampleCoverage
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_BIN
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.statistics.model.Boo
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.FileReader
import java.nio.file.Path
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * BigWig-specific peak calling tests.
 *
 * BigWig tracks carry continuous signal that may be stored on very different scales
 * (raw counts, RPKM, CPM, fold change, ...). [org.jetbrains.bio.omnipeak.coverage.BigWigBinnedCoverageQuery]
 * normalizes the signal internally and computes range scores as a genuine count over
 * the range, so that the Poisson peak scoring stays valid and FDR-sensitive regardless
 * of the input scale.
 */
class OmnipeakBigWigLongTest {

    @Before
    fun setUp() {
        OmnipeakCLA.ignoreConfigurePaths = true
        Sampling.RANDOM_DATA_GENERATOR.randomGenerator.setSeed(1234L)
        System.setProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true")
    }

    @After
    fun tearDown() {
        OmnipeakCLA.ignoreConfigurePaths = false
        // we might have unfinished tracked tasks which will never be complete, let's drop them
        MultitaskProgress.clear()
        System.setProperty(JOPTSIMPLE_SUPPRESS_EXIT, "false")
    }

    @Test
    fun analyzeSampledBigWigEnrichment() {
        withTempFile("track", ".bed.gz") { path ->
            val (enrichedRegions, zeroRegions) = markRegions()
            sampleCoverage(
                path,
                TO,
                OMNIPEAK_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            val genome = Genome["to1"]
            val genomeQuery = genome.toQuery()
            val coverage = ReadsQuery(genomeQuery, path, null).get()

            withTempDirectory("work") { dir ->
                val bigWigPath = dir / "data.bw"
                BigWigCoverageWriter.write(
                    { cr -> coverage.getBothStrandsCoverage(cr).toDouble() },
                    genomeQuery,
                    OMNIPEAK_DEFAULT_BIN,
                    bigWigPath
                )
                val bedPath = dir / "result.bed"
                OmnipeakCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", genome.chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "--peaks", bedPath.toString(),
                        "-t", bigWigPath.toString()
                    )
                )
                // Check created bed file
                checkLocations(bedPath)
                // Check correct log file name
                val logPath = Configuration.logsPath / "${bedPath.stem}.log"
                assertTrue(logPath.exists, "Log file not found")
                val log = FileReader(logPath.toFile()).use { it.readText() }
                assertIn("Signal mean:", log)
                assertIn("Noise mean:", log)
                assertIn("Signal to noise:", log)
                assertTrue(log.substringAfter("Signal to noise:").substringBefore("\n").trim().toDouble() > 5)
            }
        }
    }

    @Test
    fun analyzeBigWigScaleInvariance() {
        // BigWig may store raw counts, RPKM, CPM, etc. The peak caller normalizes the
        // input scale internally (see BigWigBinnedCoverageQuery.treatmentScale), so the
        // number of called peaks must stay stable across a wide range of input scales.
        withTempFile("track", ".bed.gz") { path ->
            val (enrichedRegions, zeroRegions) = markRegions()
            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, enrichedRegions, zeroRegions, goodQuality = true)
            val genome = Genome["to1"]
            val genomeQuery = genome.toQuery()
            val coverage = ReadsQuery(genomeQuery, path, null).get()

            withTempDirectory("work") { dir ->
                fun callPeaks(factor: Double): Int {
                    val bigWigPath = dir / "data_$factor.bw"
                    BigWigCoverageWriter.write(
                        { cr -> coverage.getBothStrandsCoverage(cr).toDouble() * factor },
                        genomeQuery,
                        OMNIPEAK_DEFAULT_BIN,
                        bigWigPath
                    )
                    val bedPath = dir / "result_$factor.bed"
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", genome.chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "--peaks", bedPath.toString(),
                            "-t", bigWigPath.toString()
                        )
                    )
                    // The expected enriched regions must be detected regardless of scale
                    checkLocations(bedPath)
                    return bedPath.toFile().readLines().count { it.isNotBlank() }
                }

                // Cover a 1000x span of input scales (down-scaled, raw and up-scaled)
                val counts = listOf(0.01, 0.1, 1.0, 10.0, 100.0).associateWith { callPeaks(it) }
                println("BigWig scale invariance peak counts: $counts")
                val min = counts.values.min()
                val max = counts.values.max()
                assertTrue(min > 0, "No peaks called at some scale: $counts")
                assertTrue(
                    max.toDouble() / min <= 1.5,
                    "Peak count is not scale-invariant across input scales: $counts"
                )
            }
        }
    }

    @Test
    fun analyzeBigWigFdrSensitivity() {
        // BigWig signal is not count data, so peak scoring does not use a Poisson enrichment
        // test; significance comes from the model posterior error probabilities scored via
        // the PEP-based FDR (see PeakScorer.coverageIsPoisson / Fdr.qvalidatePEPs). The
        // resulting scores must still span a meaningful range so that the number of called
        // peaks responds to the FDR threshold - a relaxed FDR must yield more peaks than a
        // strict one.
        withTempFile("track", ".bed.gz") { path ->
            val (enrichedRegions, zeroRegions) = markRegions()
            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, enrichedRegions, zeroRegions, goodQuality = true)
            val genome = Genome["to1"]
            val genomeQuery = genome.toQuery()
            val coverage = ReadsQuery(genomeQuery, path, null).get()

            withTempDirectory("work") { dir ->
                val bigWigPath = dir / "data.bw"
                BigWigCoverageWriter.write(
                    { cr -> coverage.getBothStrandsCoverage(cr).toDouble() },
                    genomeQuery,
                    OMNIPEAK_DEFAULT_BIN,
                    bigWigPath
                )

                fun callPeaks(fdr: Double, checkExpected: Boolean = true): Int {
                    val bedPath = dir / "result_$fdr.bed"
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", genome.chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "--peaks", bedPath.toString(),
                            "-t", bigWigPath.toString(),
                            "--fdr", fdr.toString()
                        )
                    )
                    if (checkExpected) checkLocations(bedPath)
                    return bedPath.toFile().readLines().count { it.isNotBlank() }
                }

                val relaxed = callPeaks(0.05)
                val strict = callPeaks(1e-10, checkExpected = false)
                val nightmare = callPeaks(1e-100, checkExpected = false)
                println("BigWig FDR sensitivity: fdr=0.05 -> $relaxed, fdr=1e-10 -> $strict, fdr=1e-100 -> $nightmare")
                assertTrue(relaxed > strict)
                assertTrue(strict > nightmare)

            }
        }
    }

    @Test
    fun writeBigWigAfterCoverageCleanup() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = false)

                    val peaksPath = dir / "peaks.bed"
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--peaks", peaksPath.toString(),
                            "--bigwig",
                        )
                    )

                    val bigWigPath = (peaksPath.toString() + ".bw").toPath()
                    assertEquals(0, Configuration.cachesPath.glob("coverage_${path.stemGz}_unique#*.npz").size)
                    assertEquals(0, Configuration.cachesPath.glob("coverage_${control.stemGz}_unique#*.npz").size)
                    assertTrue(peaksPath.exists, "Peaks were not created at $peaksPath")
                    assertTrue(peaksPath.size.isNotEmpty(), "Peaks file $peaksPath is empty")
                    assertTrue(bigWigPath.exists, "BigWig was not created at $bigWigPath")
                    assertTrue(bigWigPath.size.isNotEmpty(), "BigWig file $bigWigPath is empty")
                }
            }
        }
    }

    @Test
    fun checkNegativeBigWig() {
        withTempFile("track", ".bed.gz") { path ->
            val (enrichedRegions, zeroRegions) = markRegions()
            sampleCoverage(
                path,
                TO,
                OMNIPEAK_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            val genome = Genome["to1"]
            val genomeQuery = genome.toQuery()
            val coverage = ReadsQuery(genomeQuery, path, null).get()

            withTempDirectory("work") { dir ->
                val bigWigPath = dir / "data.bw"
                BigWigCoverageWriter.write(
                    { cr -> -coverage.getBothStrandsCoverage(cr).toDouble() },
                    genomeQuery,
                    OMNIPEAK_DEFAULT_BIN,
                    bigWigPath
                )
                val bedPath = dir / "result.bed"
                val (out, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", genome.chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "-t", bigWigPath.toString(),
                            "--peaks", bedPath.toString(),
                        )
                    )
                }
                assert("negative values detected" in out)
            }
        }
    }

    private fun markRegions(): Pair<GenomeMap<BitSet>, GenomeMap<BitSet>> {
        val enrichedRegions = genomeMap(TO) {
            val enriched = BitSet()
            if (it.name == "chr1") {
                enriched.set(100, 150)
                enriched.set(200, 300)
                enriched.set(500, 700)
            }
            enriched
        }

        val zeroRegions = genomeMap(TO) {
            val zeroes = BitSet()
            if (it.name == "chr1") {
                zeroes[3000] = 4000
            }
            zeroes
        }
        return Pair(enrichedRegions, zeroRegions)
    }

    private fun checkLocations(bedPath: Path) {
        val chromosome = TO.get().first()
        val peaks = LocationsMergingList.load(TO, bedPath)
        listOf(
            Location(100 * OMNIPEAK_DEFAULT_BIN, 150 * OMNIPEAK_DEFAULT_BIN, chromosome),
            Location(200 * OMNIPEAK_DEFAULT_BIN, 300 * OMNIPEAK_DEFAULT_BIN, chromosome),
            Location(500 * OMNIPEAK_DEFAULT_BIN, 700 * OMNIPEAK_DEFAULT_BIN, chromosome),
        ).forEach {
            assertTrue(peaks.intersects(it), "Expected location $it not found in called peaks")
        }
    }

    companion object {
        internal val TO = GenomeQuery(Genome["to1"])
        internal const val THREADS = 1
    }
}
