package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.Tests.assertMatches
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.omnipeak.coverage.BigWigCoverageWriter
import org.jetbrains.bio.omnipeak.coverage.CoverageSampler.sampleCoverage
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_BIN
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FDR
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.omnipeak.fit.OmnipeakDataPaths
import org.jetbrains.bio.omnipeak.fit.OmnipeakModelFitExperiment
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.FileReader
import java.text.DecimalFormatSymbols
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue

class OmnipeakLongTest {

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
    fun emptyArgs() {
        val (_, err) = Logs.captureLoggingOutput {
            OmnipeakCLA.main(arrayOf())
        }
        assertTrue("ERROR: No command given." in err)
    }

    @Test
    fun illegalArgs() {
        val (_, err) = Logs.captureLoggingOutput {
            OmnipeakCLA.main(arrayOf("foobar"))
        }
        assertTrue("ERROR: Unknown command: foobar." in err)
    }

    @Test
    fun quietError() {
        val (_, err) = Logs.captureLoggingOutput {
            OmnipeakCLA.main(arrayOf("foobar", "quiet"))
        }
        assertTrue("ERROR: Unknown command: foobar." in err)
    }


    @Test
    fun checkVersion() {
        val (out, _) = Logs.captureLoggingOutput {
            OmnipeakCLA.main(arrayOf("--version"))
        }
        val version = out.trim()
        if (version == "@VERSION@.@BUILD@ built on @DATE@") {
            return
        }
        // the test is sometimes launched in the assembled JAR, where the @@ tokens have already been substituted
        assertMatches(
            version,
            Regex("^[0-9]+(\\.[0-9]+)+(\\.build)? built on [A-Z][a-z]* [0-9]+, [0-9]{4}")
        )
    }


    @Test
    fun compareSameTestOrganismTracks() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val peaksPath = it / "peaks.bed"
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "compare",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t1", path.toString(),
                            "-t2", path.toString(),
                            "--peaks", peaksPath.toString(),
                            "--fdr", OMNIPEAK_DEFAULT_FDR.toString(),
                            "--threads", THREADS.toString()
                        )
                    )
                }

                assertTrue(
                    peaksPath.size.isEmpty(),
                    "Found differential peaks in identical signals."
                )

                assertIn(
                    """
COMMAND:
LOG:
WORKING DIR: $it
THREADS: $THREADS
TREATMENT1: $path
CONTROL1: none
TREATMENT2: $path
CONTROL2: none
CHROM.SIZES: $chromsizes
FRAGMENT: auto
BIN: $OMNIPEAK_DEFAULT_BIN
FDR: $OMNIPEAK_DEFAULT_FDR
PEAKS: $peaksPath
""", out
                )
                assertIn("Saved result to $peaksPath", out)

                // Check model fit has progress:
                // XXX: Not so important to make to types of tests for US and EU locales
                val ds = DecimalFormatSymbols(Locale.getDefault()).decimalSeparator
                assertIn("0${ds}00% (0/${OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS}), Elapsed time", out)
                assertIn("100${ds}00% (", out)
            }
        }
    }

    @Test
    fun compareSameTestOrganismTracksReplicates() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val bedPath = it / "peaks.bed"
                OmnipeakCLA.main(
                    arrayOf(
                        "compare",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", it.toString(),
                        "-b", OMNIPEAK_DEFAULT_BIN.toString(),
                        "-fragment", FRAGMENT.toString(),
                        "-t1", "$path,$path",
                        "-t2", "$path,$path,$path",
                        "--peaks", bedPath.toString(),
                        "--fdr", OMNIPEAK_DEFAULT_FDR.toString()
                    )
                )
                assertTrue(
                    bedPath.size.isEmpty(),
                    "Found differential peaks in identical signals."
                )
            }
        }
    }

    @Test
    fun testModelFitFile() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--keep-cache",
                        )
                    )
                }
                assertIn(
                    "NO peaks path given, process model fitting only.\n" +
                            "Fdr, sensitivity, gap, clip options are ignored.", out
                )
                assertIn(": done in ", out)
                assertIn("Model saved: ", out)
                assertFalse("Loading model" in out)
                assertFalse("Completed loading model" in out)
            }
        }
    }

    @Test
    fun testBadTrackQualityWarning() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, 200, goodQuality = false)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--bin", "200",
                            "--debug"
                        )
                    )
                }
                assertTrue("Emission's parameter p in LOW state" in out)
                assertTrue("Low quality of data detected during fitting the model." in out)
            }
        }
    }

    @Test
    fun testUnrecognizedFileFormat() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".foobar.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, err) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--format", "BED",
                            "--keep-cache"
                        )
                    )
                }
                assertIn("Done computing data model", out)
                assertEquals("", err)
            }
        }
    }


    @Test
    fun testWrongFileFormat() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (_, err) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--format", "BAM",
                        )
                    )
                }
                assertIn(
                    "Error parsing text SAM file.",
                    err
                )
            }
        }
    }

    @Test
    fun testQuietMode() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = false)
            println("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                /* We can't test the stdout, because the Omnipeak quiet mode redefines the JVM [System.out].
                 * But we can restore the System.out to the original value using [captureLoggingOutput].
                 */
                Logs.captureLoggingOutput {
                    val oldSystemOut = System.out
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "-q"
                        )
                    )
                    /* the best we can do is to check whether any redirection took place */
                    assertNotEquals(
                        oldSystemOut, System.out,
                        "Omnipeak quiet mode didn't redirect System.out"
                    )
                }
                // Log file
                val id = OmnipeakAnalyzeFitInformation.generateId(
                    listOf(OmnipeakDataPaths(path, null)),
                    AutoFragment,
                    OMNIPEAK_DEFAULT_BIN,
                    unique = true, regressControl = true
                )
                val logPath = Configuration.logsPath.glob("${id}*.log").first()
                assertTrue(logPath.exists, "Log file not found")
                assertTrue(logPath.size.isNotEmpty(), "Log file is empty")
            }
        }
    }


    @Test
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                        )
                    )

                    // Model test
                    val id = OmnipeakAnalyzeFitInformation.generateId(
                        listOf(OmnipeakDataPaths(path, null)),
                        AutoFragment,
                        OMNIPEAK_DEFAULT_BIN,
                        unique = true, regressControl = true
                    )
                    assertEquals(0, Configuration.experimentsPath.glob("${id}*.omni").size)
                    // Log file
                    assertEquals(0, Configuration.logsPath.glob("${id}*.log").size)
                    // Genome Coverage test
                    assertEquals(0, Configuration.cachesPath.glob("coverage_${path.stemGz}_unique#*.npz").size)
                    assertEquals(0, Configuration.cachesPath.glob("coverage_${control.stemGz}_unique#*.npz").size)
                }
            }
        }
    }


    @Test
    fun testFilesCreatedByAnalyzeKeepCache() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--keep-cache"
                        )
                    )

                    // Model test
                    val id = OmnipeakAnalyzeFitInformation.generateId(
                        listOf(OmnipeakDataPaths(path, control)),
                        AutoFragment,
                        OMNIPEAK_DEFAULT_BIN,
                        unique = true, regressControl = true
                    )
                    assertEquals(1, Configuration.experimentsPath.glob("${id}*.omni").size)
                    // Log file
                    assertEquals(1, Configuration.logsPath.glob("${id}*.log").size)
                    // Genome Coverage test
                    assertEquals(1, Configuration.cachesPath.glob("coverage_${path.stemGz}_unique#*.npz").size)
                    assertEquals(1, Configuration.cachesPath.glob("coverage_${control.stemGz}_unique#*.npz").size)
                }
            }
        }
    }


    @Test
    fun testCustomModelPath() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    val modelPath = dir / "custom" / "path" / "model.omni"
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--model", modelPath.toString()
                        )
                    )
                    assertTrue(modelPath.exists, "Model was not created at $modelPath")
                    assertTrue(modelPath.size.isNotEmpty(), "Model file $modelPath is empty")
                    val (reloadOut, reloadErr) = Logs.captureLoggingOutput {
                        OmnipeakCLA.main(
                            arrayOf(
                                "analyze",
                                "--workdir", dir.toString(),
                                "--threads", THREADS.toString(),
                                "--model", modelPath.toString()
                            )
                        )
                    }
                    assertIn("MODEL: $modelPath", reloadOut)
                    assertEquals("", reloadErr)

                    val (_, invalidErr) = Logs.captureLoggingOutput {
                        OmnipeakCLA.main(
                            arrayOf(
                                "analyze",
                                "--workdir", dir.toString(),
                                "--threads", THREADS.toString(),
                                "--model", modelPath.toString(),
                                "--bin", "137"
                            )
                        )
                    }
                    assertIn(
                        "bin size (${OMNIPEAK_DEFAULT_BIN}) differs from the command line argument (137)",
                        invalidErr
                    )
                }
            }
        }
    }

    @Test
    fun testCustomLogPath() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    val logPath = dir / "custom" / "logs" / "log.txt"
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--log", logPath.toString()
                        )
                    )
                    assertTrue(logPath.exists, "Log was not created at $logPath")
                    assertTrue(logPath.size.isNotEmpty(), "Log file $logPath is empty")
                }
            }
        }
    }


    @Test
    fun testAnalyzePeaksOrModelOrKeepCacheRequired() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (_, invalidErr) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                        )
                    )
                }
                assertIn(
                    "ERROR: At least one of the parameters is required: --peaks, --model or --keep-cache.",
                    invalidErr
                )
            }
        }
    }

    @Test
    fun testComparePeaksOrKeepCacheRequired() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (_, invalidErr) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "compare",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "--t1", path.toString(),
                            "--t2", path.toString(),
                        )
                    )
                }
                assertIn(
                    "ERROR: At least one of the parameters is required: --peaks or --keep-cache.",
                    invalidErr
                )
            }
        }
    }


    @Test
    fun testTypeOnlyValidInExperimentalMode() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)

                val chromsizes = Genome["to1"].chromSizesPath.toString()

                val (_, wrongErr) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--type", "nbhmm"
                        )
                    )
                }
                assertIn(
                    "ERROR: type is not a recognized option",
                    wrongErr
                )
            }
        }
    }


    @Test
    fun testAnalyze() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, OMNIPEAK_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val peaksPath = path.parent / "${path.stem}.peak"
                val (out, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze", "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--peaks", peaksPath.toString()
                        )
                    )
                }

                val ds =
                    DecimalFormatSymbols(Locale.getDefault()).decimalSeparator // XXX: Not so important to make to types of tests for US and EU locales
                assertIn(
                    """Omnipeak
COMMAND:
LOG:
WORKING DIR: $it
THREADS: $THREADS
TREATMENT: $path
CONTROL: none
CHROM.SIZES: $chromsizes
FRAGMENT: auto
MAX ITERATIONS: $OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS
CONVERGENCE THRESHOLD: $OMNIPEAK_DEFAULT_FIT_THRESHOLD
EXTENDED MODEL INFO: false
Library: ${path.fileName}, Depth:
100${ds}00% (
File: $peaksPath
FRIP: 
Reads: single-ended, Fragment size: 2 bp (cross-correlation estimate)
""", out
                )
                assertFalse(
                    "NO peaks path given, process model fitting only.\n" +
                            "Labels, fdr, background sensitivity, clip options are ignored." in out
                )

                /* Check that coverage is being generated */
                val format = BedFormat.from("bed6+3")
                assertTrue(
                    format.parse(peaksPath) { parser ->
                        parser.all { entry ->
                            val coverage = entry.unpack(6).extraFields?.get(0)
                            return@all coverage != null && coverage != "0.0"
                        }
                    },
                    "Peak value is reported as 0.0, although the coverage cache is present"
                )
                assertIn("Signal mean: ", out)
                assertIn("Noise mean: ", out)
                assertIn("Signal to noise: ", out)
            }
        }
    }


    @Test
    fun analyzeSampledEnrichment() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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
            sampleCoverage(
                path,
                TO,
                OMNIPEAK_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            withTempDirectory("work") { dir ->
                val bedPath = dir / "result.bed"
                OmnipeakCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", OMNIPEAK_DEFAULT_FDR.toString(),
                        "-t", path.toString()
                    )
                )
                // Check created bed file
                assertTrue(
                    Location(
                        1100 * OMNIPEAK_DEFAULT_BIN,
                        1900 * OMNIPEAK_DEFAULT_BIN,
                        TO.get().first()
                    )
                            in LocationsMergingList.load(TO, bedPath),
                    "Expected location not found in called peaks"
                )
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
    fun testAnalyzeEnrichmentWithControl() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempDirectory("track") { dir ->
            val path = dir / "track.bed.gz"
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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

            sampleCoverage(
                path, TO, OMNIPEAK_DEFAULT_BIN, enrichedRegions, zeroRegions,
                goodQuality = true,
                totalCoverage = (TO.get().first().length * 0.01).toLong()
            )
            println("Saved sampled track file: $path")
            val control = dir / "control.bed.gz"
            sampleCoverage(
                control, TO, OMNIPEAK_DEFAULT_BIN, enrichedRegions, zeroRegions,
                goodQuality = false,
                totalCoverage = (TO.get().first().length * 0.1).toLong()
            )
            println("Saved sampled control track file: $control")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val modelPath = it / "model.omni"
                val peaksPath = it / "${path.stem}.peak"
                val (out, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze", "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--model", modelPath.toString(),
                            "--peaks", peaksPath.toString(),
                        )
                    )
                    // Check created bed file
                    assertTrue(
                        Location(
                            1100 * OMNIPEAK_DEFAULT_BIN,
                            1900 * OMNIPEAK_DEFAULT_BIN,
                            TO.get().first()
                        )
                                in LocationsMergingList.load(TO, peaksPath),
                        "Expected location not found in called peaks"
                    )

                    val fitResults = OmnipeakModelFitExperiment.loadResults(Genome["to1"].toQuery(), modelPath)
                    assertTrue(fitResults.fitInfo.regressControl)

                    val data = fitResults.fitInfo.dataQuery.apply(TO.get().first())
                    val scores = data.sliceAsInt(data.labels.first())
                    assertTrue((1100 until 1900).map { scores[it] }.all { it < 100 })
                }

                assertIn("CONTROL: $control", out)
                assertIn("NO CONTROL REGRESSION: false", out)
            }
        }
    }

    @Test
    fun testAnalyzeEnrichmentWithControlNoControlRegression() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempDirectory("track") { dir ->
            val path = dir / "track.bed.gz"
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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

            sampleCoverage(
                path, TO, OMNIPEAK_DEFAULT_BIN, enrichedRegions, zeroRegions,
                goodQuality = true, totalCoverage = (TO.get().first().length * 0.01).toLong()
            )
            println("Saved sampled track file: $path")
            val control = dir / "control.bed.gz"
            sampleCoverage(
                control, TO, OMNIPEAK_DEFAULT_BIN, enrichedRegions, zeroRegions,
                goodQuality = false,
                totalCoverage = (TO.get().first().length * 0.1).toLong()
            )
            println("Saved sampled control track file: $control")

            val chromsizes = Genome["to1"].chromSizesPath.toString()
            val peaksPath = dir / "result.peak"
            val modelPath = dir / "model.omni"
            val (out, _) = Logs.captureLoggingOutput {
                OmnipeakCLA.main(
                    arrayOf(
                        "analyze", "-cs", chromsizes,
                        "--workdir", dir.toString(),
                        "-t", path.toString(),
                        "-c", control.toString(),
                        "--threads", THREADS.toString(),
                        "--model", modelPath.toString(),
                        "--peaks", peaksPath.toString(),
                        "--no-control-regression",
                    )
                )
            }
            // Check created bed file
            assertTrue(
                Location(
                    1100 * OMNIPEAK_DEFAULT_BIN,
                    1900 * OMNIPEAK_DEFAULT_BIN,
                    TO.get().first()
                )
                        in LocationsMergingList.load(TO, peaksPath),
                "Expected location not found in called peaks"
            )

            val fitResults = OmnipeakModelFitExperiment.loadResults(Genome["to1"].toQuery(), modelPath)
            assertFalse(fitResults.fitInfo.regressControl)

            val data = fitResults.fitInfo.dataQuery.apply(TO.get().first())
            val scores = data.sliceAsInt(data.labels.first())
            assertTrue((1100 until 1900).map { scores[it] }.all { it == 100 })

            assertIn("NO CONTROL REGRESSION: true", out)
        }
    }


    @Test
    fun analyzeSampledBigWigEnrichment() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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
                        "-fdr", OMNIPEAK_DEFAULT_FDR.toString(),
                        "-t", bigWigPath.toString()
                    )
                )
                // Check created bed file
                assertTrue(
                    Location(
                        1100 * OMNIPEAK_DEFAULT_BIN,
                        1900 * OMNIPEAK_DEFAULT_BIN,
                        TO.get().first()
                    )
                            in LocationsMergingList.load(TO, bedPath),
                    "Expected location not found in called peaks"
                )
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
    fun checkNegativeBigWig() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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


    @Test
    fun analyzeSampledEnrichmentReplicated() {
        val enrichedRegions = genomeMap(TO) {
            val enriched = BitSet()
            if (it.name == "chr1") {
                enriched.set(1000, 2000)
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
        withTempFile("track1", ".bed.gz") { coveragePath1 ->
            withTempFile("track2", ".bed.gz") { coveragePath2 ->
                sampleCoverage(
                    coveragePath1,
                    TO,
                    OMNIPEAK_DEFAULT_BIN,
                    enrichedRegions,
                    zeroRegions,
                    goodQuality = true
                )
                println("Saved sampled track file 1: $coveragePath1")

                sampleCoverage(
                    coveragePath2,
                    TO,
                    OMNIPEAK_DEFAULT_BIN,
                    enrichedRegions,
                    zeroRegions,
                    goodQuality = true
                )
                println("Saved sampled track file 2: $coveragePath2")

                withTempDirectory("work") { dir ->
                    val peaksPath = dir / "peaks_rep.bed"
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "--peaks", peaksPath.toString(),
                            "-fdr", OMNIPEAK_DEFAULT_FDR.toString(),
                            "-t", "$coveragePath1,$coveragePath2"
                        )
                    )
                    // Check created bed file
                    val peaksLocations = LocationsMergingList.load(TO, peaksPath)
                    assertTrue(
                        Location(
                            1100 * OMNIPEAK_DEFAULT_BIN,
                            1900 * OMNIPEAK_DEFAULT_BIN,
                            TO.get().first()
                        ) in peaksLocations,
                        "Expected location not found in called peaks"
                    )
                }
            }
        }
    }


    @Test
    fun analyzeSampledEnrichmentReusingModel() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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
            sampleCoverage(
                path,
                TO,
                OMNIPEAK_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            withTempDirectory("work") { dir ->
                val bedPath = dir / "result.bed"
                val modelPath = dir / "model.omni"
                OmnipeakCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "-m", modelPath.toString(),
                        "-t", path.toString(),
                        "-keep-cache"
                    )
                )
                OmnipeakCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-m", modelPath.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", OMNIPEAK_DEFAULT_FDR.toString()
                    )
                )
                // Check created bed file
                assertTrue(
                    Location(
                        1100 * OMNIPEAK_DEFAULT_BIN,
                        1900 * OMNIPEAK_DEFAULT_BIN,
                        TO.get().first()
                    ) in LocationsMergingList.load(TO, bedPath),
                    "Expected location not found in called peaks"
                )
            }
        }
    }

    @Test
    fun analyzeSampledEnrichmentReusingModelNoKeepCache() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
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
            sampleCoverage(
                path,
                TO,
                OMNIPEAK_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            withTempDirectory("work") { dir ->
                val (out1, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "-t", path.toString(),
                            "-d"
                        )
                    )
                }
                assertIn("Model is not saved", out1)
                val (out2, _) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "-t", path.toString(),
                        )
                    )
                }
                assertIn("Model is not saved", out2)
            }
        }
    }


    @Test
    fun analyzeEmptyCoverage() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) { BitSet() }

            val zeroRegions = genomeMap(TO) {
                val zeroes = BitSet()
                zeroes.set(0, it.length / OMNIPEAK_DEFAULT_BIN)
                zeroes
            }



            sampleCoverage(
                path,
                TO,
                OMNIPEAK_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                val (out, err) = Logs.captureLoggingOutput {
                    OmnipeakCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "-t", path.toString(),
                        )
                    )
                }
                val id = OmnipeakAnalyzeFitInformation.generateId(
                    listOf(OmnipeakDataPaths(path, null)),
                    AutoFragment,
                    OMNIPEAK_DEFAULT_BIN,
                    unique = true, regressControl = true
                )
                // Log file
                assertEquals(1, Configuration.logsPath.glob("${id}*.log").size)
                val log = FileReader(Configuration.logsPath.glob("${id}*.log").first().toFile())
                    .use { it.readText() }
                val errorMessage = "Model can't be trained on empty coverage, exiting."
                assertIn(errorMessage, log)
                assertIn(errorMessage, out)
                assertIn(errorMessage, err)
            }
        }
    }

    @Test
    fun analyzePartiallyEmptyCoverage() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(
                path,
                GenomeQuery(Genome["to1"], "chr1", "chr2"),
                OMNIPEAK_DEFAULT_BIN,
                goodQuality = true
            )
            println("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                OmnipeakCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "-t", path.toString()
                    )
                )
            }
        }
    }

    companion object {
        internal val TO = GenomeQuery(Genome["to1"])
        internal const val THREADS = 1
        private const val FRAGMENT = 200
    }
}
