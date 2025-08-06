package org.jetbrains.bio.omnipeak

import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.experiment.configurePaths
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.omnipeak.OmnipeakCLA.LOG
import org.jetbrains.bio.omnipeak.OmnipeakCLA.checkGenomeInFitInformation
import org.jetbrains.bio.omnipeak.SpanResultsAnalysis.doDeepAnalysis
import org.jetbrains.bio.omnipeak.coverage.BigWigCoverageWriter
import org.jetbrains.bio.omnipeak.fit.*
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FRAGMENTATION_HARD
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FRAGMENTATION_LIGHT
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FRAGMENTATION_SPEED
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_MULTIPLE_TEST_CORRECTION
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.printConstants
import org.jetbrains.bio.omnipeak.peaks.MultipleTesting
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks
import org.jetbrains.bio.omnipeak.peaks.OmnipeakResult
import org.jetbrains.bio.omnipeak.peaks.Peak
import org.jetbrains.bio.util.*
import org.slf4j.event.Level
import java.nio.file.Path


object OmnipeakCLAAnalyze {

    internal fun analyze(params: Array<String>) {
        with(OmnipeakCLA.getOptionParser()) {
            acceptsAll(
                listOf("t", "treatment"),
                "ChIP-seq treatment file. If multiple files are given, treated as replicates"
            )
                .requiredUnless("model")
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("c", "control"),
                "Control file. Single control file or separate file per each treatment file required"
            )
                .availableIf("treatment")
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.noCheck())

            accepts("bigwig", "Create beta-control corrected counts per million normalized track")
                .availableIf("peaks")

            accepts(
                "ext",
                """
                    Save extended states information to model file.
                    Required for model visualization in JBR Genome Browser
                    """.trimIndent()
            )

            accepts(
                "deep-analysis",
                "Deep analysis of model including analysis of coverage / candidates / peaks"
            )

            accepts(
                "hmm-snr",
                "Fraction of coverage to estimate and guard signal to noise ratio, 0 to disable"
            )
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR)

            accepts(
                "hmm-low",
                "Minimal low state mean threshold, guards against too broad peaks, 0 to disable"
            )
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD)

            accepts(
                "f-light",
                "Lightest fragmentation threshold to enable gap compensation"
            )
                .availableUnless("gap")
                .availableUnless("summits")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(OMNIPEAK_DEFAULT_FRAGMENTATION_LIGHT)

            accepts(
                "f-hard",
                "Hardest fragmentation threshold to stop gap compensation"
            )
                .availableUnless("gap")
                .availableUnless("summits")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(OMNIPEAK_DEFAULT_FRAGMENTATION_HARD)


            accepts(
                "f-speed",
                "Minimal fragmentation speed threshold for compensation"
            )
                .availableUnless("gap")
                .availableUnless("summits")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(OMNIPEAK_DEFAULT_FRAGMENTATION_SPEED)

            parse(params) { options ->
                if ("quiet" in options) {
                    Logs.quiet()
                } else {
                    Logs.addConsoleAppender(if ("debug" in options) Level.DEBUG else Level.INFO)
                }
                LOG.info("Omnipeak ${OmnipeakCLA.version()}")
                LOG.info(
                    "COMMAND: analyze ${params.joinToString(" ")}"
                )

                OmnipeakCLA.checkMemory()
                val peaksPath = options.valueOf("peaks") as Path?
                val modelPath = options.valueOf("model") as Path?
                val deepAnalysis = "deep-analysis" in options
                val blackListPath = options.valueOf("blacklist") as Path?
                val keepCacheFiles = "keep-cache" in options
                checkOrFail(peaksPath != null || modelPath != null || keepCacheFiles) {
                    "At least one of the parameters is required: --peaks, --model or --keep-cache."
                }

                val modelId = peaksPath?.stemGz ?: modelPath?.stem ?: OmnipeakAnalyzeFitInformation.generateId(
                    prepareAndCheckTreatmentControlPaths(options),
                    OmnipeakCLA.getFragment(options),
                    OmnipeakCLA.getBin(options),
                    OmnipeakCLA.getUnique(options),
                )
                val fdr = options.valueOf("fdr") as Double
                require(0 < fdr && fdr < 1) { "Illegal fdr: $fdr, expected range: (0, 1)" }
                val sensitivity = if (options.has("sensitivity")) options.valueOf("sensitivity") as Double else null
                val gap = if (options.has("gap")) options.valueOf("gap") as Int else null
                require(gap == null || gap >= 0) { "Illegal gap: $gap, expected >= 0" }

                val workingDir = options.valueOf("workdir") as Path
                val id = peaksPath?.stemGz ?: reduceIds(
                    listOfNotNull(
                        modelId,
                        fdr.toString(),
                        sensitivity?.toString(),
                        gap?.toString()
                    )
                )
                var logPath = options.valueOf("log") as Path?
                val chromSizesPath = options.valueOf("chrom.sizes") as Path?

                // Configure working directories
                LOG.info("WORKING DIR: $workingDir")
                if (!OmnipeakCLA.ignoreConfigurePaths) {
                    configurePaths(workingDir, genomesPath = chromSizesPath?.parent, logsPath = logPath?.parent)
                }
                // Configure logging to file
                if (logPath == null) {
                    logPath = Configuration.logsPath / "$id.log"
                }
                Logs.addLoggingToFile(logPath)
                LOG.info("LOG: $logPath")

                // Call now to preserve params logging order
                val lazyResults = logParametersAndPrepareLazyResults(options)

                val bigWig = options.contains("bigwig")
                LOG.info("BIGWIG: $bigWig")


                val multipleTesting = if ("multiple" in params)
                    MultipleTesting.valueOf(options.valueOf("multiple") as String)
                else
                    OMNIPEAK_DEFAULT_MULTIPLE_TEST_CORRECTION
                LOG.info("MULTIPLE TEST CORRECTION: ${multipleTesting.description}")

                val clip = options.valueOf("clip") as Double
                val summits = "summits" in options
                LOG.info("SUMMITS: $summits")
                val fragmentationLight = when {
                    summits || gap != null -> 0.0
                    options.has("f-light") -> options.valueOf("f-light") as Double
                    else -> OMNIPEAK_DEFAULT_FRAGMENTATION_LIGHT
                }
                val fragmentationHard = when {
                    summits || gap != null -> 0.0
                    options.has("f-hard") -> options.valueOf("f-hard") as Double
                    else -> OMNIPEAK_DEFAULT_FRAGMENTATION_HARD
                }

                val fragmentationSpeed = when {
                    summits || gap != null -> 0.0
                    options.has("f-speed") ->
                        options.valueOf("f-speed") as Double

                    else -> OMNIPEAK_DEFAULT_FRAGMENTATION_SPEED
                }

                if (peaksPath != null) {
                    LOG.info("FDR: $fdr")
                    if (sensitivity != null) {
                        LOG.info("SENSITIVITY: $sensitivity")
                    }
                    if (gap != null) {
                        LOG.info("GAP: $gap")
                    }
                    if (gap == null) {
                        LOG.info("FRAGMENTATION MIN THRESHOLD: $fragmentationLight")
                        LOG.info("FRAGMENTATION MAX THRESHOLD: $fragmentationHard")
                        LOG.info("FRAGMENTATION SPEED THRESHOLD: $fragmentationSpeed")
                    }
                    LOG.info("CLIP: $clip")
                    LOG.info("PEAKS: $peaksPath")
                } else {
                    LOG.info("NO peaks path given, process model fitting only.")
                    LOG.info("Fdr, sensitivity, gap, clip options are ignored.")
                }

                val threads = options.valueOf("threads") as Int? ?: Runtime.getRuntime().availableProcessors()
                check(threads > 0) {
                    "Negative threads value: $threads"
                }
                check(threads <= Runtime.getRuntime().availableProcessors()) {
                    "Too big threads value $threads > ${Runtime.getRuntime().availableProcessors()}"
                }
                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")

                // Print all the constants, which are not configured using command line
                if (LOG.isDebugEnabled) {
                    printConstants()
                }

                // Finally get results
                val (actualModelPath, results) = lazyResults.value
                val fitInfo = results.fitInfo
                check(fitInfo is AbstractOmnipeakAnalyzeFitInformation) {
                    "Expected ${OmnipeakAnalyzeFitInformation::class.java.simpleName}, got ${fitInfo::class.java.name}"
                }
                val aboutModel = results.modelInformation(actualModelPath)
                LOG.info(aboutModel.joinToString("\n") { (k, v) ->
                    "${k.name}: ${k.render(v)}"
                })
                LOG.debug(results.model.toString())

                val genomeQuery = fitInfo.genomeQuery()
                val fragment = fitInfo.fragment
                val bin = fitInfo.binSize

                if (bigWig) {
                    if (fitInfo !is OmnipeakAnalyzeFitInformation) {
                        LOG.warn("Bigwig coverage is possible only for analyze command")
                    } else {
                        val bigWigPath = (peaksPath!!.toString() + ".bw").toPath()
                        BigWigCoverageWriter.write(results, genomeQuery, bigWigPath, blackListPath)
                    }
                }

                if (deepAnalysis) {
                    if (fitInfo !is OmnipeakAnalyzeFitInformation) {
                        LOG.warn("Deep analysis is possible only for analyze command")
                    }
                }
                if (peaksPath != null) {
                    val peaks =
                        OmnipeakModelToPeaks.getPeaks(
                            results, genomeQuery, fdr, multipleTesting,
                            sensitivity, gap, summits,
                            fragmentationLight, fragmentationHard, fragmentationSpeed,
                            clip = clip,
                            blackListPath = blackListPath,
                            name = peaksPath.fileName.stem,
                        )
                    val peaksList = processBlackList(genomeQuery, peaks, blackListPath)

                    LOG.info("${peaksPath.fileName.stem} format chromosome, start, end, name, score, strand, signal, -log(p), -log(q)")
                    Peak.savePeaks(
                        peaksList, peaksPath,
                        "peak${if (fragment is FixedFragment) "_$fragment" else ""}_" +
                                "${bin}_${fdr}_${peaks.sensitivity}_${peaks.gap}"
                    )
                    LOG.info("Peaks saved to $peaksPath")
                    val aboutPeaks = PeaksInfo.compute(
                        genomeQuery,
                        peaksList.map { it.location }.stream(),
                        peaksPath.toUri(),
                        fitInfo.paths.map { it.treatment }
                    )
                    LOG.info(aboutPeaks.joinToString("\n") { (k, v) ->
                        "${k.name}: ${k.render(v)}"
                    })

                    if (deepAnalysis && fitInfo is OmnipeakAnalyzeFitInformation) {
                        fitInfo.prepareData()
                        doDeepAnalysis(
                            actualModelPath,
                            results,
                            fitInfo,
                            genomeQuery,
                            fdr,
                            sensitivity, gap,
                            fragmentationLight, fragmentationHard,
                            fragmentationSpeed,
                            blackListPath,
                            peaksList,
                            peaksPath
                        )
                    }
                }
                if (modelPath == null && !keepCacheFiles) {
                    LOG.debug("Clean coverage caches")
                    fitInfo.cleanCaches()
                }
            }
        }
    }

    fun processBlackList(
        genomeQuery: GenomeQuery,
        peaks: OmnipeakResult,
        blacklistPath: Path?
    ): List<Peak> {
        var peaksList = peaks.toList()
        if (blacklistPath != null) {
            LOG.info("Filter out blacklisted regions")
            val blackList = LocationsMergingList.load(genomeQuery, blacklistPath)
            peaksList = peaksList.filter { !blackList.intersects(it.location) }
        }
        return peaksList
    }


    /**
     * Retrieves the paths (treatment, optional control)
     * either from command-line options or from the stored fit information.
     *
     * If both are available, checks that they are consistent.
     */
    internal fun prepareAndCheckTreatmentControlPaths(
        options: OptionSet, fitInformation: AbstractOmnipeakAnalyzeFitInformation? = null, log: Boolean = false
    ): List<OmnipeakDataPaths> {
        val commandLineTreatmentPaths = options.valuesOf("treatment") as List<Path>
        val commandLineControlPaths = options.valuesOf("control") as List<Path>

        var paths = OmnipeakCLA.matchTreatmentsAndControls(
            commandLineTreatmentPaths, commandLineControlPaths
        )

        val fitInfoPaths = fitInformation?.paths
        if (fitInfoPaths != null) {
            if (paths != null) {
                check(paths == fitInfoPaths) {
                    "Stored treatment-control pairs ${fitInfoPaths.joinToString()} differ from the ones inferred " +
                            "from the command line arguments: ${paths!!.joinToString()}"
                }
            } else {
                paths = fitInfoPaths
            }
        }

        checkNotNull(paths) {
            "No treatment files and no existing model file provided, exiting."
        }

        if (log) {
            LOG.info("TREATMENT: ${paths.map { it.treatment }.joinToString(", ", transform = Path::toString)}")
            paths.mapNotNull { it.control }.let {
                if (it.isNotEmpty()) {
                    LOG.info("CONTROL: ${it.joinToString(", ", transform = Path::toString)}")
                } else {
                    LOG.info("CONTROL: none")
                }
            }
        }
        return paths
    }

    /**
     * Log parameters and optionally computation of Omnipeak model, represented by [OmnipeakFitResults].
     * Doesn't fit the model; this happens only if .value is invoked
     */
    private fun logParametersAndPrepareLazyResults(options: OptionSet): Lazy<Pair<Path, OmnipeakFitResults>> {
        val modelPath = options.valueOf("model") as Path?
        if (modelPath != null) {
            LOG.info("MODEL: $modelPath")
        }
        if (modelPath != null && modelPath.exists && modelPath.size.isNotEmpty()) {
            LOG.debug(
                "Model file {} exists and is not empty, Omnipeak will use it to substitute the missing " +
                        "command line arguments and verify the provided ones.",
                modelPath
            )
            val results = OmnipeakModelFitExperiment.loadResults(modelPath = modelPath)
            check(results.fitInfo is AbstractOmnipeakAnalyzeFitInformation) {
                "Expected ${OmnipeakAnalyzeFitInformation::class.java.simpleName}, got ${results.fitInfo::class.java.name}"
            }
            prepareAndCheckTreatmentControlPaths(options, results.fitInfo, log = true)
            val blacklistPath = options.valueOf("blacklist") as Path?
            if (blacklistPath != null) {
                LOG.info("BLACKLIST FILE: $blacklistPath")
            }
            val workingDir = options.valueOf("workdir") as Path
            LOG.info("WORKING DIR: $workingDir")
            val chromSizesPath = options.valueOf("chrom.sizes") as Path?
            LOG.info("CHROM.SIZES: $chromSizesPath")
            val chromosomesToProcess = if (options.has("chromosomes"))
                options.valueOf("chromosomes").toString().split(',', ' ')
            else
                null
            if (chromosomesToProcess != null) {
                LOG.info("CHROMOSOMES: ${chromosomesToProcess.joinToString(", ")}")
            }
            if (chromSizesPath != null) {
                checkGenomeInFitInformation(chromSizesPath, results.fitInfo)
            }
            OmnipeakCLA.getBin(options, results.fitInfo, log = true)
            OmnipeakCLA.getFragment(options, results.fitInfo, log = true)
            OmnipeakCLA.getUnique(options, results.fitInfo, log = true)
            val keepCacheFiles = "keep-cache" in options
            LOG.info("KEEP-CACHE: $keepCacheFiles")
            if (!keepCacheFiles) {
                LOG.warn("Keep cache files setting (false) is discarded, since fitted model already exists.")
            }
            return lazyOf(modelPath to results)        // Create fake lazy of already computed results
        } else {
            val paths = prepareAndCheckTreatmentControlPaths(options, log = true)
            val blacklistPath = options.valueOf("blacklist") as Path?
            if (blacklistPath != null) {
                LOG.info("BLACKLIST FILE: $blacklistPath")
            }
            val workingDir = options.valueOf("workdir") as Path
            LOG.info("WORKING DIR: $workingDir")
            val chromSizesPath = options.valueOf("chrom.sizes") as Path?
            LOG.info("CHROM.SIZES: $chromSizesPath")
            val chromosomesToProcess = if (options.has("chromosomes"))
                options.valueOf("chromosomes").toString().split(',', ' ')
            else
                null
            if (chromosomesToProcess != null) {
                LOG.info("CHROMOSOMES: ${chromosomesToProcess?.joinToString(", ")}")
            }
            val explicitFormat: InputFormat? = OmnipeakCLA.readsFormat(options, log = true)
            val fragment = OmnipeakCLA.getFragment(options, log = true)
            val unique = OmnipeakCLA.getUnique(options, log = true)
            val bin = OmnipeakCLA.getBin(options, log = true)
            val fitThreshold = OmnipeakCLA.getFitThreshold(options, log = true)
            val fitMaxIterations = OmnipeakCLA.getFitMaxIteration(options, log = true)
            val hmmEstimateSNR = options.valueOf("hmm-snr") as Double
            LOG.info("HMM ESTIMATE SNR: $hmmEstimateSNR")
            val hmmLow = options.valueOf("hmm-low") as Double
            LOG.info("HMM LOW STATE MIN: $hmmLow")
            val saveExtendedInfo = options.has("ext")
            LOG.info("EXTENDED MODEL INFO: $saveExtendedInfo")
            val keepCacheFiles = "keep-cache" in options
            LOG.info("KEEP-CACHE: $keepCacheFiles")
            return lazy {
                val genomeQuery = if (chromosomesToProcess != null)
                    GenomeQuery(Genome[chromSizesPath!!], *chromosomesToProcess.toTypedArray())
                else
                    GenomeQuery(Genome[chromSizesPath!!])
                val experiment = OmnipeakPeakCallingExperiment.getExperiment(
                    genomeQuery,
                    paths, explicitFormat,
                    fragment, unique, bin,
                    hmmEstimateSNR, hmmLow,
                    modelPath, fitThreshold, fitMaxIterations,
                    saveExtendedInfo, keepCacheFiles
                )
                experiment.modelPath to experiment.results
            }
        }
    }

}
