package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiment.Experiment
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.omnipeak.SPAN2
import org.jetbrains.bio.omnipeak.SPAN2.toOmnipeak
import org.jetbrains.bio.omnipeak.fit.OmnipeakModelFitExperiment.Companion.loadResults
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.toFloatArray
import org.jetbrains.bio.util.*
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.Callable


/**
 * A generic class for Omnipeak a tool for analyzing and comparing ChIP-Seq data.
 *
 * The end result of the experiment is the [results] property. It's lazy (won't do any calculation until
 * actually accessed) and cached (if possible, will be loaded from the previously created file, if not, will
 * be saved to a file after the computation). The results are saved in a TAR file.
 *
 * @param fixedModelPath    If not null, the experiment will use this path for saving/loading the results.
 *                          Otherwise, [defaultModelPath] will be used (it usually depends on [fitInformation] id).
 *
 */
abstract class OmnipeakModelFitExperiment<
        out Model : ClassificationModel, out FitInfo : OmnipeakFitInformation, State : Any
        > protected constructor(
    val fitInformation: FitInfo,
    private val modelFitter: Fitter<Model>,
    private val modelClass: Class<out Model>,
    private val availableStates: Array<State>,
    private val nullHypothesis: NullHypothesis<State>,
    private val threshold: Double,
    private val maxIterations: Int,
    private val fixedModelPath: Path?,
    private val saveExtendedInfo: Boolean,
    private val keepModelFile: Boolean
) : Experiment(null) {

    val genomeQuery = fitInformation.genomeQuery()
    private val dataQuery = fitInformation.dataQuery

    /**
     * Preprocessed data by chromosomes, chromosomes are sorted by name.
     */
    private val preprocessedData: List<Preprocessed<DataFrame>> by lazy {
        genomeQuery.get().sortedBy { it.name }.map { Preprocessed.of(dataQuery.apply(it)) }
    }

    val results: OmnipeakFitResults by lazy { getOrLoadResults() }

    fun getStatesDataFrame(chromosome: Chromosome): DataFrame = sliceStatesDataFrame(statesDataFrame, chromosome)

    override fun doCalculations() {
        results.logNullMemberships
    }

    /**
     * A unique path used to store/load results if [fixedModelPath] is null.
     *
     * This property is normally implemented through [fitInformation] id.
     */
    abstract val defaultModelPath: Path

    /**
     * We use "get" because we need the actual value of [defaultModelPath] implemented in the descendant class.
     */
    val modelPath get() = fixedModelPath ?: defaultModelPath

    private val saveModel get() = fixedModelPath != null || keepModelFile

    private fun calculateModel(): Model {
        return modelFitter.fit(
            preprocessedData,
            title = dataQuery.id,
            threshold = threshold,
            maxIterations = maxIterations
        )
    }

    private fun calculateStatesDataFrame(model: Model): DataFrame {
        val empty = DataFrame()
        val dataFrames = Array(preprocessedData.size) { empty }
        val tasks = preprocessedData.mapIndexed { index, preprocessed ->
            Callable {
                val logMemberships = model.evaluate(preprocessed)
                var df = DataFrame()
                availableStates.forEachIndexed { j, state ->
                    val f64Array = logMemberships.V[j]
                    // Convert [Double] to [Float] to save space, see #1163
                    df = df.with(state.toString(), f64Array.toFloatArray())
                }
                df = df.with(
                    "state", model.predict(preprocessed)
                        .map { availableStates[it].toString() }.toTypedArray()
                )
                dataFrames[index] = df
            }
        }
        tasks.await(parallel = true)
        return DataFrame.rowBind(dataFrames)
    }

    private val statesDataFrame: DataFrame by lazy {
        @Suppress("UNCHECKED_CAST")
        calculateStatesDataFrame(results.model as Model)
    }

    /**
     * The [statesDataFrame] contains aggregated information for the whole genome,
     * this method returns data chunk for given [chromosome].
     * @return dataframe slice for given chromosome.
     */
    private fun sliceStatesDataFrame(statesDataFrame: DataFrame, chromosome: Chromosome): DataFrame {
        val (start, end) = fitInformation.getChromosomesIndices(chromosome)
        return statesDataFrame.iloc[start until end]
    }

    /**
     * Return map state -> f64 array of log probabilities for each position to be in given hidden state.
     * Membership = Probability here.
     */
    private fun getLogMemberships(chromosomeStatesDF: DataFrame): Map<State, F64Array> =
        availableStates.associateBy({ it }) { chromosomeStatesDF.f64Array(it.toString()) }

    /**
     * Compute and save [OmnipeakFitResults], i.e. fit information, trained model and null hypothesis probabilities.
     * If already processed, load them [loadResults].
     *
     * IMPORTANT!
     * We take care not to access any of the lazy properties here, since they depend on this method for initialization.
     */
    private fun getOrLoadResults(): OmnipeakFitResults {
        var computedResults: OmnipeakFitResults? = null
        modelPath.checkOrRecalculate("Model fit", ignoreEmptyFile = !saveModel) { (p) ->
            withTempDirectory(modelPath.stem) { dir ->
                val modelPath = dir / MODEL_JSON
                LOG.info("Computing data model...")
                val model = calculateModel()
                LOG.info("Done computing data model")
                if (saveModel) {
                    LOG.info("Saving model information")
                    model.save(modelPath)
                    LOG.debug("Model saved to {}", modelPath)
                } else {
                    LOG.debug("Model is not saved")
                }
                val informationPath = dir / INFORMATION_JSON
                if (saveModel) {
                    fitInformation.save(informationPath)
                    LOG.debug("Fit information saved to {}", informationPath)
                } else {
                    LOG.debug("Fit information is not saved")
                }

                LOG.info("Analyzing model states")
                val statesDataFrame = calculateStatesDataFrame(model)
                LOG.debug("Computing probabilities")
                val chromosomeToDataFrameMap = genomeMap(genomeQuery, parallel = true) {
                    val logMemberships = getLogMemberships(sliceStatesDataFrame(statesDataFrame, it))
                    val logNullMemberships = nullHypothesis.apply(logMemberships)
                    // Convert [Double] to [Float] to save space, see #1163
                    DataFrame().with(NULL, logNullMemberships.toFloatArray())
                }
                val logNullMembershipsDF = fitInformation.merge(
                    genomeQuery.get().associate { it.name to chromosomeToDataFrameMap[it] }
                )
                LOG.debug("Done null hypothesis log memberships")

                val logNullMembershipsPath = dir / NULL_NPZ
                logNullMembershipsDF.save(logNullMembershipsPath)
                LOG.debug("LogNullMemberships saved to {}", logNullMembershipsPath)
                var statesDataFrameMap: Map<String, DataFrame>? = null

                if (saveModel) {
                    if (saveExtendedInfo) {
                        LOG.debug("Saving full states dataframe")
                        val statesDataFramePath = dir / "states.npz"
                        statesDataFrame.save(statesDataFramePath)
                        LOG.debug("States saved to {}", statesDataFramePath)
                        Tar.compress(
                            p,
                            *listOf(
                                modelPath,
                                informationPath,
                                logNullMembershipsPath,
                                statesDataFramePath
                            ).map(Path::toFile)
                                .toTypedArray()
                        )
                        statesDataFrameMap = fitInformation.split(statesDataFrame, genomeQuery)
                    } else {
                        Tar.compress(
                            p,
                            *listOf(modelPath, informationPath, logNullMembershipsPath).map(Path::toFile)
                                .toTypedArray()
                        )
                    }
                }

                val logNullMembershipsMap = fitInformation.split(logNullMembershipsDF, genomeQuery)
                computedResults = OmnipeakFitResults(fitInformation, model, logNullMembershipsMap, statesDataFrameMap)
            }
        }
        return if (computedResults != null) {
            if (saveModel) {
                LOG.info("Model saved: $modelPath")
            } else {
                // Clean empty file
                modelPath.deleteIfExists()
                LOG.info("Model is not saved")
            }
            computedResults!!
        } else {
            val loadedResults = loadResults(genomeQuery, modelPath)
            val loadedFitInfo = loadedResults.fitInfo
            check(loadedFitInfo == fitInformation) {
                "Different model information!\n${loadedFitInfo.difference(fitInformation)}"
            }
            loadedResults
        }
    }

    override fun toString(): String {
        return """
            OmnipeakModelFitExperiment:
                fitInformation: $fitInformation
                modelPath: $modelPath
                threshold: $threshold
                maxIterations: $maxIterations
                saveModel: $saveModel
                saveExtendedInfo: $saveExtendedInfo
        """.trimIndent()
    }

    companion object {
        private const val INFORMATION_JSON = "information.json"
        private const val MODEL_JSON = "model.json"
        private const val NULL_NPZ = "null.npz"
        const val NULL = "null"

        val LOG: Logger = LoggerFactory.getLogger(OmnipeakModelFitExperiment::class.java)

        /**
         * Retain only the chromosomes for which at least one treatment file has at least one read on them.
         */
        fun filterGenomeQueryWithData(
            genomeQuery: GenomeQuery,
            paths: List<OmnipeakDataPaths>,
            explicitFormat: ReadsFormat?,
            fragment: Fragment,
            unique: Boolean = true
        ): GenomeQuery {
            val chromosomes = genomeQuery.get()
            val nonEmptyChromosomes = hashSetOf<Chromosome>()
            paths.forEach { (t, c) ->
                val coverage =
                    ReadsQuery(genomeQuery, t, explicitFormat, unique, fragment, showLibraryInfo = false).get()
                if (c != null) {
                    // we have to be sure that the control coverage cache is calculated for the full genome query,
                    // otherwise we can get some very hard-to-catch bugs later
                    ReadsQuery(genomeQuery, c, explicitFormat, unique, fragment, showLibraryInfo = false).get()
                }
                nonEmptyChromosomes.addAll(
                    chromosomes.filter { coverage.getBothStrandsCoverage(it.chromosomeRange) > 0 }
                )
            }

            if (nonEmptyChromosomes.isEmpty()) {
                val errMessage = "Model can't be trained on empty coverage, exiting."
                LOG.error(errMessage)
                throw IllegalStateException(errMessage)
            }

            val emptyChromosomes = chromosomes.filter { it !in nonEmptyChromosomes }
            if (emptyChromosomes.isNotEmpty()) {
                LOG.info("Chromosomes with no reads detected are ignored. Use --debug for details.")
                LOG.debug("Ignored chromosomes: ${emptyChromosomes.joinToString(",") { it.name }}")
            }

            return GenomeQuery(genomeQuery.genome, *nonEmptyChromosomes.map { it.name }.toTypedArray())
        }


        fun loadResults(
            genomeQuery: GenomeQuery? = null,
            modelPath: Path
        ): OmnipeakFitResults {
            LOG.info("Loading model: $modelPath")
            return withTempDirectory(modelPath.stem) { dir ->
                LOG.debug("Started model file decompress: ${modelPath.stem}")
                Tar.decompress(modelPath, dir.toFile())

                LOG.debug("Completed model file decompress and started loading: ${modelPath.stem}")
                val info = OmnipeakFitInformation.load<OmnipeakFitInformation>(
                    toOmnipeak((dir / INFORMATION_JSON).bufferedReader().readText())
                )
                checkNotNull(info) { "Failed to load model information" }
                // Check genome build
                genomeQuery?.let { info.checkGenome(it.genome) }
                // Load model and PEPs
                val model = ClassificationModel.load<ClassificationModel>(
                    toOmnipeak((dir / MODEL_JSON).bufferedReader().readText())
                )
                checkNotNull(model) { "Failed to load model from file" }
                val logNullMembershipsDF = DataFrame.load(dir / NULL_NPZ)
                val logNullMembershipsMap = info.split(logNullMembershipsDF, genomeQuery)
                var statesDfMap: Map<String, DataFrame>? = null
                val statesPath = dir / "states.npz"
                LOG.info("Loading states data frame ${modelPath.stem}")
                if (statesPath.exists) {
                    statesDfMap = info.split(DataFrame.load(statesPath), genomeQuery)
                }
                LOG.info("Completed loading model: ${modelPath.stem}")
                return@withTempDirectory OmnipeakFitResults(info, model, logNullMembershipsMap, statesDfMap)
            }
        }
    }
}

