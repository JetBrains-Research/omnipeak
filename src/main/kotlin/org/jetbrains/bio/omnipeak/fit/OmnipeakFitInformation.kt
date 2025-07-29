package org.jetbrains.bio.omnipeak.fit

import com.google.common.math.IntMath
import com.google.gson.*
import com.google.gson.reflect.TypeToken
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.omnipeak.coverage.BinnedCoverageQuery
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation.Companion.chromSizes
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation.Companion.load
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.toPath
import java.lang.reflect.Type
import java.math.RoundingMode
import java.nio.file.Path
import java.util.*

/**
 * The most common interface for all fit information classes.
 *
 * [OmnipeakFitInformation] instance is designed to contain all information necessary to uniquely identify the input
 * of a Omnipeak-like model fitting experiment. For example, [OmnipeakAnalyzeFitInformation] completely describes
 * the input of the classical `omnipeak analyze` command.
 *
 * [OmnipeakFitInformation] object is a part of [OmnipeakFitResults], and its type is type parameter
 * of [OmnipeakModelFitExperiment].
 *
 * All Omnipeak-like experiments produce a single squashed float array of log null probabilities ("null.npz").
 * This interface contains methods to squash ([merge]) and unsquash ([split]) the chromosome-wise dataframes.
 * It can also generate bin start [offsets] for a single chromosome.
 */
interface OmnipeakFitInformation {

    /**
     * A unique string identifier (include some kind of object hash if you compress identifiers). It's used
     * to generate the model file name if it's not provided. [reduceIds] is a recommended way to implement this property.
     */
    val id: String

    /** Genome build. */
    val build: String

    /** Bin size in base pairs. */
    val binSize: Int

    /** A map of chromosome name -> chromosome length entries. */
    val chromosomesSizes: LinkedHashMap<String, Int>

    /** A query that returns a dataframe for each chromosome to serve as model input. */
    val dataQuery: Query<Chromosome, DataFrame>


    /**
     * Prepares data queries for [score] function.
     */
    fun prepareData()

    fun scoresAvailable(): Boolean

    /**
     * Computes range score, either coverage for analyze experiment or log2 fold change for difference.
     * Call [prepareData] beforehand!
     * Can be slow in multi-thread usage, due to synchronized lazy on NormalizedCoverageQuery#treatmentReads
     * Consider replacing to direct calls to ReadsCoverage
     */
    fun score(chromosomeRange: ChromosomeRange): Double

    fun isControlAvailable(): Boolean

    /**
     * Returns scaled control score or null if not available.
     * Call [prepareData] beforehand!
     * Can be slow in multi-thread usage, due to synchronized lazy on NormalizedCoverageQuery#treatmentReads
     * Consider replacing to direct calls to ReadsCoverage
     */
    fun controlScore(chromosomeRange: ChromosomeRange): Double

    fun cleanCaches()

    fun genomeQuery(): GenomeQuery =
        GenomeQuery(Genome[build, chromosomesSizes], *chromosomesSizes.keys.toTypedArray())

    fun checkGenome(genome: Genome) {
        check(this.build == genome.build) {
            "Wrong genome build, expected: ${this.build}, got: ${genome.build}"
        }
    }

    private fun checkChromosome(chromosome: Chromosome) {
        check(containsChromosomeInfo(chromosome)) {
            "Missing chromosome in ${chromosomesSizes.keys.toList()}: ${chromosome.name}"
        }
        check(chromosome.length == chromosomesSizes[chromosome.name]) {
            "Wrong chromosome ${chromosome.name} size, " +
                    "expected ${chromosomesSizes[chromosome.name]}, got: ${chromosome.length}"
        }
    }

    /**
     * Offsets map is used to store all chromosome information into single array.
     */
    private fun offsetsMap(): IntArray =
        (listOf(0) + chromosomesSizes.keys.sorted().map {
            IntMath.divide(chromosomesSizes[it]!!, binSize, RoundingMode.CEILING)
        }).toIntArray().let {
            Arrays.parallelPrefix(it) { a, b -> a + b }; it
        }


    /**
     * Creates binned offsets for [chromosome] using [binSize].
     */
    fun offsets(chromosome: Chromosome): IntArray {
        checkGenome(chromosome.genome)
        checkChromosome(chromosome)
        return chromosome.range.slice(binSize).mapToInt { it.startOffset }.toArray()
    }

    /**
     * Returns boundaries of the squashed dataframe region corresponding to the given [chromosome].
     *
     * See also: [merge] and [split].
     */
    fun getChromosomesIndices(chromosome: Chromosome): Pair<Int, Int> {
        checkChromosome(chromosome)
        val offsetsMap = offsetsMap()
        val index = chromosomesSizes.keys.sorted().indexOf(chromosome.name)
        return offsetsMap[index] to offsetsMap[index + 1]
    }

    /**
     * Returns boundaries of the squashed dataframe region corresponding to the given [chromosome].
     *
     * See also: [merge] and [split].
     */
    fun getChromosomesIndices(chromosome: String): Pair<Int, Int> {
        check(chromosome in chromosomesSizes) {
            "Missing chromosome in ${chromosomesSizes.keys.toList()}: $chromosome"
        }
        val offsetsMap = offsetsMap()
        val index = chromosomesSizes.keys.sorted().indexOf(chromosome)
        return offsetsMap[index] to offsetsMap[index + 1]
    }

    /**
     * Merges (row-binds) the chromosome-wise dataframes in an unambiguous way.
     *
     * Inverse of [split].
     *
     * @param statesDataFrame a map of chromosome name -> dataframe entries. Must contain a dataframe with
     * row number equal to the number of bins on the appropriate chromosome for each chromosome in [chromSizes].
     */
    fun merge(statesDataFrame: Map<String, DataFrame>): DataFrame {
        return DataFrame.rowBind(chromosomesSizes.keys.sorted().map { statesDataFrame[it]!! }.toTypedArray())
    }

    /**
     * Splits the squashed dataframe into chromosome-wise parts.
     *
     * Inverse of [merge].
     *
     * @param genomeQuery Optional smaller [GenomeQuery]. If provided, only requested chromosome-wise dataframes
     * are returned.
     */
    fun split(dataFrame: DataFrame, genomeQuery: GenomeQuery?): Map<String, DataFrame> {
        return if (genomeQuery != null) {
            checkGenome(genomeQuery.genome)
            genomeQuery.get()
                .filter { containsChromosomeInfo(it) }.associate { chromosome ->
                    val (start, end) = getChromosomesIndices(chromosome)
                    chromosome.name to dataFrame.iloc[start until end]
                }
        } else {
            chromosomesSizes.keys.associateWith { chromosome ->
                val (start, end) = getChromosomesIndices(chromosome)
                dataFrame.iloc[start until end]
            }
        }
    }

    fun containsChromosomeInfo(chromosome: Chromosome) = chromosome.name in chromosomesSizes

    /**
     * Save the [OmnipeakFitInformation] object at the given path as JSON.
     *
     * Inverse of [load].
     */
    fun save(path: Path) {
        path.parent.createDirectories()
        path.bufferedWriter().use { GSON.toJson(this, it) }
    }

    fun difference(loadedFitInfo: OmnipeakFitInformation): String?

    companion object {

        /**
         * Generate [OmnipeakFitInformation.chromosomesSizes] instance from a [GenomeQuery]
         */
        fun chromSizes(genomeQuery: GenomeQuery) =
            LinkedHashMap<String, Int>().apply {
                genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
            }

        object FragmentTypeAdapter : JsonSerializer<Fragment>, JsonDeserializer<Fragment> {

            override fun serialize(
                src: Fragment, typeOfSrc: Type,
                context: JsonSerializationContext
            ): JsonElement = context.serialize(src.toString())

            override fun deserialize(
                json: JsonElement, typeOfT: Type,
                context: JsonDeserializationContext
            ): Fragment {
                val str = context.deserialize<String>(json, object : TypeToken<String>() {}.type)
                try {
                    return Fragment.fromString(str)
                } catch (e: NumberFormatException) {
                    throw IllegalStateException("Failed to deserialize $str", e)
                }
            }
        }

        object PathTypeAdapter : JsonSerializer<Path>, JsonDeserializer<Path> {

            override fun serialize(
                src: Path, typeOfSrc: Type,
                context: JsonSerializationContext
            ): JsonElement = JsonPrimitive(src.toString())

            override fun deserialize(
                json: JsonElement, typeOfT: Type,
                context: JsonDeserializationContext
            ): Path {
                try {
                    return json.asString.toPath()
                } catch (e: NumberFormatException) {
                    throw IllegalStateException("Failed to deserialize ${json.asString}", e)
                }
            }
        }

        private val GSON: Gson = GsonBuilder()
            .setPrettyPrinting()
            .setFieldNamingStrategy(GSONUtil.NO_MY_UNDERSCORE_NAMING_STRATEGY)
            .registerTypeAdapterFactory(
                GSONUtil.classSpecificFactory(OmnipeakFitInformation::class.java) { gson, factory ->
                    GSONUtil.classAndVersionAdapter(gson, factory, "fit.information.fqn", "version")
                })
            .registerTypeAdapter(object : TypeToken<Fragment>() {}.type, FragmentTypeAdapter)
            .registerTypeHierarchyAdapter(Path::class.java, PathTypeAdapter)
            .create()

        /**
         * Loads a [OmnipeakFitInformation] instance from a JSON file.
         *
         * Inverse of [OmnipeakFitInformation.save]. Since "save" stores the fully-qualified class name,
         * "load" instantiates the correct class. If this class cannot be cast to [T],
         * [IllegalStateException] is thrown.
         */
        @Suppress("unchecked_cast")
        fun <T : OmnipeakFitInformation> load(json: String): T? {
            return GSON.fromJson(json, OmnipeakAnalyzeFitInformation::class.java) as T?
        }

    }
}

/**
 * Calculates the binned coverage DataFrame for a list of normalized coverage queries.
 * Using clip percentile may slightly help to avoid out-of-range during model fit.
 */
fun List<BinnedCoverageQuery>.binnedCoverageDataFrame(
    chromosome: Chromosome,
    labels: Array<String>,
): DataFrame {
    var res = DataFrame()
    forEachIndexed { d, inputQuery ->
        val binnedCoverage = inputQuery.apply(chromosome)
        res = res.with(labels[d], binnedCoverage)
    }
    return res
}