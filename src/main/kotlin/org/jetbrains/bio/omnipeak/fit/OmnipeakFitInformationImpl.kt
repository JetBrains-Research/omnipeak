package org.jetbrains.bio.omnipeak.fit

import com.google.common.math.IntMath
import com.google.gson.*
import com.google.gson.reflect.TypeToken
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.toPath
import java.lang.reflect.Type
import java.math.RoundingMode
import java.nio.file.Path
import java.util.*

fun OmnipeakFitInformation.checkGenome(genome: Genome) {
    check(this.build == genome.build) {
        "Wrong genome build, expected: ${this.build}, got: ${genome.build}"
    }
}

fun OmnipeakFitInformation.checkChromosome(chromosome: Chromosome) {
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
fun OmnipeakFitInformation.offsetsMap(): IntArray =
    (listOf(0) + chromosomesSizes.keys.sorted().map {
        IntMath.divide(chromosomesSizes[it]!!, binSize, RoundingMode.CEILING)
    }).toIntArray().let {
        Arrays.parallelPrefix(it) { a, b -> a + b }; it
    }

/**
 * Creates binned offsets for [chromosome] using [binSize].
 */
fun OmnipeakFitInformation.offsets(chromosome: Chromosome): IntArray {
    checkGenome(chromosome.genome)
    checkChromosome(chromosome)
    return chromosome.range.slice(binSize).mapToInt { it.startOffset }.toArray()
}

/**
 * Returns boundaries of the squashed dataframe region corresponding to the given [chromosome].
 *
 * See also: [merge] and [split].
 */
fun OmnipeakFitInformation.getChromosomesIndices(chromosome: Chromosome): Pair<Int, Int> {
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
fun OmnipeakFitInformation.getChromosomesIndices(chromosome: String): Pair<Int, Int> {
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
 * row number equal to the number of bins on the appropriate chromosome for each chromosome.
 */
fun OmnipeakFitInformation.merge(statesDataFrame: Map<String, DataFrame>): DataFrame {
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
fun OmnipeakFitInformation.split(dataFrame: DataFrame, genomeQuery: GenomeQuery? = null): Map<String, DataFrame> {
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

/**
 * Save the [OmnipeakFitInformation] object at the given path as JSON.
 *
 * Inverse of [load].
 */
fun OmnipeakFitInformation.save(path: Path) {
    path.parent.createDirectories()
    path.bufferedWriter().use { OmnipeakFitInformationGSON.GSON.toJson(this, it) }
}

/**
 * Generate [OmnipeakFitInformation.chromosomesSizes] instance from a [GenomeQuery]
 */
fun OmnipeakFitInformation.Companion.chromSizes(genomeQuery: GenomeQuery): LinkedHashMap<String, Int> =
    LinkedHashMap<String, Int>().apply {
        genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
    }

/**
 * Loads a [OmnipeakFitInformation] instance from a JSON file.
 *
 * Inverse of [OmnipeakFitInformation.save]. Since "save" stores the fully-qualified class name,
 * "load" instantiates the correct class. If this class cannot be cast to [T],
 * [IllegalStateException] is thrown.
 */
inline fun <reified T : OmnipeakFitInformation> OmnipeakFitInformation.Companion.load(json: String): T? {
    val info = OmnipeakFitInformationGSON.GSON.fromJson(json, OmnipeakFitInformation::class.java)
    if (info != null) {
        check(info is T) {
            "Information loaded from is not of expected type, expected: ${T::class.java.name}, got: ${info::class.java.name}"
        }
    }
    return info as T?
}

object OmnipeakFitInformationGSON {

    private class FragmentTypeAdapter : JsonSerializer<Fragment>, JsonDeserializer<Fragment> {

        override fun serialize(
            src: Fragment,
            typeOfSrc: Type,
            context: JsonSerializationContext
        ) = JsonPrimitive(src.toString())

        override fun deserialize(
            json: JsonElement,
            typeOfT: Type,
            context: JsonDeserializationContext
        ) = Fragment.fromString(json.asJsonPrimitive.asString)
    }

    private class PathTypeAdapter : JsonSerializer<Path>, JsonDeserializer<Path> {

        override fun serialize(
            src: Path,
            typeOfSrc: Type,
            context: JsonSerializationContext
        ) = JsonPrimitive(src.toAbsolutePath().toString())

        override fun deserialize(
            json: JsonElement,
            typeOfT: Type,
            context: JsonDeserializationContext
        ) = json.asJsonPrimitive.asString.toPath()
    }

    val GSON: Gson = GsonBuilder()
        .setFieldNamingStrategy(GSONUtil.NO_MY_UNDERSCORE_NAMING_STRATEGY)
        .registerTypeAdapterFactory(
            GSONUtil.classSpecificFactory(OmnipeakFitInformation::class.java) { gson, factory ->
                GSONUtil.classAndVersionAdapter(
                    gson, factory, "fit.information.fqn", "version"
                )
            }
        )
        .registerTypeAdapter(
            object : TypeToken<Fragment>() {}.type,
            FragmentTypeAdapter()
        )
        .registerTypeHierarchyAdapter(
            Path::class.java,
            PathTypeAdapter()
        )
        .setPrettyPrinting()
        .create()
}
