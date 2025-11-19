package org.jetbrains.bio.omnipeak.fit

import com.google.common.math.DoubleMath
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.omnipeak.InputFormat
import org.jetbrains.bio.omnipeak.coverage.BinnedCoverageQuery
import org.jetbrains.bio.util.reduceIds

data class OmnipeakCompareFitInformation(
    override val build: String,
    val data1: List<OmnipeakDataPaths>,
    val data2: List<OmnipeakDataPaths>,
    val labels1: List<String>,
    val labels2: List<String>,
    override val explicitFormat: InputFormat?,
    override val fragment: Fragment,
    override val unique: Boolean,
    override val binSize: Int,
    override val regressControl: Boolean,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : OmnipeakFitInformation {

    override val id
        get() = OmnipeakAnalyzeFitInformation.generateId(data1 + data2, fragment, binSize, unique, regressControl)

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            check(scoresAvailable()) { "Scores are not available!" }
            val queries = binnedCoverageQueries1!! + binnedCoverageQueries2!!
            val labels = labels1 + labels2
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return queries.binnedCoverageDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(
                        queries.zip(labels).flatMap { (s, l) -> listOf(s.id, l) }
                                + listOf(binSize.toString())
                    )
            }
        }

    @Transient
    var binnedCoverageQueries1: List<BinnedCoverageQuery>? = null

    @Transient
    var binnedCoverageQueries2: List<BinnedCoverageQuery>? = null

    override fun prepareData() {
        if (binnedCoverageQueries1 == null) {
            binnedCoverageQueries1 = data1.map {
                BinnedCoverageQuery.create(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    explicitFormat,
                    fragment,
                    unique,
                    binSize,
                    regressControl,
                    forceCaching = true,
                    showLibraryInfo = false
                )
            }
        }
        if (binnedCoverageQueries2 == null) {
            binnedCoverageQueries2 = data2.map {
                BinnedCoverageQuery.create(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    explicitFormat,
                    fragment,
                    unique,
                    binSize,
                    regressControl,
                    forceCaching = true,
                    showLibraryInfo = false
                )
            }
        }
    }

    override fun scoresAvailable(): Boolean {
        prepareData()
        return binnedCoverageQueries1 != null && binnedCoverageQueries2 != null
    }

    override fun isControlAvailable(): Boolean {
        return scoresAvailable()
    }

    /**
     * Return log2 fold change of average summary coverage across data
     */
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(scoresAvailable()) { "Scores are not available!" }
        return if (binnedCoverageQueries1!!.all { it.areCachesPresent() } &&
            binnedCoverageQueries2!!.all { it.areCachesPresent() }) {
            val score1 = binnedCoverageQueries1!!.sumOf { it.controlNormalizedScore(chromosomeRange) }
                .toDouble() / binnedCoverageQueries1!!.size
            val score2 = binnedCoverageQueries2!!.sumOf { it.controlNormalizedScore(chromosomeRange) }
                .toDouble() / binnedCoverageQueries2!!.size
            if (score2 != 0.0) DoubleMath.log2(score1) - DoubleMath.log2(score2) else Double.MAX_VALUE
        } else {
            0.0
        }
    }

    override fun controlScore(chromosomeRange: ChromosomeRange): Double {
        check(isControlAvailable()) { "Control is not available" }
        return binnedCoverageQueries2!!.sumOf { it.score(chromosomeRange) } /
                binnedCoverageQueries2!!.size
    }

    override fun cleanCaches() {
        binnedCoverageQueries1?.forEach {
            it.cleanCaches()
        }
        binnedCoverageQueries2?.forEach {
            it.cleanCaches()
        }
    }

    override fun toString(): String {
        return """
            OmnipeakCompareFitInformation
                build=$build,
                data1=$data1,
                data2=$data2,
                fragment=$fragment,
                unique=$unique,
                binSize=$binSize,
                chromosomesSizes=$chromosomesSizes
        """.trimIndent()
    }

    override fun difference(loadedFitInfo: OmnipeakFitInformation): String? {
        if (loadedFitInfo !is OmnipeakCompareFitInformation) {
            return "Incompatible fit information type: ${loadedFitInfo::class.java.simpleName} " +
                    "instead of ${this::class.java.simpleName}"
        }
        if (binSize != loadedFitInfo.binSize) {
            return "Incompatible bin size: $binSize vs ${loadedFitInfo.binSize}"
        }
        if (fragment != loadedFitInfo.fragment) {
            return "Incompatible fragment: $fragment vs ${loadedFitInfo.fragment}"
        }
        if (genomeQuery().id != loadedFitInfo.genomeQuery().id) {
            return "Incompatible genome: ${genomeQuery().id} vs ${loadedFitInfo.genomeQuery().id}"
        }
        if (data1 != loadedFitInfo.data1) {
            return "Incompatible data1: $data1 vs ${loadedFitInfo.data1}"
        }
        if (data2 != loadedFitInfo.data2) {
            return "Incompatible data2: $data2 vs ${loadedFitInfo.data2}"
        }
        return null
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 5

        fun effective(
            genomeQuery: GenomeQuery,
            paths1: List<OmnipeakDataPaths>,
            paths2: List<OmnipeakDataPaths>,
            labels1: List<String>,
            labels2: List<String>,
            explicitFormat: InputFormat?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int,
            regressControl: Boolean
        ): OmnipeakCompareFitInformation {
            return OmnipeakCompareFitInformation(
                genomeQuery.build, paths1, paths2, labels1, labels2,
                explicitFormat, fragment, unique, binSize, regressControl,
                OmnipeakFitInformation.chromSizes(
                    OmnipeakModelFitExperiment.filterGenomeQueryWithData(
                        genomeQuery, paths1 + paths2, explicitFormat, fragment, unique, regressControl
                    )
                )
            )
        }
    }
}