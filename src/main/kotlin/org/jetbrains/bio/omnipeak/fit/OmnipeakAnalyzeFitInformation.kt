package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.omnipeak.InputFormat
import org.jetbrains.bio.omnipeak.coverage.BigWigBinnedCoverageQuery
import org.jetbrains.bio.omnipeak.coverage.BinnedCoverageQuery
import org.jetbrains.bio.omnipeak.coverage.NormalizedBinnedCoverageQuery
import org.jetbrains.bio.util.isAccessible
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz

/**
 * Since all the chromosomes are squashed in [OmnipeakModelFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 *
 * [labels] refer to the coverage dataframe column labels, not to the supervised learning annotations.
 */
data class OmnipeakAnalyzeFitInformation(
    override val build: String,
    val paths: List<OmnipeakDataPaths>,
    val labels: List<String>,
    override val explicitFormat: InputFormat?,
    override val fragment: Fragment,
    override val unique: Boolean,
    override val binSize: Int,
    override val regressControl: Boolean,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : OmnipeakFitInformation {

    constructor(
        genomeQuery: GenomeQuery,
        paths: List<OmnipeakDataPaths>,
        labels: List<String>,
        explicitFormat: InputFormat?,
        fragment: Fragment,
        unique: Boolean,
        binSize: Int,
        regressControl: Boolean
    ) : this(
        genomeQuery.build, paths, labels, explicitFormat,
        fragment, unique, binSize, regressControl,
        OmnipeakFitInformation.chromSizes(genomeQuery)
    )

    override val id: String
        get() = generateId(paths, fragment, binSize, unique, regressControl)


    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            check(scoresAvailable()) { "Scores are not available!" }
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return binnedCoverageQueries!!.binnedCoverageDataFrame(
                        input, labels.toTypedArray()
                    )
                }

                override val id: String
                    get() = reduceIds(
                        binnedCoverageQueries!!.zip(labels)
                            .flatMap { (s, l) -> listOf(s.id, l) }
                    )
            }
        }

    @Transient
    var binnedCoverageQueries: List<BinnedCoverageQuery>? = null

    override fun prepareData() {
        if (binnedCoverageQueries == null) {
            val pathsPresent = paths.all { it.treatment.isAccessible() && it.control?.isAccessible() != false }
            if (pathsPresent) {
                binnedCoverageQueries = paths.map {
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
                        showLibraryInfo = true
                    )
                }
            }
        }
    }

    override fun scoresAvailable(): Boolean {
        prepareData()
        return binnedCoverageQueries != null
    }

    /**
     * Returns average coverage by tracks
     */
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(scoresAvailable()) { "Scores are not available!" }
        return binnedCoverageQueries!!.sumOf { it.score(chromosomeRange) } /
                binnedCoverageQueries!!.size
    }

    override fun isControlAvailable(): Boolean =
        binnedCoverageQueries?.all { it.controlAvailable() } ?: false

    override fun controlScore(chromosomeRange: ChromosomeRange): Double {
        check(isControlAvailable()) { "Control is not available" }
        return binnedCoverageQueries!!.sumOf { it.controlScore(chromosomeRange) } /
                binnedCoverageQueries!!.size
    }

    override fun cleanCaches() {
        binnedCoverageQueries?.forEach {
            it.cleanCaches()
        }
    }

    override fun toString(): String {
        return """
            OmnipeakAnalyzeFitInformation
                build='$build',
                paths=$paths,
                fragment=$fragment,
                unique=$unique,
                binSize=$binSize,
                chromosomesSizes=$chromosomesSizes
        """.trimIndent()
    }

    override fun difference(loadedFitInfo: OmnipeakFitInformation): String? {
        if (loadedFitInfo !is OmnipeakAnalyzeFitInformation) {
            return "Incompatible fit information type: ${loadedFitInfo::class.java.simpleName} " +
                    "instead of ${this::class.java.simpleName}"
        }
        if (binSize != loadedFitInfo.binSize) {
            return "Incompatible bin size: $binSize vs ${loadedFitInfo.binSize}"
        }
        if (fragment != loadedFitInfo.fragment) {
            return "Incompatible fragment: $fragment vs ${loadedFitInfo.fragment}"
        }
        if (unique != loadedFitInfo.unique) {
            return "Incompatible unique: $unique vs ${loadedFitInfo.unique}"
        }
        if (chromosomesSizes != loadedFitInfo.chromosomesSizes) {
            return "Incompatible chromosomes sizes: $chromosomesSizes vs ${loadedFitInfo.chromosomesSizes}"
        }
        if (paths != loadedFitInfo.paths) {
            return "Incompatible paths: $paths vs ${loadedFitInfo.paths}"
        }
        if (genomeQuery().id != loadedFitInfo.genomeQuery().id) {
            return "Incompatible genome: ${genomeQuery().id} vs ${loadedFitInfo.genomeQuery().id}"
        }
        return null
    }

    fun getTreatmentComputable(): ((ChromosomeRange) -> Double)? {
        if (!scoresAvailable()) {
            return null
        }
        when (binnedCoverageQueries!!.first()) {
            is NormalizedBinnedCoverageQuery -> {
                val bncq = this.binnedCoverageQueries!! as List<NormalizedBinnedCoverageQuery>
                if (!bncq.all { it.areCachesPresent() }) {
                    return null
                }

                // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
                // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
                val treatmentCoverages = bncq.map { it.ncq.treatmentReads.get() }
                return { chromosomeRange: ChromosomeRange ->
                    treatmentCoverages.sumOf { it.getBothStrandsCoverage(chromosomeRange) }.toDouble() /
                            treatmentCoverages.size
                }
            }

            is BigWigBinnedCoverageQuery -> {
                val bncq = this.binnedCoverageQueries!! as List<BigWigBinnedCoverageQuery>
                return { chromosomeRange: ChromosomeRange ->
                    bncq.sumOf { it.score(chromosomeRange) } / bncq.size
                }
            }

            else -> {
                throw IllegalStateException(
                    "Unexpected binned coverage query: " +
                            "${binnedCoverageQueries!!.first()::class.java.simpleName}"
                )
            }
        }
    }

    fun getTreatmentControlComputable(): ((ChromosomeRange) -> Pair<Double, Double>)? {
        if (!scoresAvailable()) {
            return null
        }
        when (val query = this.binnedCoverageQueries!!.first()) {
            is NormalizedBinnedCoverageQuery -> {
                val bncq = this.binnedCoverageQueries!! as List<NormalizedBinnedCoverageQuery>
                if (!bncq.all { it.areCachesPresent() }) {
                    return null
                }

                // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
                // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
                val treatmentCoverages = bncq.map { it.ncq.treatmentReads.get() }
                val controlCoverages = bncq.map { it.ncq.controlReads?.get() }
                val controlScales = bncq.map { it.ncq.coveragesNormalizedInfo.controlScale }
                return { chromosomeRange: ChromosomeRange ->
                    val score =
                        treatmentCoverages.sumOf {
                            it.getBothStrandsCoverage(chromosomeRange)
                        }.toDouble() / treatmentCoverages.size
                    val controlScore =
                        controlCoverages.zip(controlScales).sumOf { (c, s) ->
                            if (c != null) c.getBothStrandsCoverage(chromosomeRange) * s else 0.0
                        } / controlCoverages.size
                    score to controlScore
                }
            }

            is BigWigBinnedCoverageQuery -> {
                val bncq = this.binnedCoverageQueries!! as List<BigWigBinnedCoverageQuery>
                return { chromosomeRange: ChromosomeRange ->
                    val score = bncq.sumOf { it.score(chromosomeRange) } / bncq.size
                    val controlScore = bncq.sumOf {
                        if (it.controlAvailable()) it.controlScore(chromosomeRange) else 0.0
                    } / bncq.size
                    score to controlScore
                }
            }

            else -> {
                throw IllegalStateException("Unexpected binned coverage query: ${query::class.java.simpleName}")
            }
        }
    }


    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 5

        /**
         * Generates a model ID based on the provided parameters.
         *
         * @param paths A list of [OmnipeakDataPaths] representing treatment and control pairs.
         * @param binSize The bin size.
         * @param fragment The fragment type.
         * @param unique Indicates whether the experiment is unique.
         * @return The generated ID
         */
        fun generateId(
            paths: List<OmnipeakDataPaths>,
            fragment: Fragment,
            binSize: Int,
            unique: Boolean,
            regressControl: Boolean
        ) = reduceIds(
            paths.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOfNotNull(fragment.nullableInt, binSize).map { it.toString() } +
                    listOfNotNull(if (unique) "unique" else null) +
                    listOfNotNull(if (!regressControl) "no-regress-control" else null)
        )

        fun createFitInformation(
            genomeQuery: GenomeQuery,
            paths: List<OmnipeakDataPaths>,
            labels: List<String>,
            explicitFormat: InputFormat?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int,
            regressControl: Boolean
        ): OmnipeakAnalyzeFitInformation {
            val genomeQueryWithData =
                OmnipeakModelFitExperiment.filterGenomeQueryWithData(
                    genomeQuery,
                    paths,
                    explicitFormat,
                    fragment,
                    unique,
                    regressControl
                )
            return OmnipeakAnalyzeFitInformation(
                genomeQueryWithData.build,
                paths,
                labels,
                explicitFormat,
                fragment,
                unique,
                binSize,
                regressControl,
                OmnipeakFitInformation.chromSizes(genomeQueryWithData)
            )
        }
    }
}
