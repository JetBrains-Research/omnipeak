package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import java.util.*

/**
 * The common interface for all fit information classes.
 *
 * [OmnipeakFitInformation] instance is designed to contain all information necessary to uniquely identify the input
 * of an Omnipeak-like model fitting experiment.
 *
 * [OmnipeakFitInformation] object is a part of [OmnipeakFitResults].
 */
interface OmnipeakFitInformation {

    /**
     * A unique string identifier.
     */
    val id: String

    /** Genome build. */
    val build: String

    /** Fragment size */
    val fragment: Fragment

    /** Keep unique reads only. */
    val unique: Boolean

    /** Bin size in base pairs. */
    val binSize: Int

    /** Whether to regress out control signal from the data. */
    val regressControl: Boolean

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
     */
    fun score(chromosomeRange: ChromosomeRange): Double

    fun isControlAvailable(): Boolean

    /**
     * Returns scaled control score or null if not available.
     * Call [prepareData] beforehand!
     */
    fun controlScore(chromosomeRange: ChromosomeRange): Double

    fun cleanCaches()

    fun difference(loadedFitInfo: OmnipeakFitInformation): String?

    fun genomeQuery(): GenomeQuery =
        GenomeQuery(Genome[build, chromosomesSizes], *chromosomesSizes.keys.toTypedArray())

    fun containsChromosomeInfo(chromosome: Chromosome): Boolean =
        chromosome.name in chromosomesSizes

    companion object
}
