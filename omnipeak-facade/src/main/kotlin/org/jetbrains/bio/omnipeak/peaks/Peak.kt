package org.jetbrains.bio.omnipeak.peaks

import org.jetbrains.bio.genome.*

/**
 * @param value Foldchange or coverage
 * @param mlogpvalue -log10(pvalue)
 * @param mlogqvalue -log10(qvalue)
 * @param score A score between 0 and 1000
 */
data class Peak(
    val chromosome: Chromosome,
    val startOffset: Int,
    val endOffset: Int,
    val mlogpvalue: Double,
    val mlogqvalue: Double,
    var value: Double = 0.0,
    val score: Int
) : Comparable<Peak>, LocationAware {
    val range: Range
        get() = Range(startOffset, endOffset)

    override val location: Location
        get() = Location(startOffset, endOffset, chromosome, Strand.PLUS)

    override fun compareTo(other: Peak): Int =
        compareValuesBy(this, other, { it.chromosome.name }, { it.startOffset }, { it.endOffset })

    companion object
}
