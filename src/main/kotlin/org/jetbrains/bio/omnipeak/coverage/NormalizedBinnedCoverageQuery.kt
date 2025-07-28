package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.query.CachingQuery

class NormalizedBinnedCoverageQuery(val ncq: NormalizedCoverageQuery, val binSize: Int) :
    BinnedCoverageQuery, CachingQuery<Chromosome, IntArray>() {

    override val id: String
        get() = "${ncq.id}_binned"

    override val description: String
        get() = super<CachingQuery>.description + " binned to $binSize bp"

    override fun getUncached(input: Chromosome): IntArray {
        return input.range.slice(binSize).mapToInt { range ->
            ncq.apply(range.on(input))
        }.toArray()
    }

    override fun score(chromosomeRange: ChromosomeRange): Double =
        ncq.score(chromosomeRange)

    override fun controlAvailable(): Boolean =
        ncq.isControlAvailable()

    override fun controlScore(chromosomeRange: ChromosomeRange): Double =
        ncq.controlScore(chromosomeRange)

    override fun controlNormalizedScore(chromosomeRange: ChromosomeRange): Int =
        ncq.controlNormalizedScore(chromosomeRange)

    override fun areCachesPresent(): Boolean =
        ncq.areCachesPresent()

    override fun cleanCaches() =
        ncq.cleanCaches()

}