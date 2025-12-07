package org.jetbrains.bio.omnipeak.coverage

import org.apache.commons.math3.stat.descriptive.rank.Percentile
import org.jetbrains.bio.big.BigSummary
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.name
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import org.slf4j.LoggerFactory
import java.nio.file.Path
import kotlin.math.max
import kotlin.math.min

/**
 * A BinnedQuery implementation that reads coverage data from BigWig files.
 */
class BigWigBinnedCoverageQuery(
    private val genomeQuery: GenomeQuery,
    private val treatmentPath: Path,
    private val controlPath: Path?,
    private val binSize: Int,
    private val regressControl: Boolean
) : BinnedCoverageQuery {

    constructor(
        genomeQuery: GenomeQuery,
        treatmentPath: Path,
        controlPath: Path?,
        binSize: Int,
        regressControl: Boolean,
        showLibraryInfo: Boolean = true
    ) : this(genomeQuery, treatmentPath, controlPath, binSize, regressControl) {
        if (showLibraryInfo) {
            showLibraryInfo()
        }
    }

    // Lazy initialization of BigWig files
    private val treatmentBigWig by lazy { BigWigFile.read(treatmentPath) }

    private val treatmentTotalCoverage by lazy { treatmentBigWig.totalSummary.sum }

    private val treatmentTopPercentile by lazy {
        val matches = genomeQuery.get().mapNotNull {
            // Check if chromosome exists in BigWig file - 0 if chromosome not found
            val match = findMatchedChromosome(treatmentBigWig, it.name)
            if (match != null) {
                it to match
            } else null
        }
        if (matches.isEmpty()) return@lazy 0.0
        val (chromosome, matchedBfChr) = matches.maxBy { it.first.length }
        // Calculate number of bins
        val binsCount = chromosome.length / binSize + 1
        // Use summarize to get coverage data for the entire chromosome
        val data = treatmentBigWig.summarize(
            matchedBfChr, 0, chromosome.length, binsCount
        ).map { it.value() }.filter { it > 0 }.toDoubleArray()
        // Do not use StatUtils.percentile(scores.toArray(), XX) to avoid redundant
        //   score array copying
        val percentile = object : Percentile(TRIM_PERCENTILE_MAX) {
            // force Percentile not to copy scores
            override fun getWorkArray(values: DoubleArray?, begin: Int, length: Int) = data
        }.evaluate(data)
        LOG.debug("Treatment top $TRIM_PERCENTILE_MAX percentile: ${"%.3f".format(percentile)}")
        return@lazy percentile
    }

    // We don't know the scale of bigWig, we want to keep it in ranges
    private val treatmentScale by lazy {
        val s = if (treatmentTopPercentile < TOP_SIGNAL_MIN_AVG * binSize)
            TOP_SIGNAL_MIN_AVG * binSize / treatmentTopPercentile
        else if (treatmentTopPercentile > TOP_SIGNAL_MAX_AVG * binSize)
            TOP_SIGNAL_MAX_AVG * binSize / treatmentTopPercentile
        else 1.0
        LOG.debug("Treatment scale: ${"%.3f".format(s)}")
        return@lazy s
    }

    private val controlScale by lazy { treatmentTotalCoverage * treatmentScale / controlTotalCoverage }

    private val controlBigWig by lazy { controlPath?.let { BigWigFile.read(it) } }

    private val controlTotalCoverage by lazy { controlBigWig?.totalSummary?.sum ?: 0.0 }

    override fun score(chromosomeRange: ChromosomeRange): Double {
        val chromosome = chromosomeRange.chromosome
        val chrName = chromosome.name

        // Check if chromosome exists in BigWig file - return 0 if chromosome not found
        val matchedBfChr =
            findMatchedChromosome(treatmentBigWig, chrName) ?: return 0.0

        // Use summarize to get coverage data for the range - one bin for the entire range
        val summaries = treatmentBigWig.summarize(
            matchedBfChr, chromosomeRange.startOffset, chromosomeRange.endOffset, 1
        )

        return if (summaries.isNotEmpty()) treatmentScore(summaries[0].sum) else 0.0
    }

    override fun controlAvailable(): Boolean {
        return controlBigWig != null
    }

    override fun controlScore(chromosomeRange: ChromosomeRange): Double {
        require(controlBigWig != null) { "Control BigWig file is not available" }

        val chromosome = chromosomeRange.chromosome
        val chrName = chromosome.name

        // Check if chromosome exists in BigWig file - 0 if chromosome not found
        val matchedBfChr =
            findMatchedChromosome(controlBigWig!!, chrName) ?: return 0.0

        // Use summarize to get coverage data for the range - one bin for the entire range
        val summaries = controlBigWig!!.summarize(
            matchedBfChr, chromosomeRange.startOffset, chromosomeRange.endOffset, 1
        )

        return if (summaries.isNotEmpty()) controlScore(summaries[0].sum) else 0.0
    }

    override fun controlNormalizedScore(chromosomeRange: ChromosomeRange): Int {
        if (!regressControl || !controlAvailable()) {
            return score(chromosomeRange).toInt()
        }
        // Scale control to treatment
        return max(0.0, score(chromosomeRange) - controlScore(chromosomeRange)).toInt()
    }

    override fun areCachesPresent(): Boolean {
        // BigWig files don't use caches in the same way as other query types
        return true
    }

    override fun cleanCaches() {
        // No caches to clean for BigWig files
    }

    override val id: String
        get() = reduceIds(
            listOfNotNull(
                treatmentPath.stemGz, controlPath?.stemGz, binSize.toString()
            )
        )

    override fun apply(t: Chromosome): IntArray {
        val chrName = t.name
        val chrLength = t.length

        // Check if chromosome exists in BigWig file - return empty array if chromosome not found
        val matchedBfChr = findMatchedChromosome(treatmentBigWig, chrName)
            ?: return IntArray(chrLength / binSize + 1)

        // Calculate number of bins
        val binsCount = chrLength / binSize + 1

        // Use summarize to get coverage data for the entire chromosome
        val treatmentCoverage = treatmentBigWig.summarize(
            matchedBfChr, 0, chrLength, binsCount
        )
        checkNonNegative(treatmentCoverage, "Treatment")

        // Convert summaries to IntArray
        if (!regressControl || !controlAvailable()) {
            return IntArray(binsCount) { i ->
                if (i < treatmentCoverage.size) {
                    treatmentScore(treatmentCoverage[i].value()).toInt()
                } else {
                    0
                }
            }
        }
        // Scale control to treatment
        val controlCoverage = controlBigWig!!.summarize(
            matchedBfChr, 0, chrLength, binsCount
        )
        checkNonNegative(controlCoverage, "Control")

        return IntArray(binsCount) { i ->
            if (i < treatmentCoverage.size) {
                max(0.0,
                    treatmentScore(treatmentCoverage[i].value()) - controlScore(controlCoverage[i].value())
                ).toInt()
            } else {
                0
            }
        }
    }


    private fun treatmentScore(value: Double): Double = min(value, treatmentTopPercentile) * treatmentScale

    private fun controlScore(value: Double): Double = value * controlScale

    private fun checkNonNegative(coverage: List<BigSummary>, title: String) {
        coverage.forEach {
            check(it.sum >= 0) { "$title: negative values detected" }
        }
    }

    /**
     * Helper function to find a matching chromosome name in the BigWig file.
     * Handles cases where chromosome names might differ (e.g., "chr1" vs "1").
     */
    private fun findMatchedChromosome(bigWigFile: BigWigFile, chrName: String): String? {
        return if (bigWigFile.chromosomes.containsValue(chrName)) {
            chrName
        } else {
            (genomeQuery.genome.chromosomeNamesToAltNamesMap[chrName] ?: emptyList()).firstOrNull {
                bigWigFile.chromosomes.containsValue(it)
            }
        }
    }

    private fun showLibraryInfo() {
        LOG.info("Library: ${treatmentPath.name}, Depth: ${"%,d".format(treatmentTotalCoverage.toLong())}")
        if (controlPath != null) {
            LOG.info("Library: ${controlPath.name}, Depth: ${"%,d".format(controlTotalCoverage.toLong())}")
            LOG.info("Control normalization: ${"%.3f".format(treatmentTotalCoverage / controlTotalCoverage)}")
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(BigWigBinnedCoverageQuery::class.java)

        // Equal to AbstractBedTrackView.TRIM_PERCENTILE_MAX = 99.0
        const val TRIM_PERCENTILE_MAX = 99.0

        const val TOP_SIGNAL_MAX_AVG = 2.0
        const val TOP_SIGNAL_MIN_AVG = 0.2
    }
}

fun BigSummary.value() = this.sum / max(1, this.count)
