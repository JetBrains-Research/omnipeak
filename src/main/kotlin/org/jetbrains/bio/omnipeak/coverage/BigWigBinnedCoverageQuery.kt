package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.big.BigSummary
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_BIN
import org.jetbrains.bio.util.name
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import org.slf4j.LoggerFactory
import java.nio.file.Path
import kotlin.math.max

/**
 * A BinnedQuery implementation that reads coverage data from BigWig files.
 */
class BigWigBinnedCoverageQuery : BinnedCoverageQuery {

    private val genomeQuery: GenomeQuery
    private val treatmentPath: Path
    private val controlPath: Path?
    private val binSize: Int

    constructor(
        genomeQuery: GenomeQuery,
        treatmentPath: Path,
        controlPath: Path?,
        binSize: Int = OMNIPEAK_DEFAULT_BIN,
        showLibraryInfo: Boolean = true
    ) {
        this.genomeQuery = genomeQuery
        this.treatmentPath = treatmentPath
        this.controlPath = controlPath
        this.binSize = binSize
        if (showLibraryInfo) {
            showLibraryInfo()
        }
    }

    // Lazy initialization of BigWig files
    private val treatmentBigWig by lazy { BigWigFile.read(treatmentPath) }
    private val treatmentTotalCoverage by lazy { treatmentBigWig.totalSummary.sum }

    private val controlBigWig by lazy { controlPath?.let { BigWigFile.read(it) } }
    private val controlTotalCoverage by lazy { controlBigWig?.totalSummary?.sum ?: 0.0 }

    override fun score(chromosomeRange: ChromosomeRange): Double {
        val chromosome = chromosomeRange.chromosome
        val chrName = chromosome.name

        // Check if chromosome exists in BigWig file
        val matchedBfChr = findMatchedChromosome(treatmentBigWig, chrName)
            ?: return 0.0 // Return 0 if chromosome not found

        // Use summarize to get coverage data for the range
        val summaries = treatmentBigWig.summarize(
            matchedBfChr,
            chromosomeRange.startOffset,
            chromosomeRange.endOffset,
            1 // Just one bin for the entire range
        )

        return if (summaries.isNotEmpty()) {
            summaries[0].sum
        } else {
            0.0
        }
    }

    override fun controlAvailable(): Boolean {
        return controlBigWig != null
    }

    override fun controlScore(chromosomeRange: ChromosomeRange): Double {
        require(controlBigWig != null) { "Control BigWig file is not available" }

        val chromosome = chromosomeRange.chromosome
        val chrName = chromosome.name

        // Check if chromosome exists in BigWig file
        val matchedBfChr = findMatchedChromosome(controlBigWig!!, chrName)
            ?: return 0.0 // Return 0 if chromosome not found

        // Use summarize to get coverage data for the range
        val summaries = controlBigWig!!.summarize(
            matchedBfChr,
            chromosomeRange.startOffset,
            chromosomeRange.endOffset,
            1 // Just one bin for the entire range
        )

        return if (summaries.isNotEmpty()) {
            summaries[0].sum
        } else {
            0.0
        }
    }

    override fun controlNormalizedScore(chromosomeRange: ChromosomeRange): Int {
        if (!controlAvailable()) {
            return score(chromosomeRange).toInt()
        }
        // Scale control to treatment
        val controlScale = treatmentTotalCoverage / controlTotalCoverage
        return max(0.0, score(chromosomeRange) - controlScore(chromosomeRange) * controlScale).toInt()
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

        // Check if chromosome exists in BigWig file
        val matchedBfChr = findMatchedChromosome(treatmentBigWig, chrName)
            ?: return IntArray(chrLength / binSize + 1) // Return empty array if chromosome not found

        // Calculate number of bins
        val binsCount = chrLength / binSize + 1

        // Use summarize to get coverage data for the entire chromosome
        val treatmentCoverage = treatmentBigWig.summarize(
            matchedBfChr,
            0,
            chrLength,
            binsCount
        )
        checkNonNegative(treatmentCoverage, "Treatment")

        // Convert summaries to IntArray
        if (!controlAvailable()) {
            return IntArray(binsCount) { i ->
                if (i < treatmentCoverage.size) {
                    treatmentCoverage[i].sum.toInt()
                } else {
                    0
                }
            }
        }
        // Scale control to treatment
        val controlScale = treatmentTotalCoverage / controlTotalCoverage
        val controlCoverage = controlBigWig!!.summarize(
            matchedBfChr,
            0,
            chrLength,
            binsCount
        )
        checkNonNegative(controlCoverage, "Control")

        return IntArray(binsCount) { i ->
            if (i < treatmentCoverage.size) {
                max(0.0, treatmentCoverage[i].value() - controlCoverage[i].value() * controlScale).toInt()
            } else {
                0
            }
        }
    }

    private fun checkNonNegative(coverage: List<BigSummary>, title: String) {
        coverage.forEach {
            check(it.sum >= 0) {"$title: negative values detected"}
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
    }
}

fun BigSummary.value() = this.sum / max(1, this.count)
