package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.omnipeak.InputFormat
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_BIN
import java.nio.file.Path

interface BinnedCoverageQuery : Query<Chromosome, IntArray> {

    fun score(chromosomeRange: ChromosomeRange): Double

    fun controlAvailable(): Boolean

    fun controlScore(chromosomeRange: ChromosomeRange): Double

    fun controlNormalizedScore(chromosomeRange: ChromosomeRange): Int

    fun areCachesPresent(): Boolean

    fun cleanCaches()


    companion object {

        fun create(
            genomeQuery: GenomeQuery,
            treatmentPath: Path,
            controlPath: Path?,
            explicitFormat: InputFormat?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int,
            regressControl: Boolean,
            forceCaching: Boolean = false,
            showLibraryInfo: Boolean = true,
        ): BinnedCoverageQuery {
            val format = explicitFormat ?: InputFormat.guess(treatmentPath)
            if (format == InputFormat.BIGWIG) {
                return BigWigBinnedCoverageQuery(
                    genomeQuery,
                    treatmentPath,
                    controlPath,
                    binSize,
                    regressControl,
                    showLibraryInfo
                )
            }
            val ncq = NormalizedCoverageQuery(
                genomeQuery,
                treatmentPath,
                controlPath,
                explicitFormat?.internalReadsFormat,
                fragment,
                unique,
                binSize,
                regressControl,
                showLibraryInfo
            )
            if (forceCaching) {
                ncq.treatmentReads.get()
                ncq.controlReads?.get()
            }
            return NormalizedBinnedCoverageQuery(ncq, binSize)
        }

    }
}
