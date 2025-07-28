package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.CompressionType
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.big.FixedStepSection
import org.jetbrains.bio.big.WigSection
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.util.await
import org.jetbrains.bio.util.time
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.slf4j.event.Level
import java.nio.file.Path
import java.util.concurrent.Callable

object BigWigCoverageWriter {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun write(
        coverageCompute: (ChromosomeRange) -> Int,
        genomeQuery: GenomeQuery,
        binSize: Int,
        bigWigPath: Path
    ) {
        val chromosomes = genomeQuery.get()

        val collectDataProgress = Progress {
            title = "Creating bigwig file"
        }.bounded(chromosomes.size.toLong())

        val fakeSection = FixedStepSection("chrFake", binSize, span = binSize)
        val wigSections = Array<WigSection>(chromosomes.size) {
            fakeSection
        }

        try {
            chromosomes.mapIndexed { i, chromosome ->
                Callable {
                    val section = FixedStepSection(
                        chromosome.name, 0,
                        step = binSize, span = binSize
                    ).apply {
                        chromosome.range.slice(binSize).forEach { range ->
                            this.add(coverageCompute(range.on(chromosome)).toFloat())
                        }
                    }
                    wigSections[i] = section
                    collectDataProgress.report()
                }
            }.await(true)
        } finally {
            collectDataProgress.done()
        }

        LOG.time(Level.INFO, "Saving bigwig scores: $bigWigPath") {
            BigWigFile.write(
                wigSections.asIterable(),
                chromosomes.map { it.name to it.length },
                bigWigPath,
                // Old compression for compatibility with IGV 2.3.92 browser
                compression = CompressionType.DEFLATE
            ) {}
        }

    }
}