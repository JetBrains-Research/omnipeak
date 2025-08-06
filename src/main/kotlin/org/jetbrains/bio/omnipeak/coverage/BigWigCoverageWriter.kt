package org.jetbrains.bio.omnipeak.coverage

import org.jetbrains.bio.CompressionType
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.big.FixedStepSection
import org.jetbrains.bio.big.WigSection
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.util.await
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.Callable

object BigWigCoverageWriter {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun write(
        results: OmnipeakFitResults,
        genomeQuery: GenomeQuery,
        bigWigPath: Path,
        blackListPath: Path? = null
    ) {
        val fitInfo = results.fitInfo
        check(fitInfo is OmnipeakAnalyzeFitInformation) {
            "Cannot create bigwig coverage is only available for analyze command"
        }
        fitInfo.prepareData()
        check(fitInfo.binnedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        check(fitInfo.binnedCoverageQueries!!.all { it.areCachesPresent() }) {
            "Coverage information is not available"
        }

        val blackList = if (blackListPath != null) {
            LOG.info("Loading blacklist regions: $blackListPath")
            LocationsMergingList.load(genomeQuery, blackListPath)
        } else null

        val chromosomes = genomeQuery.get()
        val scoresSumProgress = Progress {
            title = "Estimating total scores"
        }.bounded(chromosomes.size.toLong())
        var scoresSum = 0L
        val data: Map<Chromosome, IntArray> = chromosomes.associateWith { chromosome ->
            val data = fitInfo.dataQuery.apply(chromosome)
            val scores = data.sliceAsInt(data.labels.first())
            scoresSum += scores.sum().toLong()
            scoresSumProgress.report()
            scores
        }

        val cpmScale = 1e6 / scoresSum
        LOG.debug("Total scores: ${scoresSum}, scale: $cpmScale")

        val bin = fitInfo.binSize
        if (blackList != null) {
            chromosomes.forEach { chromosome ->
                val coverage: IntArray = data[chromosome] as IntArray
                blackList[chromosome, Strand.PLUS].forEach { (startOffset, endOffset) ->
                    for (b in startOffset / bin until endOffset / bin) {
                        coverage[b] = 0
                    }
                }
                blackList[chromosome, Strand.MINUS].forEach { (startOffset, endOffset) ->
                    for (b in startOffset / bin until endOffset / bin) {
                        coverage[b] = 0
                    }
                }
            }
        }

        write(
            { (start, end, chr) ->
                val coverage: IntArray = data[chr] as IntArray
                (start / bin until end / bin).sumOf { coverage[it] } * cpmScale
            },
            genomeQuery, fitInfo.binSize, bigWigPath
        )
    }

    fun write(
        scoreCompute: (ChromosomeRange) -> Double,
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
                            this.add(scoreCompute(range.on(chromosome)).toFloat())
                        }
                    }
                    wigSections[i] = section
                    collectDataProgress.report()
                }
            }.await(true)
        } finally {
            collectDataProgress.done()
        }

        BigWigFile.write(
            wigSections.asIterable(),
            chromosomes.map { it.name to it.length },
            bigWigPath,
            // Old compression for compatibility with IGV 2.3.92 browser
            compression = CompressionType.DEFLATE
        ) {}
        LOG.info("Successfully created bigwig file: $bigWigPath")
    }
}