package org.jetbrains.bio.omnipeak.peaks

import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_FRAGMENTATION_MAX_GAP_BP
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_SCORE_BLOCKS
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.omnipeak.fit.OmnipeakModelFitExperiment
import org.jetbrains.bio.omnipeak.modes.HartiganDip
import org.jetbrains.bio.omnipeak.modes.KdeModeFinder
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.estimateCandidatesNumberLens
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.estimateGap
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.estimateSensitivity
import org.jetbrains.bio.omnipeak.peaks.Signal.clipPeakBySignal
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.await
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.ln
import kotlin.math.max
import kotlin.math.min


object OmnipeakModelToPeaks {

    val LOG: Logger = LoggerFactory.getLogger(OmnipeakModelToPeaks.javaClass)

    /**
     * Main method to compute peaks from the model.
     *
     * 1) Estimate posterior error probabilities PEP
     * 2) Find sensitivity setting as stable point wrt candidates number and average length.
     * 3) Compute required gap to compensate for extra high fragmentation, i.e. fraction of peaks when increasing gap
     * 4) Assign p-value to each peak based on combined p-values for cores (consequent foreground bins).
     *    In case when control track is present, we use Poisson CDF to estimate log P-value;
     *    otherwise, an average log PEP (posterior error probability) for bins in blocks is used.
     *    N% top significant blocks scores are aggregated using length-weighted average as P for peak.
     * 5) Compute qvalues by peaks p-values, filter by alpha.
     * 6) Fine-tuning boundaries of point-wise peaks according to the local signal.
     */
    fun getPeaks(
        omnipeakFitResults: OmnipeakFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        multipleTesting: MultipleTesting,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        summits: Boolean,
        fragmentationThreshold: Int,
        clip: Double,
        blackListPath: Path? = null,
        parallel: Boolean = true,
        name: String? = null,
        cancellableState: CancellableState? = null,
    ): OmnipeakResult {
        val fitInfo = omnipeakFitResults.fitInfo
        // Prepare fit information for scores computations
        fitInfo.prepareData()

        // Collect candidates from model
        val logNullMembershipsMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                // Cannot be null because of GenomeMap and empty F64Array is not supported
                return@genomeMap FAKE_ARRAY
            }
            getLogNulls(omnipeakFitResults, chromosome)
        }
        val bitList2reuseMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            BitList(logNullMembershipsMap[chromosome].length)
        }

        // Estimate sensitivity
        val sensitivity = if (sensitivityCmdArg != null) {
            sensitivityCmdArg
        } else {
            estimateSensitivity(
                genomeQuery, omnipeakFitResults, logNullMembershipsMap, bitList2reuseMap,
                parallel, name, cancellableState
            )
        }
        LOG.info("${name ?: ""} Selecting candidates with sensitivity: $sensitivity")

        // Estimate gap
        val gap = when {
            gapCmdArg != null -> gapCmdArg
            summits -> 0
            else -> {
                LOG.info("${name ?: ""} Analysing fragmentation...")
                val candidateGapNs = IntArray(OMNIPEAK_FRAGMENTATION_MAX_GAP_BP / fitInfo.binSize) {
                    estimateCandidatesNumberLens(
                        genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                        sensitivity, it
                    ).n
                }
                estimateGap(candidateGapNs, name, fitInfo.binSize, fragmentationThreshold)
            }
        }
        LOG.info("${name ?: ""} Candidates selection with gap: $gap")

        val candidatesMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(
                chromosome, logNullMemberships, bitList2reuse,
                sensitivity, gap
            )
        }

        val singleModeLength =
            estimateSingleModeLength(genomeQuery, fitInfo, candidatesMap, parallel)

        val peaks = candidatesToPeaks(
            genomeQuery,
            fitInfo,
            logNullMembershipsMap,
            candidatesMap,
            summits,
            singleModeLength,
            clip,
            fdr,
            multipleTesting,
            blackListPath,
            name,
            parallel,
            cancellableState
        )
        return OmnipeakResult(fdr, sensitivity, gap, peaks)
    }

    fun candidatesToPeaks(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        candidatesMap: GenomeMap<List<Range>>,
        summits: Boolean,
        singleModeLength: Int,
        clip: Double,
        fdr: Double,
        testing: MultipleTesting,
        blackListPath: Path?,
        name: String?,
        parallel: Boolean,
        cancellableState: CancellableState?
    ): List<Peak> {
        var candidatesOrSummitsMap = candidatesMap
        if (summits) {
            LOG.info("${name ?: ""} Processing summits")
            if (singleModeLength > 0) {
                candidatesOrSummitsMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
                    cancellableState?.checkCanceled()
                    if (!fitInfo.containsChromosomeInfo(chromosome)) {
                        return@genomeMap emptyList()
                    }
                    val chrCandidates = candidatesMap[chromosome]
                    getChromosomeSummits(chromosome, fitInfo, chrCandidates, singleModeLength)
                }
            }
        }

        // Load blacklisted file
        val blackList = if (blackListPath != null) {
            LOG.info("${name ?: ""} Loading blacklisted regions: $blackListPath")
            LocationsMergingList.load(genomeQuery, blackListPath)
        } else null

        val peakScorer = PeakScorer.create(fitInfo, logNullMembershipsMap)

        val (avgSignalDensity, avgNoiseDensity) = peakScorer.analyzeSignalAndNoise(
            genomeQuery,
            fitInfo,
            candidatesOrSummitsMap,
            false
        )
        LOG.debug("${name ?: ""} Signal density $avgSignalDensity, noise density $avgNoiseDensity")

        // Compute candidate locations and blocks of significance
        val candidateCoordinatesMap = genomeMap(genomeQuery, parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList()
            }
            val chrCandidates = candidatesOrSummitsMap[chromosome]
            val offsets = fitInfo.offsets(chromosome)
            val logNullMemberships = logNullMembershipsMap[chromosome]
            processCandidateCoordinates(
                chromosome, chrCandidates, fitInfo, logNullMemberships, offsets,
                blackList, avgSignalDensity, avgNoiseDensity, clip, summits,
                cancellableState
            )
        }

        // Estimate total number of candidate locations and compute indices by chromosome
        var candidateLocationsTotal = 0
        val chromosomeLocationsOffsets = hashMapOf<Chromosome, Pair<Int, Int>>()
        genomeQuery.get().forEach { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@forEach
            }
            val n = candidateCoordinatesMap[chromosome].size
            chromosomeLocationsOffsets[chromosome] = candidateLocationsTotal to candidateLocationsTotal + n
            candidateLocationsTotal += n
        }

        // Finally create peaks list
        val peaks: List<Peak>

        // Return empty list when nothing found
        if (candidateLocationsTotal != 0) {
            // Estimate global pvalues
            val logPValsMap = genomeMap(genomeQuery, parallel) { chromosome ->
                cancellableState?.checkCanceled()
                if (!fitInfo.containsChromosomeInfo(chromosome)) {
                    return@genomeMap FAKE_ARRAY // Return fake array to keep typing
                }
                return@genomeMap estimateCandidatesLogPs(
                    chromosome, candidateCoordinatesMap[chromosome], peakScorer, cancellableState
                )
            }

            // Adjust pvalues globally
            val genomeLogPVals = F64Array(candidateLocationsTotal)
            genomeQuery.get().forEach { chromosome ->
                if (!fitInfo.containsChromosomeInfo(chromosome)) {
                    return@forEach
                }
                val (start, end) = chromosomeLocationsOffsets[chromosome]!!
                if (start == end) return@forEach
                val chrLogPVals = logPValsMap[chromosome]
                check(chrLogPVals.length == end - start)
                for (i in start until end) genomeLogPVals[i] = chrLogPVals[i - start]
            }

            // Adjust globally log pValues -> log qValues
            LOG.info("${name ?: ""} Adjusting pvalues ${testing.description}, N=${genomeLogPVals.length}")
            val genomeLogQVals = if (testing == MultipleTesting.BH)
                Fdr.qvalidate(genomeLogPVals, logResults = true)
            else
                F64Array(genomeLogPVals.length) { genomeLogPVals[it] + ln(genomeLogPVals.length.toDouble()) }

            // Collect peaks together from all chromosomes
            // Utilizing block scores allow for:
            // 1) Return broad peaks in case of broad modifications even for strict FDR settings
            // 2) Mitigate the problem when the number of peaks for strict FDR is much bigger than for relaxed FDR
            peaks = genomeQuery.get().map { chromosome ->
                cancellableState?.checkCanceled()
                if (!fitInfo.containsChromosomeInfo(chromosome)) {
                    return@map emptyList()
                }
                val candidates = candidateCoordinatesMap[chromosome]
                val (start, end) = chromosomeLocationsOffsets[chromosome]!!
                val logPVals = genomeLogPVals.slice(start, end)
                val logQVals = genomeLogQVals.slice(start, end)
                candidates.mapIndexedNotNull { i, candidate ->
                    val lnFdr = ln(fdr)
                    val logPValue = logPVals[i]
                    val logQValue = logQVals[i]
                    if (logPValue > lnFdr || logQValue > lnFdr) {
                        return@mapIndexedNotNull null
                    }
                    val blocks = candidate.blocks
                    val vals = blocks.map {
                        peakScorer.valueScore(ChromosomeRange(it.startOffset, it.endOffset, chromosome))
                    }
                    val peakValue = lengthWeightedScores(blocks, vals)
                    Peak(
                        chromosome = chromosome,
                        startOffset = candidate.startOffset,
                        endOffset = candidate.endOffset,
                        mlogpvalue = -logPValue / LOG_10,
                        mlogqvalue = -logQValue / LOG_10,
                        value = peakValue,
                        // Score should be proportional original q-value
                        score = min(1000.0, -logQValue / LOG_10).toInt()
                    )
                }
            }.flatten()
        } else {
            peaks = emptyList()
        }
        return peaks
    }


    fun getLogNullPvals(
        genomeQuery: GenomeQuery,
        omnipeakFitResults: OmnipeakFitResults,
        blackList: LocationsMergingList?
    ): DoubleArray {
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { omnipeakFitResults.logNullMemberships.containsKey(it.name) }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val bin = omnipeakFitResults.fitInfo.binSize
        val totalBins = limitedQuery.get()
            .sumOf { omnipeakFitResults.logNullMemberships[it.name]!!.f64Array(OmnipeakModelFitExperiment.NULL).length }
        val result = DoubleArray(totalBins) { 0.0 }
        var i = 0
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            val logNullPeps =
                omnipeakFitResults.logNullMemberships[chr.name]!!.f64Array(OmnipeakModelFitExperiment.NULL)
            for (j in 0 until logNullPeps.length) {
                val start = j * bin
                val end = (j + 1) * bin
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                result[i++] = logNullPeps[j]
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, totalBins)
        }
        return result
    }

    fun getChromosomeCandidates(
        chromosome: Chromosome,
        logNullMemberships: F64Array,
        bitList2reuse: BitList?,
        sensitivity: Double,
        gap: Int
    ): List<Range> {
        if ('_' in chromosome.name ||
            "random" in chromosome.name.lowercase() ||
            "un" in chromosome.name.lowercase()
        ) {
            LOG.trace("Ignore ${chromosome.name}: chromosome name looks like contig")
            return emptyList()
        }
        check(bitList2reuse == null || logNullMemberships.length == bitList2reuse.size())
        val bins = if (bitList2reuse != null)
            bitList2reuse
        else
            BitList(logNullMemberships.length)
        for (i in 0 until logNullMemberships.length) {
            bins.set(i, logNullMemberships[i] <= sensitivity)
        }
        return bins.aggregate(gap)
    }


    fun getChromosomeSummits(
        chromosome: Chromosome,
        fitInfo: OmnipeakFitInformation,
        candidates: List<Range>,
        avgModeLen: Int,
    ): List<Range> {
        val data = fitInfo.dataQuery.apply(chromosome)
        val scores = data.sliceAsInt(data.labels.first())
        val scoresDouble = DoubleArray(scores.size) { scores[it].toDouble() }
        val workArray = DoubleArray(avgModeLen * 100)
        val result = arrayListOf<Range>()
        candidates.forEach { candidate ->
            val start = candidate.startOffset
            val end = candidate.endOffset
            val summits = KdeModeFinder.detectModes(
                scoresDouble, avgModeLen, start, end, workArray = workArray
            )
            if (summits.isNotEmpty()) {
                summits.forEach {
                    result.add(Range(start + it.startOffset, start + it.endOffset))
                }
            } else {
                result.add(candidate)
            }
        }
        return result
    }


    fun getLogNulls(
        omnipeakFitResults: OmnipeakFitResults,
        chromosome: Chromosome
    ) = omnipeakFitResults.logNullMemberships[chromosome.name]!!.f64Array(OmnipeakModelFitExperiment.NULL)

    data class CandidateCoordinates(val startOffset: Int, val endOffset: Int, val blocks: List<Range>)

    fun processCandidateCoordinates(
        chromosome: Chromosome,
        candidates: List<Range>,
        fitInfo: OmnipeakFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        blackList: LocationsMergingList?,
        avgSignalDensity: Double?,
        avgNoiseDensity: Double?,
        clip: Double,
        summits: Boolean,
        cancellableState: CancellableState?
    ): List<CandidateCoordinates> {
        // Compute candidates coordinates and blocks to estimate significance
        // 1) Return broad peaks in case of broad modifications even for strict FDR settings
        // 2) Mitigate the problem when the number of peaks for strict FDR is much bigger than for relaxed FDR
        if (candidates.isEmpty()) {
            return emptyList()
        }
        val canEstimateScore = fitInfo is OmnipeakAnalyzeFitInformation &&
                fitInfo.binnedCoverageQueries?.all { it.areCachesPresent() } ?: false
        return candidates.mapNotNull { (from, to) ->
            cancellableState?.checkCanceled()
            var start = offsets[from]
            var end = if (to < offsets.size) offsets[to] else chromosome.length
            if (blackList != null && blackList.intersects(Location(start, end, chromosome))) {
                return@mapNotNull null
            }
            if (canEstimateScore && avgSignalDensity != null && avgNoiseDensity != null) {
                val (clippedStart, clippedEnd) = clipPeakBySignal(
                    chromosome, start, end, fitInfo, avgSignalDensity, avgNoiseDensity, clipSignal = clip,
                )
                start = clippedStart
                end = clippedEnd
            }
            var blocks =
                if (summits)
                    listOf(Range(start, end))
                else
                    // Compute significant blocks within candidate
                    candidateScoreBlocks(logNullMemberships, from, to).mapNotNull { (blockFrom, blockTo) ->
                        val blockStart = max(offsets[blockFrom], start)
                        val blockEnd = min(if (blockTo < offsets.size) offsets[blockTo] else chromosome.length, end)
                        return@mapNotNull if (blockStart < blockEnd) {
                            Range(blockStart, blockEnd)
                        } else
                            null
                    }
            if (blocks.isEmpty()) {
                blocks = listOf(Range(start, end))
            }
            return@mapNotNull CandidateCoordinates(start, end, blocks)
        }
    }


    /**
     * Using blocks for scores computation allow to hold true two invariants:
     * 1) The stricter FDR, the fewer peaks with smaller average length
     * 2) Peaks should not disappear when relaxing FDR
     * Peak score is computed as length-weighted average p-value in its consequent enriched bins.
     */
    fun estimateCandidatesLogPs(
        chromosome: Chromosome,
        candidateCoordinates: List<CandidateCoordinates>,
        peakScorer: PeakScorer,
        cancellableState: CancellableState?
    ): F64Array {
        if (candidateCoordinates.isEmpty()) {
            return FAKE_ARRAY
        }
        return F64Array(candidateCoordinates.size) { idx ->
            cancellableState?.checkCanceled()
            val blocks = candidateCoordinates[idx].blocks
            val blocksLogPs = blocks.map { (blockStart, blockEnd) ->
                val logPValue = peakScorer.logPValue(ChromosomeRange(blockStart, blockEnd, chromosome))
                check(!logPValue.isNaN()) { "P-value is nan" }
                return@map logPValue
            }
            return@F64Array lengthWeightedScores(blocks, blocksLogPs)
        }
    }

    /**
     * Compute candidate score blocks to get final score of the peak
     */
    fun candidateScoreBlocks(logNullMemberships: F64Array, start: Int, end: Int): List<Range> {
        val p = StatUtils.percentile(
            logNullMemberships.slice(start, end).toDoubleArray(), OMNIPEAK_SCORE_BLOCKS * 100
        )
        return BitList(end - start) { logNullMemberships[it + start] <= p }
            .aggregate(3).map { Range(start + it.startOffset, start + it.endOffset) }
    }

    fun estimateSingleModeLength(
        genomeQuery: GenomeQuery,
        fitInfo: OmnipeakFitInformation,
        candidates: GenomeMap<List<Range>>,
        parallel: Boolean,
        dipTestN: Int = 100,
        dipTestPValue: Double = 0.05,
        estimateNumber: Int = 100
    ): Int {
        val singleModeCandidates = AtomicInteger(0)
        val singleModeLength = AtomicInteger(0)

        genomeQuery.get().map { chromosome ->
            Callable {
                // We don't need every single mode candidate for estimation
                if (singleModeCandidates.get() >= estimateNumber) {
                    return@Callable
                }
                if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                    return@Callable
                }
                if (!fitInfo.scoresAvailable()) {
                    return@Callable
                }
                val data = fitInfo.dataQuery.apply(chromosome)
                val scores = data.sliceAsInt(data.labels.first())

                val chrCandidates = candidates[chromosome]
                chrCandidates.forEach { (from, to) ->
                    // We don't need every single mode candidate for estimation
                    if (singleModeCandidates.get() >= estimateNumber) {
                        return@Callable
                    }
                    val slice = scores.slice(from until to)
                    val signal = DoubleArray(to - from) { slice[it].toDouble() }
                    // Ignore multiple test adjustment here
                    if (HartiganDip.dipTest(signal, dipTestN).pValue > dipTestPValue) {
                        singleModeCandidates.incrementAndGet()
                        singleModeLength.addAndGet(to - from)
                    }
                }
            }
        }.await(parallel)
        val n = singleModeCandidates.get()
        return if (n == 0) 0 else singleModeLength.get() / n
    }


    /**
     * Summary length-weighted score for a peak.
     * Use length weighted mean to take into account the difference in block lengths
     */
    fun lengthWeightedScores(
        blocks: List<Range>,
        scores: List<Double>,
    ): Double {
        require(blocks.size == scores.size) { "Different lengths of blocks and scores lists" }
        if (blocks.isEmpty()) {
            println()
        }
        require(blocks.isNotEmpty()) { "Empty blocks list" }
        if (blocks.size == 1) {
            return scores.first()
        }
        val sum = KahanSum()
        var l = 0L
        blocks.zip(scores).sortedBy { it.second }
            .forEach { (b, p) ->
                val blockLen = b.endOffset - b.startOffset
                check(blockLen > 0) { "Not positive block length" }
                sum += p * blockLen
                l += blockLen
            }
        val value = sum.result() / l
        check(l > 0) { "Empty summary length" }
        check(!value.isNaN()) { "Length weighted mean is nan" }
        return value
    }

    val LOG_10 = ln(10.0)

    private val FAKE_ARRAY = F64Array.full(1, 0.0)

}
