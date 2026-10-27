package org.jetbrains.bio.omnipeak.peaks

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.omnipeak.OmnipeakBridgeImpl
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.getChromosomeCandidates
import org.jetbrains.bio.omnipeak.peaks.CandidatesScores.Companion.effectiveBins
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.SensitivityGapResult
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.math.max
import kotlin.math.min
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CandidatesTest {

    @Test
    fun testSensitivityGapResultSensitivity() {
        val res = SensitivityGapResult(sensitivity = -2.5, gap = 5)
        assertEquals(-2.5, res.sensitivity)
        assertEquals(5, res.gap)
    }

    @Test
    fun testGetChromosomeCandidates() {
        val chr = Genome["to1"].chromosomes.first()
        val logNulls = F64Array(20) { i ->
            when (i) {
                in 2..4 -> -5.0
                in 5..6 -> -1.0
                in 7..9 -> -5.0
                else -> 0.0
            }
        }
        val bitList = BitList(20)
        val candidates = getChromosomeCandidates(chr, logNulls, bitList, -0.5, 0)
        assertEquals(listOf(Range(2, 10)), candidates)

        val strictCandidates = getChromosomeCandidates(chr, logNulls, bitList, -3.0, 0)
        assertEquals(listOf(Range(2, 5), Range(7, 10)), strictCandidates)

        val strictCandidatesMerged = getChromosomeCandidates(chr, logNulls, bitList, -3.0, 2)
        assertEquals(listOf(Range(2, 10)), strictCandidatesMerged)
    }

    @Test
    fun testEffectiveBinsCalculation() {
        assertEquals(0.0, effectiveBins(0, 0.5), 1e-6)
        // 1 bin always has effective bins = 1.0
        assertEquals(1.0, effectiveBins(1, 0.0), 1e-6)
        assertEquals(1.0, effectiveBins(1, 0.5), 1e-6)
        assertEquals(1.0, effectiveBins(1, 0.8), 1e-6)
        assertEquals(1.0, effectiveBins(1, 1.0), 1e-6)

        // rho = 0.0 (independent): linear scaling N_eff = 1 + 4 * 1.0 = 5.0
        assertEquals(5.0, effectiveBins(5, 0.0), 1e-6)
        // rho = 0.5: N_eff = 1 + 4 * 0.5 = 3.0
        assertEquals(3.0, effectiveBins(5, 0.5), 1e-6)
        // rho = 0.8: N_eff = 1 + 4 * 0.2 = 1.8
        assertEquals(1.8, effectiveBins(5, 0.8), 1e-6)
        // rho = 1.0 (fully correlated): N_eff = 1.0
        assertEquals(1.0, effectiveBins(5, 1.0), 1e-6)

        // Clamping checks
        assertEquals(5.0, effectiveBins(5, -0.5), 1e-6)
        assertEquals(1.0, effectiveBins(5, 1.5), 1e-6)
    }

    @Test
    fun testPeakScorerModelScoreFallback() {
        val chr = Genome["to1"].chromosomes.first()
        val gq = GenomeQuery(chr.genome, chr.name)
        val fitInfo = OmnipeakAnalyzeFitInformation(
            gq, emptyList(), emptyList(), null, AutoFragment, true, 200, true
        )
        val logNulls: F64Array = F64Array.of(-10.0, -20.0, -30.0, -40.0, -50.0)
        val logNullMembershipsMap = genomeMap<F64Array>(gq) { logNulls }

        val cr1 = ChromosomeRange(0, 200, chr) // 1 bin: index 0 -> -10.0
        val cr3 = ChromosomeRange(200, 800, chr) // 3 bins: indices 1, 2, 3 -> (-20 + -30 + -40)/3 = -30.0

        val scorer = CandidatesScores.create(fitInfo, logNullMembershipsMap, avgAutocorrelation = 0.5)
        assertEquals(-10.0, scorer.logP(cr1), 1e-6)
        assertEquals(-60.0, scorer.logP(cr3), 1e-6)
        assertEquals(10.0, scorer.valueScore(cr1), 1e-6)
        assertEquals(60.0, scorer.valueScore(cr3), 1e-6)
    }

    @Test
    fun testComputeMergedSensitivityNarrowMarks() {
        val minAdditionalSensitivity = -5.0
        val logFdr = -20.0 // stricter FDR
        // For narrow marks (avgAutocorrelation = 0.0), merged sensitivity matches FDR
        val merged = SensitivityGap.computeMergedSensitivity(
            0.0, minAdditionalSensitivity, logFdr
        )
        assertEquals(-20.0, merged, 1e-6)
    }

    @Test
    fun testComputeMergedSensitivityBroadMarks() {
        val minAdditionalSensitivity = -5.0
        val logFdr = -20.0
        // For broad marks (avgAutocorrelation = 1.0), merged sensitivity matches minAdditionalSensitivity
        val merged = SensitivityGap.computeMergedSensitivity(
            1.0, minAdditionalSensitivity, logFdr
        )
        assertEquals(-5.0, merged, 1e-6)
    }

    @Test
    fun testComputeMergedSensitivityIntermediateScaling() {
        val minAdditionalSensitivity = -5.0
        val logFdr = -25.0
        // Intermediate autocorrelation = 0.5: -25.0 + (-5.0 - (-25.0)) * 0.5 = -15.0
        val merged = SensitivityGap.computeMergedSensitivity(
            0.5, minAdditionalSensitivity, logFdr
        )
        assertEquals(-15.0, merged, 1e-6)
    }

    @Test
    fun testComputeMergedSensitivityLogFdrMorePermissive() {
        val minAdditionalSensitivity = -5.0
        val logFdr = -2.0 // more permissive than minAdditionalSensitivity
        val merged = SensitivityGap.computeMergedSensitivity(
            0.0, minAdditionalSensitivity, logFdr
        )
        assertEquals(-5.0, merged, 1e-6)
    }

    @Test
    fun testComputeMergedSensitivityEdgeCases() {
        assertEquals(-5.0, SensitivityGap.computeMergedSensitivity(0.5, null, -5.0), 1e-6)
    }

    @Test
    fun testInterpolateSensitivity() {
        // avgAutocorrelation = 0.0 -> stricterSensitivity
        assertEquals(-20.0, SensitivityGap.interpolateSensitivity(0.0, -20.0, -5.0), 1e-6)
        // avgAutocorrelation = 1.0 -> stableSensitivity
        assertEquals(-5.0, SensitivityGap.interpolateSensitivity(1.0, -20.0, -5.0), 1e-6)
        // avgAutocorrelation = 0.5 -> intermediate value -12.5
        assertEquals(-12.5, SensitivityGap.interpolateSensitivity(0.5, -20.0, -5.0), 1e-6)
    }

    @Test
    fun testPeakScorerHybridLocalBackgroundAutocorrelation() {
        val chr = Genome["to1"].chromosomes.first()
        val cr = ChromosomeRange(5000, 5200, chr) // 200 bp candidate

        // Setup coverage computable: 50 reads in peak [5000, 5200), background outside
        val scorerLowNoise = CandidatesScores(
            coverageControlComputable = null,
            coverageComputable = { range ->
                // Return 50 inside peak, 0.001 reads/bp outside
                val peakOverlap = max(0, min(range.endOffset, 5200) - max(range.startOffset, 5000))
                val nonPeakLength = (range.endOffset - range.startOffset) - peakOverlap
                (peakOverlap / 200.0) * 50.0 + nonPeakLength * 0.001
            },
            binSize = 200,
            avgAutocorrelation = 0.5,
            failures = 5.0
        ) { 0.0 }

        val scorerHighNoise = CandidatesScores(
            coverageControlComputable = null,
            coverageComputable = { range ->
                // Return 50 inside peak, 0.2 reads/bp outside (high background)
                val peakOverlap = max(0, min(range.endOffset, 5200) - max(range.startOffset, 5000))
                val nonPeakLength = (range.endOffset - range.startOffset) - peakOverlap
                (peakOverlap / 200.0) * 50.0 + nonPeakLength * 0.2
            },
            binSize = 200,
            avgAutocorrelation = 0.5,
            failures = 5.0
        ) { 0.0 }

        val logPLow = scorerLowNoise.logP(cr)
        val logPHigh = scorerHighNoise.logP(cr)

        assertTrue(logPLow < -20.0, "Low noise background should yield highly significant logP: $logPLow")
        assertTrue(logPHigh > logPLow, "High local background noise must penalize significance: low=$logPLow, high=$logPHigh")
    }

    @Test
    fun testPeakScorerNegativeBinomialVsPoisson() {
        val chr = Genome["to1"].chromosomes.first()
        val cr = ChromosomeRange(1000, 1200, chr)

        // 20 reads in peak, expected background rate ~ 5 reads
        val scorerPoisson = CandidatesScores(
            coverageControlComputable = null,
            coverageComputable = { range ->
                if (range.startOffset == 1000 && range.endOffset == 1200) 20.0 else 5.0 * ((range.endOffset - range.startOffset) / 200.0)
            },
            binSize = 200,
            avgNoiseDensity = 5.0 / 200.0
        ) { 0.0 } // Poisson


        val scorerNegBin = CandidatesScores(
            coverageControlComputable = null,
            coverageComputable = { range ->
                if (range.startOffset == 1000 && range.endOffset == 1200) 20.0 else 5.0 * ((range.endOffset - range.startOffset) / 200.0)
            },
            binSize = 200,
            avgNoiseDensity = 5.0 / 200.0,
            failures = 2.0
        ) { 0.0 } // Overdispersed


        val logPPoisson = scorerPoisson.logP(cr)
        val logPNegBin = scorerNegBin.logP(cr)

        // Overdispersion accounts for extra variance, leading to conservative log p-value
        assertTrue(logPNegBin > logPPoisson, "Negative binomial with overdispersion should be more conservative than Poisson: NB=$logPNegBin, Poisson=$logPPoisson")
    }

    @Test
    fun testBridgeCandidatesToPeaksPropagatesFailures() {
        val chr = Genome["to1"].chromosomes.first()
        val gq = GenomeQuery(chr.genome, chr.name)
        val fitInfo = OmnipeakAnalyzeFitInformation(
            gq, emptyList(), emptyList(), null, AutoFragment, true, 200, true
        )
        val logNulls = F64Array(10) { if (it in 2..4) -20.0 else 0.0 }
        val logNullMap = genomeMap(gq) { logNulls }
        val candidates = genomeMap(gq) { listOf(Range(2, 5)) }
        val bridge = OmnipeakBridgeImpl()

        val peaksPoisson = bridge.candidatesToPeaks(
            gq, fitInfo, logNullMap, candidates,
            sensitivity = -10.0, avgAutocorrelation = 0.0, summits = false, singleModeLength = 0,
            clip = 0.0, fdr = 0.5, failures = null, multipleTesting = MultipleTesting.BH,
            blacklistPath = null, name = "test", parallel = false, cancellableState = CancellableState.fake()
        )
        val peaksNegBin = bridge.candidatesToPeaks(
            gq, fitInfo, logNullMap, candidates,
            sensitivity = -10.0, avgAutocorrelation = 0.0, summits = false, singleModeLength = 0,
            clip = 0.0, fdr = 0.5, failures = 1.0, multipleTesting = MultipleTesting.BH,
            blacklistPath = null, name = "test", parallel = false, cancellableState = CancellableState.fake()
        )
        assertEquals(1, peaksPoisson.size)
        assertEquals(1, peaksNegBin.size)
    }

    @Test
    fun testPeakScorerControlCoverageSingleQueryAndPerformance() {
        val chr = Genome["to1"].chromosomes.first()
        val cr = ChromosomeRange(1000, 1200, chr)
        var queryCount = 0

        val scorer = CandidatesScores(
            coverageControlComputable = { range ->
                queryCount++
                if (range.startOffset == 1000 && range.endOffset == 1200) {
                    25.0 to 2.0
                } else {
                    val scale = (range.endOffset - range.startOffset) / 200.0
                    (5.0 * scale) to (2.0 * scale)
                }
            },
            coverageComputable = null,
            binSize = 200,
            avgAutocorrelation = 0.5,
            failures = 2.0
        ) { 0.0 }

        val logP = scorer.logP(cr)
        assertTrue(logP < -5.0, "LogP should be significant: $logP")
        // cr query + 1 local window query
        assertEquals(2, queryCount, "Expected exactly 2 queries (candidate range + local window)")

        // Benchmark speed for 10_000 calls
        val iterations = 10_000
        val startTime = System.nanoTime()
        for (i in 0 until iterations) {
            scorer.logP(cr)
        }
        val elapsedMs = (System.nanoTime() - startTime) / 1_000_000.0
        assertTrue(elapsedMs < 500.0, "10k logP calls took ${elapsedMs}ms, should complete in <500ms")
    }

    @Test
    fun testPeakScorerMultiBinNegativeBinomialScaling() {
        val chr = Genome["to1"].chromosomes.first()
        // Broad domain spanning 50 bins (10,000 bp at binSize = 200)
        val cr = ChromosomeRange(10_000, 20_000, chr)

        // 2x fold enrichment: 50 bins * 10 reads/bin = 500 reads vs background rate of 5 reads/bin = 250 reads
        val scorerBroad = CandidatesScores(
            coverageControlComputable = null,
            coverageComputable = { range ->
                val overlap = max(0, min(range.endOffset, 20_000) - max(range.startOffset, 10_000))
                val nonOverlap = (range.endOffset - range.startOffset) - overlap
                (overlap / 200.0) * 10.0 + (nonOverlap / 200.0) * 5.0
            },
            binSize = 200,
            avgNoiseDensity = 5.0 / 200.0,
            avgAutocorrelation = 0.5,
            failures = 2.0
        ) { 0.0 }

        val logP = scorerBroad.logP(cr)
        // With correct N^2 / N_eff dispersion scaling, 500 reads vs 250 expected in a 10kb broad domain must be highly significant
        assertTrue(logP < -10.0, "Broad 10kb domain with 2x enrichment should be highly significant, got logP = $logP")
    }

    @Test
    fun testPeakScorerBroadPeakNearSharpSpikeWithControl() {
        val chr = Genome["to1"].chromosomes.first()
        // Broad domain [1000, 3000) (2000 bp) with moderate 3x enrichment (60 reads vs 20 control)
        val cr = ChromosomeRange(1000, 3000, chr)

        // Nearby sharp spike at [3200, 3400) with 500 reads in treatment, but control is flat (2 reads)
        val scorer = CandidatesScores(
            coverageControlComputable = { range ->
                var treat = 0.0
                var ctrl = 0.0
                // Candidate overlap [1000, 3000)
                val candOverlap = max(0, min(range.endOffset, 3000) - max(range.startOffset, 1000))
                treat += (candOverlap / 2000.0) * 60.0
                ctrl += (candOverlap / 2000.0) * 20.0

                // Spike overlap [3200, 3400)
                val spikeOverlap = max(0, min(range.endOffset, 3400) - max(range.startOffset, 3200))
                treat += (spikeOverlap / 200.0) * 500.0
                ctrl += (spikeOverlap / 200.0) * 2.0

                // Background
                val otherLength = (range.endOffset - range.startOffset) - candOverlap - spikeOverlap
                treat += otherLength * 0.01
                ctrl += otherLength * 0.01

                treat to ctrl
            },
            coverageComputable = null,
            binSize = 200,
            avgNoiseDensity = 0.01,
            avgAutocorrelation = 0.5,
            failures = 5.0
        ) { 0.0 }

        val logP = scorer.logP(cr)
        // The broad peak has 60 reads vs 20 control reads; because control is used (ignoring treatment spike in flank),
        // the candidate must remain statistically significant.
        assertTrue(logP < -5.0, "Broad peak near sharp spike with control should remain significant, got logP = $logP")
    }

    @Test
    fun testPeakScorerBroadPeakNearSharpSpikeWithoutControl() {
        val chr = Genome["to1"].chromosomes.first()
        // Broad domain [1000, 3000) (2000 bp) with 80 reads vs global noise density 0.01 (20 expected reads)
        val cr = ChromosomeRange(1000, 3000, chr)

        // Nearby sharp spike at [3200, 3400) with 1000 reads in treatment
        val scorer = CandidatesScores(
            coverageControlComputable = null,
            coverageComputable = { range ->
                var treat = 0.0
                val candOverlap = max(0, min(range.endOffset, 3000) - max(range.startOffset, 1000))
                treat += (candOverlap / 2000.0) * 80.0

                val spikeOverlap = max(0, min(range.endOffset, 3400) - max(range.startOffset, 3200))
                treat += (spikeOverlap / 200.0) * 1000.0

                val otherLength = (range.endOffset - range.startOffset) - candOverlap - spikeOverlap
                treat += otherLength * 0.01

                treat
            },
            binSize = 200,
            avgNoiseDensity = 0.01,
            avgAutocorrelation = 0.5,
            failures = 5.0
        ) { 0.0 }

        val logP = scorer.logP(cr)
        // Bounding local treatment noise to 2x global noise density prevents the 1000-read spike from masking the broad peak
        assertTrue(logP < -4.0, "Broad peak near sharp spike without control should remain significant, got logP = $logP")
    }
}
