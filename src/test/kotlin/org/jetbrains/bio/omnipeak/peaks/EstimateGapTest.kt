package org.jetbrains.bio.omnipeak.peaks

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_BIN
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_FRAGMENTATION_MAX_GAP_BP
import org.jetbrains.bio.omnipeak.peaks.SensitivityGap.estimateGap
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Tests for [SensitivityGap.estimateGap] on real candidates-by-gap curves (the `*.peak.gaps.tsv`
 * diagnostics). The estimator combines normalized log-space convexity with the spatial
 * autocorrelation correlation length scale prior.
 * The arrays below are the verbatim `CandidatesN` columns (gap 0..49).
 */
class EstimateGapTest {

    /**
     * OD10 H3K27ac, narrow: gentle ~2x decline across the range.
     * Data source: https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/H3K27ac/OD10_k27ac_hg19.bed.gz
     */
    private val od10K27ac = intArrayOf(
        68936, 68607, 68166, 67578, 66843, 65919, 64849, 63677, 62596, 61389,
        60207, 59022, 57886, 56727, 55652, 54600, 53588, 52564, 51628, 50734,
        49841, 49045, 48237, 47525, 46806, 46087, 45372, 44719, 44095, 43479,
        42898, 42295, 41698, 41150, 40596, 40128, 39622, 39108, 38690, 38221,
        37763, 37354, 36913, 36519, 36101, 35681, 35304, 34914, 34529, 34157
    )

    private val od10K27acChr1 = intArrayOf(
        1588, 1581, 1573, 1563, 1542, 1510, 1485, 1448, 1412, 1378, 1338, 1293, 1261, 1238,
        1204, 1179, 1153, 1120, 1093, 1070, 1049, 1027, 1005, 986, 961, 939, 923, 899, 886,
        872, 853, 847, 832, 821, 803, 789, 779, 769, 754, 745, 736, 724, 709, 699, 687, 678, 666, 654, 641, 633)

    /**
     * OD5 H3K27ac, narrow: very gentle ~1.5x decline.
     * Data source: https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/H3K27ac/OD5_k27ac_hg19.bed.gz
     */
    private val od5K27ac = intArrayOf(
        35144, 35012, 34826, 34620, 34350, 34043, 33691, 33329, 32983, 32653,
        32321, 32008, 31707, 31423, 31108, 30808, 30546, 30276, 29996, 29768,
        29470, 29200, 28936, 28683, 28439, 28187, 27930, 27704, 27491, 27274,
        27041, 26864, 26641, 26442, 26226, 26034, 25856, 25659, 25475, 25290,
        25120, 24925, 24768, 24589, 24403, 24240, 24073, 23919, 23794, 23640
    )

    /**
     * OD1 H3K36me3, broad: steep, convex ~7.7x decline.
     * Data source: https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/H3K36me3/OD1_k36me3_hg19.bed.gz
     */
    private val od1K36me3 = intArrayOf(
        171132, 154526, 131877, 116107, 103897, 93968, 85671, 79048, 73514, 68917,
        64781, 61220, 58001, 55234, 52697, 50431, 48354, 46447, 44796, 43240,
        41828, 40470, 39253, 38010, 36911, 35874, 34944, 34028, 33194, 32414,
        31605, 30891, 30197, 29547, 28938, 28353, 27789, 27231, 26756, 26258,
        25746, 25280, 24892, 24491, 24102, 23676, 23304, 22911, 22581, 22228
    )

    /**
     * OD4 H3K36me3, broad: convex ~4x decline.
     * Data source: https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/H3K36me3/OD4_k36me3_hg19.bed.gz
     */
    private val od4K36me3 = intArrayOf(
        71755, 71030, 68531, 64890, 61110, 57543, 54239, 51403, 48854, 46615,
        44562, 42638, 40945, 39345, 37854, 36554, 35285, 34199, 33165, 32249,
        31296, 30423, 29671, 28977, 28244, 27601, 27001, 26399, 25784, 25229,
        24715, 24223, 23739, 23287, 22863, 22480, 22046, 21632, 21260, 20889,
        20537, 20217, 19906, 19596, 19297, 18994, 18726, 18475, 18236, 17978
    )

    /**
     * E065 H3K36me3, broad but bad quality: merges a lot overall (278k -> 72k) yet gently and
     * persistently - the per-step reduction dips below 5% early (~gap 7) while the curve never
     * flattens. It must get a large gap, not the false convergence at 7.
     * Data source: https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/E065-H3K36me3.tagAlign.gz
     */
    private val e065K36me3 = intArrayOf(
        278675, 268204, 248990, 230268, 214989, 202081, 190770, 181337, 172883, 165280,
        158522, 152393, 146924, 142102, 137513, 133339, 129461, 125894, 122613, 119472,
        116572, 113861, 111327, 108879, 106572, 104292, 102259, 100257, 98362, 96670,
        94923, 93292, 91742, 90285, 88825, 87423, 86101, 84798, 83592, 82384,
        81237, 80086, 78998, 78015, 77041, 76037, 75105, 74238, 73308, 72422
    )

    /**
     * GM12878 H3K27ac chr1 (100 bp bin): punctate mark with clustered active enhancers.
     * Evaluated across gaps 0..100.
     */
    private val gm12878K27ac = intArrayOf(
        4974, 4961, 4931, 4889, 4834, 4761, 4660, 4550, 4456, 4375,
        4305, 4215, 4134, 4062, 3997, 3938, 3889, 3842, 3797, 3752,
        3705, 3661, 3614, 3565, 3527, 3492, 3455, 3426, 3393, 3351,
        3321, 3293, 3272, 3247, 3220, 3193, 3162, 3138, 3117, 3089,
        3063, 3041, 3021, 3001, 2978, 2963, 2940, 2924, 2910, 2888,
        2865, 2848, 2820, 2802, 2782, 2766, 2747, 2738, 2718, 2705,
        2689, 2673, 2653, 2643, 2623, 2611, 2606, 2591, 2569, 2551,
        2541, 2534, 2527, 2511, 2494, 2474, 2458, 2451, 2435, 2422,
        2413, 2402, 2385, 2375, 2363, 2352, 2345, 2334, 2322, 2306,
        2296, 2289, 2277, 2268, 2260, 2250, 2238, 2224, 2216, 2213, 2203
    )

    @Test
    fun narrowMarksGetSmallGap() {
        assertEquals(0, estimateGap(od10K27acChr1), "OD10 H3K27ac chr1")
        assertEquals(0, estimateGap(od10K27ac), "OD10 H3K27ac")
        assertEquals(0, estimateGap(od5K27ac), "OD5 H3K27ac")
        assertEquals(0, estimateGap(gm12878K27ac), "GM12878 H3K27ac chr1")
    }

    @Test
    fun broadMarksGetLargerGap() {
        assertEquals(10, estimateGap(od1K36me3, avgAutocorrelation = 0.8), "OD1 H3K36me3")
        assertEquals(7, estimateGap(od4K36me3, avgAutocorrelation = 0.8), "OD4 H3K36me3")
        assertEquals(9, estimateGap(e065K36me3, avgAutocorrelation = 0.8), "E065 H3K36me3")
    }

    @Test
    fun worstFragmentationGetsLargeGap() {
        // Step function: many candidates at gap 0, but they all merge at gap 1.
        val n = OMNIPEAK_FRAGMENTATION_MAX_GAP_BP / OMNIPEAK_DEFAULT_BIN
        val worst = IntArray(n) { if (it == 0) 10000 else 1 }
        val gap = estimateGap(worst, avgAutocorrelation = 1.0)
        assertEquals(n - 1, gap, "Worst case gap should be ${n - 1}")
    }

    @Test
    fun narrowGapsAreSmallerThanBroadGaps() {
        val narrow = listOf(od10K27acChr1, od10K27ac, od5K27ac).maxOf { estimateGap(it) }
        val broad = listOf(od1K36me3, od4K36me3, e065K36me3).minOf { estimateGap(it, avgAutocorrelation = 0.8) }
        assertTrue(narrow < broad, "Narrow gap $narrow should be smaller than broad gap $broad")
    }

    @Test
    fun broadAutocorrelationModulatesGapThreshold() {
        val moderateConvexity = IntArray(50) { (10000.0 / (1.0 + it * 0.01)).toInt() }
        assertEquals(0, estimateGap(moderateConvexity, avgAutocorrelation = 0.0), "Default/narrow autocorrelation returns gap 0")
        assertEquals(0, estimateGap(moderateConvexity, avgAutocorrelation = 0.1), "Low autocorrelation (0.1) returns gap 0")
        assertEquals(0, estimateGap(od10K27ac, avgAutocorrelation = 0.2), "Narrow mark with low autocorrelation returns gap 0")
        val gapBroad = estimateGap(moderateConvexity, avgAutocorrelation = 0.8)
        assertTrue(gapBroad > 0, "High autocorrelation enables gap merging: $gapBroad")
    }

    @Test
    fun testAutocorrelationScalePriorGraduation() {
        val gap0 = estimateGap(od1K36me3, avgAutocorrelation = 0.0)
        val gap2 = estimateGap(od1K36me3, avgAutocorrelation = 0.2)
        val gap5 = estimateGap(od1K36me3, avgAutocorrelation = 0.5)
        val gap8 = estimateGap(od1K36me3, avgAutocorrelation = 0.8)
        val gap9 = estimateGap(od1K36me3, avgAutocorrelation = 0.9)
        val gap10 = estimateGap(od1K36me3, avgAutocorrelation = 1.0)

        assertEquals(0, gap0)
        assertTrue(gap2 in 1..3, "Gap at rho=0.2 should be small: $gap2")
        assertTrue(gap5 in 4..7, "Gap at rho=0.5 should be moderate: $gap5")
        assertTrue(gap8 >= 10, "Gap at rho=0.8 should be large: $gap8")
        assertTrue(gap9 >= gap8, "Gap should increase monotonically with autocorrelation")
        assertTrue(gap10 >= gap9, "Gap at rho=1.0 should be maximal")
    }

    @Test
    fun testWindowLengthScaling() {
        val candidateNs50 = IntArray(50) { (10000.0 / (1.0 + it * 0.05)).toInt() }
        val candidateNs25 = IntArray(25) { (10000.0 / (1.0 + it * 0.10)).toInt() }
        val gap50 = estimateGap(candidateNs50, avgAutocorrelation = 0.8)
        val gap25 = estimateGap(candidateNs25, avgAutocorrelation = 0.8)
        assertTrue(gap50 > 0)
        assertTrue(gap25 > 0)
    }

    @Test
    fun edgeCases() {
        assertEquals(0, estimateGap(intArrayOf()), "Empty candidates")
        assertEquals(0, estimateGap(intArrayOf(100)), "Single element")
        assertEquals(0, estimateGap(intArrayOf(100, 100, 100)), "Constant candidates (no decline)")
        assertEquals(0, estimateGap(intArrayOf(100, 105, 110)), "Increasing candidates")
    }

    @Test
    fun testEstimateCandidateGapNsMatchesIterative() {
        val genomeQuery = GenomeQuery(Genome["to1"], "chr1", "chr2")
        val fitInfo = OmnipeakAnalyzeFitInformation(
            genomeQuery, emptyList(), emptyList(),
            null, AutoFragment, true, 200, true
        )
        val sensitivity = -2.0
        val logNullMembershipsMap = genomeMap(genomeQuery) { chr ->
            val len = (chr.length / 200) + 1
            F64Array(len) { (it * 37 % 100).toDouble() / -20.0 }
        }
        val bitList2reuseMap = genomeMap(genomeQuery) { chr ->
            BitList(logNullMembershipsMap[chr].length)
        }
        val maxGap = 20
        val fastNs = SensitivityGap.estimateCandidateGapNs(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
            sensitivity, maxGap = maxGap, parallel = false
        )
        val slowNs = IntArray(maxGap + 1) { g ->
            SensitivityGap.estimateCandidatesNumberLens(
                genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                sensitivity, g
            ).n
        }
        assertEquals(maxGap + 1, fastNs.size)
        for (g in 0..maxGap) {
            assertEquals(slowNs[g], fastNs[g], "Mismatch at gap $g")
        }
    }
}
