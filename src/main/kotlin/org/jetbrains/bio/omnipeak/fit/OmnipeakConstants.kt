package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.omnipeak.OmnipeakCLA
import org.jetbrains.bio.omnipeak.peaks.MultipleTesting
import org.jetbrains.bio.viktor.F64Array
import java.lang.reflect.Modifier
import kotlin.math.ln

/**
 * Constants used in Omnipeak.
 * Those, configurable from command line interface, have "DEFAULT" in their names.
 */
object OmnipeakConstants {
    /**
     * Default bin size used for binned centered reads coverage aggregation.
     */
    const val OMNIPEAK_DEFAULT_BIN = 100

    /**
     * The default variable represents the default value for the FDR parameter.
     * The FDR is a measure of the expected proportion of false positive peaks in the peak calling results.
     */
    const val OMNIPEAK_DEFAULT_FDR = 0.05

    val OMNIPEAK_DEFAULT_MULTIPLE_TEST_CORRECTION = MultipleTesting.BH

    /**
     * The default step size used for beta values in the Omnipeak calculation.
     * Beta value is used to minimize correlation between treatment and control track.
     * Step is used to iterate in the [0, 1] interval.
     */
    const val OMNIPEAK_BETA_STEP = 0.01

    // Max shift to compute auto correlations
    const val OMNIPEAK_AUTOCORRELATION_MAX_SHIFT = 50

    // Technical minimal coefficient between variance and mean of Negative Binomials
    const val OMNIPEAK_HMM_NB_VAR_MEAN_MULTIPLIER = 1 + 1e-3

    // Fraction scores used for HMM signal, noise and ratio estimation, guards decent SNR in model
    const val OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR = 0.1

    // Fraction scores used for HMM noise estimation
    const val OMNIPEAK_HMM_ESTIMATE_LOW = 0.5

    // Minimal low state mean threshold, guards against too broad peaks
    const val OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD = 0.3

    // Technical threshold to limit mean to std, guards against artificial data without noise
    const val OMNIPEAK_HMM_MAX_MEAN_TO_STD = 5.0

    // General model priors and priors based on real data peaks footprint
    val OMNIPEAK_HMM_PRIORS = F64Array.of(0.75, 0.249, 0.001)

    val OMNIPEAK_HMM_TRANSITIONS = listOf(
        doubleArrayOf(0.75, 0.2499, 0.0001),
        doubleArrayOf(0.2, 0.798, 0.002),
        doubleArrayOf(0.005, 0.015, 0.98))

    /**
     * The threshold  used in the Omnipeak peak-fitting algorithm,
     * relative delta between the model log likelihood
     */
    const val OMNIPEAK_DEFAULT_FIT_THRESHOLD = 1e-4

    /**
     * The maximum number of iterations to perform during the fitting
     */
    const val OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS = 10

    /**
     * Number of points between relaxed and strict sensitivity to analyse
     */
    const val OMNIPEAK_SENSITIVITY_N = 100

    /**
     * Limit min sensitivity during candidates saturation analysis
     */
    const val OMNIPEAK_MIN_SENSITIVITY = -1e-10

    /**
     * Fraction of top significant bins to be used for peak score estimation.
     * It should be robust wrt appending blocks of low significance,
     * so take into account top N% bins into block, otherwise we'll get fdr-blinking peaks, i.e.,
     * peaks which are present for stronger fdr, but missing for more relaxed settings
     */
    const val OMNIPEAK_SCORE_BLOCKS = 0.5

    const val OMNIPEAK_SCORE_BLOCKS_GAP = 3

    const val OMNIPEAK_DEFAULT_GAP = 0

    val OMNIPEAK_DEFAULT_SENSITIVITY = ln(OMNIPEAK_DEFAULT_FDR)

    // Max gap to compute fragmentation,
    // i.e. reduction of candidate number when merging with gap
    const val OMNIPEAK_FRAGMENTATION_MAX_GAP_BP = 5000

    // Rule of thumb: max when narrow marks and ATAC-seq data are not fragmented
    // Fragmentation score is an area above the curve of relative candidates number by gap
    const val OMNIPEAK_DEFAULT_FRAGMENTATION_THRESHOLD_BP = 500

    // When calling summits, min summit length
    const val OMNIPEAK_SUMMITS_MIN_LENGTH = 3

    // When calling summits, minimal relative length between summits to merge
    const val OMNIPEAK_SUMMITS_MIN_DISTANCE = 2

    /**
     * Clipping allows to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    const val OMNIPEAK_DEFAULT_CLIP_MAX_SIGNAL = 0.4

    const val OMNIPEAK_CLIP_MAX_LENGTH = 0.8

    val OMNIPEAK_CLIP_STEPS = doubleArrayOf(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0)

    fun printConstants() {
        for (field in OmnipeakConstants::class.java.declaredFields) {
            // Ignore singleton and fields, configured with params
            if (field.name == "INSTANCE" || "DEFAULT" in field.name) {
                continue
            }
            if (Modifier.isStatic(field.modifiers) && Modifier.isFinal(field.modifiers)) {
                field.isAccessible = true
                val value = field.get(null) // null because it's a static field
                OmnipeakCLA.LOG.debug("${field.name}: ${pp(value)}")
            }
        }
    }

    fun pp(value: Any?): String {
        return when (value) {
            is List<*> -> "[${value.joinToString(", ") { pp(it) }}]"
            is IntArray -> value.contentToString()
            is DoubleArray -> value.contentToString()
            is Array<*> -> value.contentDeepToString()
            else -> value.toString()
        }
    }
}