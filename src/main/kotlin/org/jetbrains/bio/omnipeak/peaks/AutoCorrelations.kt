package org.jetbrains.bio.omnipeak.peaks

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_AUTOCORRELATION_MAX_SHIFT
import java.util.*
import kotlin.math.min

object AutoCorrelations {

    fun computeAutoCorrelations(
        values: DoubleArray,
        maxCorrelationDistance: Int = OMNIPEAK_AUTOCORRELATION_MAX_SHIFT
    ): DoubleArray {
        val distanceCorrelations = DoubleArray(min(values.size / 2, maxCorrelationDistance) + 1)
        distanceCorrelations[0] = 1.0
        for (d in 1 until distanceCorrelations.size) {
            val correlation = computeAutocorrelation(values, d)
            // Stop computing correlation after small one
            if (Precision.equals(correlation, 0.0, 1e-6)) {
                break
            }
            distanceCorrelations[d] = correlation
        }
        return distanceCorrelations
    }

    private fun computeAutocorrelation(values: DoubleArray, d: Int): Double {
        val original = DoubleArray(values.size - d - 1)
        System.arraycopy(values, 0, original, 0, values.size - d - 1)
        val shifted = DoubleArray(values.size - d - 1)
        System.arraycopy(values, d, shifted, 0, values.size - d - 1)
        Arrays.fill(shifted, values.size - d - 1, shifted.size, 0.0)
        var correlation = if (original.size >= 2)
            PearsonsCorrelation().correlation(original, shifted)
        else
            0.0
        // Safeguard against NaN from Pearson correlation for zero or small vectors,
        // See example at: https://github.com/JetBrains-Research/span/issues/34
        if (correlation.isNaN() || !correlation.isFinite()) {
            correlation = 0.0
        }
        return correlation
    }

}