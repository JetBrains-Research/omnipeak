package org.jetbrains.bio.omnipeak.peaks

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_AUTOCORRELATION_MAX_SHIFT
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks.getLogNullPvals
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

    fun avgAutocorrelation(
        fitResults: OmnipeakFitResults,
        fitInfo: OmnipeakFitInformation = fitResults.fitInfo,
        genomeQuery: GenomeQuery = fitInfo.genomeQuery(),
        blackList: LocationsMergingList? = null
    ): Double {
        val logNullPvals = getLogNullPvals(genomeQuery, fitResults, blackList)
        if (logNullPvals.size < 2) {
            return 0.0
        }
        val correlations = computeAutoCorrelations(logNullPvals)
        val avg = correlations.average()
        return if (avg.isFinite()) avg.coerceIn(0.0, 1.0) else 0.0
    }

    private fun computeAutocorrelation(values: DoubleArray, d: Int): Double {
        val len = values.size - d
        if (len < 2) return 0.0
        val original = DoubleArray(len)
        System.arraycopy(values, 0, original, 0, len)
        val shifted = DoubleArray(len)
        System.arraycopy(values, d, shifted, 0, len)
        var correlation = PearsonsCorrelation().correlation(original, shifted)
        // Safeguard against NaN from Pearson correlation for zero or small vectors,
        // See example at: https://github.com/JetBrains-Research/span/issues/34
        if (correlation.isNaN() || !correlation.isFinite()) {
            correlation = 0.0
        }
        return correlation
    }

}