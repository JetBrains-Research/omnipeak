package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.omnipeak.InputFormat
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FRAGMENTATION_THRESHOLD_BP
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_GAP
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_SENSITIVITY
import org.jetbrains.bio.omnipeak.peaks.OmnipeakModelToPeaks
import org.jetbrains.bio.omnipeak.peaks.Peak
import org.jetbrains.bio.omnipeak.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Omnipeak `compare` invocation.
 *
 * The treatment-control pairs are split into two sets that are to compare.
 *
 * For each treatment-control pair, we compute binned normalized coverage.
 * These coverages are used as the input for a five-state multidimensional negative binomial HMM.
 * For each dimension `d`, there are two negative binomial distributions, low_d and high_d.
 * - ZERO state corresponds to zero emissions for all dimensions
 * - LOW state employs `low_d` emission for each dimension `d`
 * - HIGH state employs `high_d` emission for each dimension `d`
 * - INCREASED state employs `low_d` emission for each dimension `d` from the first set and `high_d` for the second set
 * - DECREASED state employs `high_d` emission for each dimension `d` from the first set and `low_d` for the second set
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class OmnipeakDifferentialPeakCallingExperiment private constructor(
    fitInformation: OmnipeakCompareFitInformation,
    threshold: Double,
    maxIterations: Int
) : OmnipeakModelFitExperiment<ConstrainedNBZHMM, OmnipeakCompareFitInformation, ZLHID>(
    fitInformation,
    ConstrainedNBZHMM.fitter(fitInformation.data1.size, fitInformation.data2.size), ConstrainedNBZHMM::class.java,
    ZLHID.entries.toTypedArray(), NullHypothesis.of(ZLHID.same()),
    threshold, maxIterations,
    null, false, false
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.peak"

    fun computeDirectedDifferencePeaks(
        fdr: Double
    ): Pair<List<Peak>, List<Peak>> {
        val peaks = OmnipeakModelToPeaks.getPeaks(
            results,
            genomeQuery,
            fdr, OmnipeakConstants.OMNIPEAK_DEFAULT_MULTIPLE_TEST_CORRECTION,
            OMNIPEAK_DEFAULT_SENSITIVITY,
            OMNIPEAK_DEFAULT_GAP,
            false,
            OMNIPEAK_DEFAULT_FRAGMENTATION_THRESHOLD_BP,
            OMNIPEAK_DEFAULT_CLIP_MAX_SIGNAL,
        ).peaks
        val highLow = peaks.filter { it.value > 1 }
        val lowHigh = peaks.filter { it.value < 1 }
        return highLow to lowHigh
    }


    companion object {
        /**
         * Technical prefix used for generating track labels.
         */
        private const val OMNIPEAK_TRACK1_PREFIX = "track1"
        private const val OMNIPEAK_TRACK2_PREFIX = "track2"

        /**
         * Creates experiment for model-based comparison of binned coverage tracks for given queries.
         *
         * @return experiment [OmnipeakDifferentialPeakCallingExperiment]
         */
        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths1: List<OmnipeakDataPaths>,
            paths2: List<OmnipeakDataPaths>,
            explicitFormat: InputFormat?,
            bin: Int,
            fragment: Fragment,
            unique: Boolean,
            threshold: Double,
            maxIterations: Int
        ): OmnipeakDifferentialPeakCallingExperiment {
            require(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            val fitInformation = OmnipeakCompareFitInformation.effective(
                genomeQuery,
                paths1, paths2,
                MultiLabels.generate(OMNIPEAK_TRACK1_PREFIX, paths1.size).toList(),
                MultiLabels.generate(OMNIPEAK_TRACK2_PREFIX, paths2.size).toList(),
                explicitFormat, fragment, unique, bin
            )
            return OmnipeakDifferentialPeakCallingExperiment(fitInformation, threshold, maxIterations)
        }
    }
}