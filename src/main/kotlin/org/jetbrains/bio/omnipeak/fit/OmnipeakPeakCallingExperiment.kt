package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants.OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.omnipeak.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM.Companion.MODEL_EXT
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Omnipeak `analyze --type nbhmm` invocation.
 *
 * For each treatment-control pair, we compute binned normalized coverage.
 * These coverages are used as the input for a three-state multidimensional negative binomial HMM.
 * For each dimension `d`, there are two negative binomial distributions, low_d and high_d.
 * - ZERO state corresponds to zero emissions for all dimensions
 * - LOW state employs `low_d` emission for each dimension `d`
 * - HIGH state employs `high_d` emission for each dimension `d`
 *
 * @author Alexey Dievsky
 * @author Oleg Shpynov
 * @since 10/04/15
 */
class OmnipeakPeakCallingExperiment<Model : ClassificationModel> private constructor(
    fitInformation: OmnipeakAnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepCacheFiles: Boolean
) : OmnipeakModelFitExperiment<Model, OmnipeakAnalyzeFitInformation, ZLH>(
    fitInformation, modelFitter, modelClass, ZLH.entries.toTypedArray(), NullHypothesis.of(ZLH.Z, ZLH.L),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepCacheFiles
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.$MODEL_EXT"

    companion object {

        /**
         * Technical prefix used for generating track labels.
         */
        const val OMNIPEAK_TRACK_PREFIX = "track"

        /**
         * Creates experiment for model-based enrichment of binned coverage tracks (e.g. ChIP-seq tracks)
         * for given number of [paths].
         * Not restricted for single query and constrained for multiple paths.
         *
         * @return experiment [OmnipeakPeakCallingExperiment]
         */
        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths: List<OmnipeakDataPaths>,
            explicitFormat: ReadsFormat?,
            fragment: Fragment = AutoFragment,
            unique: Boolean = true,
            bin: Int,
            hmmEstimateSNR: Double = OMNIPEAK_DEFAULT_HMM_ESTIMATE_SNR,
            hmmLow: Double = OMNIPEAK_DEFAULT_HMM_LOW_THRESHOLD,
            fixedModelPath: Path? = null,
            threshold: Double = OMNIPEAK_DEFAULT_FIT_THRESHOLD,
            maxIterations: Int = OMNIPEAK_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo: Boolean = false,
            keepCacheFiles: Boolean = false
        ): OmnipeakPeakCallingExperiment<out ClassificationModel> {
            require(paths.isNotEmpty()) { "No data" }
            val fitInformation = OmnipeakAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, explicitFormat, MultiLabels.generate(OMNIPEAK_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            return if (paths.size == 1) {
                OmnipeakPeakCallingExperiment(
                    fitInformation,
                    NB2ZHMM.fitter(hmmEstimateSNR, hmmLow),
                    NB2ZHMM::class.java,
                    fixedModelPath,
                    threshold, maxIterations,
                    saveExtendedInfo,
                    keepCacheFiles
                )
            } else {
                OmnipeakPeakCallingExperiment(
                    fitInformation,
                    ConstrainedNBZHMM.fitter(paths.size),
                    ConstrainedNBZHMM::class.java,
                    fixedModelPath,
                    threshold, maxIterations,
                    saveExtendedInfo,
                    keepCacheFiles
                )
            }
        }
    }
}

