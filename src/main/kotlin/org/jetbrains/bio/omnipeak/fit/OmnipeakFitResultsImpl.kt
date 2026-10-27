package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.statistics.model.ClassificationModel
import java.nio.file.Path

class OmnipeakFitResultsImpl(
    fitInfo: OmnipeakFitInformation,
    model: ClassificationModel,
    logNullMemberships: Map<String, DataFrame>,
    statesDataFrameMap: Map<String, DataFrame>?
) : OmnipeakFitResults(fitInfo, model, logNullMemberships, statesDataFrameMap) {

    override fun modelInformation(modelPath: Path): List<TrackAboutMetricValue<*>> {
        val model = this.model
        return when (model) {
            is NB2ZHMM -> {
                val outOfSnrHitDown = model.outOfSignalToNoiseRatioRangeDown
                val outOfLowerNoise = model.outOfLowerNoise
                val statesSwitched = model.statesSwitched
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                listOf(
                    CT_MODEL_FILE to modelPath,
                    CT_MODEL_TYPE to "Negative binomial HMM 2 states with zero inflation",
                    CT_SIGNAL_MEAN to signalMean,
                    CT_NOISE_MEAN to noiseMean,
                    CT_SIGNAL_TO_NOISE to ((signalMean + 1e-10) / (noiseMean + 1e-10)),
                    CT_OUT_OF_SNR_DOWN to outOfSnrHitDown,
                    CT_STATES_SWITCHED to statesSwitched,
                    CT_OUT_OF_NOISE_DOWN to outOfLowerNoise
                )
            }
            else -> super.modelInformation(modelPath)
        }
    }
}
