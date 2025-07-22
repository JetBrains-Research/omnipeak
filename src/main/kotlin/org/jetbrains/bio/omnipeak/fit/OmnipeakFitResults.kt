package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.TrackAboutBooleanColumnType
import org.jetbrains.bio.genome.TrackAboutDoubleColumnType
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.genome.TrackAboutStringColumnType
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * Contains the results of a Span-like model-fitting experiment.
 *
 * @property fitInfo The [OmnipeakFitInformation] instance that describes the experiment input.
 * @property model The [ClassificationModel] that was fitted during the experiment.
 * @property logNullMemberships The chromosome-wise dataframes of log null probabilities, i.e.
 * the log probability of each observation under the null hypothesis. Each dataframe should at least contain
 * a column of floats or doubles labelled [OmnipeakModelFitExperiment.NULL].
 * @property statesDataFrameMap extended information which contains additional states per position mapping
 */
open class OmnipeakFitResults(
    val fitInfo: OmnipeakFitInformation,
    val model: ClassificationModel,
    val logNullMemberships: Map<String, DataFrame>,
    val statesDataFrameMap: Map<String, DataFrame>?
) {
    companion object {
        internal val LOG = LoggerFactory.getLogger(OmnipeakFitResults::class.java)

        val CT_MODEL_FILE = TrackAboutStringColumnType("Model")
        val CT_MODEL_TYPE = TrackAboutStringColumnType("Model type")
        val CT_SIGNAL_MEAN = TrackAboutDoubleColumnType("Signal mean")
        val CT_NOISE_MEAN = TrackAboutDoubleColumnType("Noise mean")
        val CT_SIGNAL_TO_NOISE = TrackAboutDoubleColumnType("Signal to noise")
        val CT_OUT_OF_SNR_DOWN = TrackAboutBooleanColumnType("Out of signal-to-noise range down")
        val CT_OUT_OF_NOISE_DOWN = TrackAboutBooleanColumnType("Out of low noise level down")
        val CT_STATES_SWITCHED = TrackAboutBooleanColumnType("States switched")
    }

    /**
     * @return Information about fit results including model and other parameters
     */
    open fun modelInformation(modelPath: Path): List<TrackAboutMetricValue<*>> {
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
            else -> emptyList()
        }
    }
}
