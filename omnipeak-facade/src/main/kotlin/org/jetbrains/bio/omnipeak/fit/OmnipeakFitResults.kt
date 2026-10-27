package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.TrackAboutBooleanColumnType
import org.jetbrains.bio.genome.TrackAboutDoubleColumnType
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.genome.TrackAboutStringColumnType
import org.jetbrains.bio.statistics.model.ClassificationModel
import java.nio.file.Path

open class OmnipeakFitResults(
    val fitInfo: OmnipeakFitInformation,
    val model: ClassificationModel,
    val logNullMemberships: Map<String, DataFrame>,
    val statesDataFrameMap: Map<String, DataFrame>?
) {
    open fun modelInformation(modelPath: Path): List<TrackAboutMetricValue<*>> = emptyList()

    companion object {
        const val NULL = "null"

        val CT_MODEL_FILE = TrackAboutStringColumnType("Model")
        val CT_MODEL_TYPE = TrackAboutStringColumnType("Model type")
        val CT_SIGNAL_MEAN = TrackAboutDoubleColumnType("Signal mean")
        val CT_NOISE_MEAN = TrackAboutDoubleColumnType("Noise mean")
        val CT_SIGNAL_TO_NOISE = TrackAboutDoubleColumnType("Signal to noise")
        val CT_OUT_OF_SNR_DOWN = TrackAboutBooleanColumnType("Out of signal-to-noise range down")
        val CT_OUT_OF_NOISE_DOWN = TrackAboutBooleanColumnType("Out of low noise level down")
        val CT_STATES_SWITCHED = TrackAboutBooleanColumnType("States switched")
    }
}
