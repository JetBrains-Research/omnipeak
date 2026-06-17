package org.jetbrains.bio.omnipeak.peaks

import org.jetbrains.bio.omnipeak.fit.OmnipeakFitResults
import org.jetbrains.bio.omnipeak.fit.ZLH
import org.jetbrains.bio.omnipeak.fit.ZLHID
import org.jetbrains.bio.omnipeak.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.hmm.MLAbstractHMM
import org.jetbrains.bio.viktor.F64Array
import org.knowm.xchart.*
import org.knowm.xchart.internal.chartpart.Chart
import org.knowm.xchart.style.Styler
import org.knowm.xchart.style.markers.SeriesMarkers
import org.knowm.xchart.style.theme.GGPlot2Theme
import org.slf4j.LoggerFactory
import java.awt.Color
import java.nio.file.Path
import kotlin.math.exp
import kotlin.math.max

object OmnipeakModelPlots {
    private val LOG = LoggerFactory.getLogger(OmnipeakModelPlots::class.java)

    // ColorBrewer "Blues" gradient (light -> dark) for heatmaps, min -> max
    private val BLUE_PALETTE = arrayOf(
        Color(0xF7FBFF),
        Color(0xC6DBEF),
        Color(0x6BAED6),
        Color(0x2171B5),
        Color(0x08306B)
    )

    fun prepareModelPlots(
        peaksPath: Path?,
        fitResults: OmnipeakFitResults
    ) {
        if (peaksPath == null) return
        val model = fitResults.model
        if (model is MLAbstractHMM) {
            LOG.info("Plotting model plots...")
            prepareModelDistributionsPlot(peaksPath, model)
            prepareInitialMatrixPlot(peaksPath, model)
            prepareTransitionMatrixPlot(peaksPath, model)
        }
    }

    private fun getStates(model: MLAbstractHMM): List<String> {
        val numStates = model.logPriorProbabilities.length
        return when (model) {
            is NB2ZHMM -> ZLH.entries.map { it.name }
            is ConstrainedNBZHMM -> {
                if (numStates == 3) ZLH.entries.map { it.name }
                else if (numStates == 5) ZLHID.entries.map { it.name }
                else (0 until numStates).map { "S$it" }
            }
            else -> (0 until numStates).map { "S$it" }
        }
    }

    private fun prepareModelDistributionsPlot(peaksPath: Path, model: MLAbstractHMM) {
        val means: F64Array
        val failures: F64Array
        val numStates = model.logPriorProbabilities.length
        val states = getStates(model)
        val negBinStates: List<String>
        
        when (model) {
            is NB2ZHMM -> {
                means = model.means
                failures = model.failures
                negBinStates = listOf("Low", "High")
            }
            is ConstrainedNBZHMM -> {
                means = model.means
                failures = model.failures
                negBinStates = if (numStates == 3) {
                    listOf("Low", "High")
                } else if (numStates == 5) {
                    listOf("Low", "Increased", "Decreased", "High")
                } else {
                    (0 until means.length).map { "NB$it" }
                }
            }
            else -> return
        }

        // Determine threshold for plotting (similar to PeakCallerModelRenderer)
        val maxMean = means.toDoubleArray().maxOrNull() ?: 10.0
        val threshold = max(10, (maxMean * 3).toInt())

        val chart = XYChartBuilder()
            .width(800)
            .height(600)
            .title("Model States Distributions")
            .xAxisTitle("Counts")
            .yAxisTitle("Probability")
            .build()

        chart.styler.theme = GGPlot2Theme()
        chart.styler.plotBackgroundColor = Color.WHITE
        chart.styler.chartBackgroundColor = Color.WHITE
        chart.styler.plotGridLinesColor = Color.LIGHT_GRAY
        chart.styler.legendPosition = Styler.LegendPosition.InsideNE
        chart.styler.defaultSeriesRenderStyle = XYSeries.XYSeriesRenderStyle.Line

        for (i in 0 until means.length) {
            val mean = means[i]
            val failure = failures[i]
            val nb = NegativeBinomialDistribution.usingMean(mean, failure)
            val xData = DoubleArray(threshold + 1) { it.toDouble() }
            val yData = DoubleArray(threshold + 1) {
                val p = exp(nb.logProbability(it))
                if (p.isNaN() || !p.isFinite()) 0.0 else p
            }
            val stateName = if (i < negBinStates.size) negBinStates[i] else "NB$i"
            val seriesName = "$stateName (Mean=${"%.2f".format(mean)}, Failures=${"%.2f".format(failure)})"
            val series = chart.addSeries(seriesName, xData, yData)
            series.marker = SeriesMarkers.NONE
        }

        savePlot(peaksPath, "model_distributions", chart)
    }

    private fun prepareInitialMatrixPlot(peaksPath: Path, model: MLAbstractHMM) {
        val states = getStates(model)
        val priors = model.priorProbabilities
        
        val chart = HeatMapChartBuilder()
            .width(800)
            .height(600)
            .title("Initial Matrix (Priors)")
            .xAxisTitle("")
            .yAxisTitle("State")
            .build()

        chart.styler.theme = GGPlot2Theme()
        chart.styler.plotBackgroundColor = Color.WHITE
        chart.styler.chartBackgroundColor = Color.WHITE
        chart.styler.plotGridLinesColor = Color.LIGHT_GRAY
        chart.styler.isLegendVisible = false
        chart.styler.setShowValue(true)
        chart.styler.heatMapValueDecimalPattern = "0.000"
        chart.styler.rangeColors = BLUE_PALETTE

        // XChart heat data is a list of [xIndex, yIndex, value] triplets
        val xData = listOf("Prior")
        val yData = states
        val zData = states.indices.map { i ->
            arrayOf<Number>(0, i, priors[i])
        }

        chart.addSeries("Priors", xData, yData, zData)
        
        savePlot(peaksPath, "model_priors", chart)
    }

    private fun prepareTransitionMatrixPlot(peaksPath: Path, model: MLAbstractHMM) {
        val states = getStates(model)
        val transitions = model.transitionProbabilities

        val chart = HeatMapChartBuilder()
            .width(800)
            .height(600)
            .title("Transition Matrix")
            .xAxisTitle("To State")
            .yAxisTitle("From State")
            .build()

        chart.styler.theme = GGPlot2Theme()
        chart.styler.plotBackgroundColor = Color.WHITE
        chart.styler.chartBackgroundColor = Color.WHITE
        chart.styler.isLegendVisible = false
        chart.styler.setShowValue(true)
        chart.styler.heatMapValueDecimalPattern = "0.000"
        chart.styler.rangeColors = BLUE_PALETTE

        // XChart heat data is a list of [xIndex, yIndex, value] triplets.
        // x = "To State" (column j), y = "From State" (row i)
        val xData = states
        val yData = states
        val zData = states.indices.flatMap { i ->
            states.indices.map { j ->
                arrayOf<Number>(j, i, transitions[i, j])
            }
        }

        chart.addSeries("Transitions", xData, yData, zData)

        savePlot(peaksPath, "model_transitions", chart)
    }

    private fun savePlot(peaksPath: Path, name: String, chart: Chart<*, *>) {
        val plotFile = "$peaksPath.$name.png"
        try {
            BitmapEncoder.saveBitmapWithDPI(chart, plotFile, BitmapEncoder.BitmapFormat.PNG, 300)
            LOG.info("See $plotFile")
        } catch (e: Exception) {
            LOG.error("Failed to save plot $plotFile", e)
        }
    }
}
