package org.jetbrains.bio.omnipeak.peaks

import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.util.bufferedWriter
import java.nio.file.Path

fun Peak.Companion.savePeaks(peaks: List<Peak>, path: Path, peakName: String = "peak") {
    CSVFormat.TDF.print(path.bufferedWriter()).use { printer ->
        peaks.sorted().forEachIndexed { i, peak ->
            /* See ENCODE broadPeak (BED 6+3) format for details */
            printer.printRecord(
                peak.chromosome.name,
                peak.range.startOffset.toString(),
                peak.range.endOffset.toString(),
                "${peakName}_${i + 1}",
                peak.score.toString(),
                ".",
                peak.value.toString(),
                peak.mlogpvalue.toString(),
                peak.mlogqvalue.toString()
            )
        }
    }
}
