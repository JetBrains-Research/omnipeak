package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.util.extension
import org.slf4j.LoggerFactory
import java.nio.file.Path

enum class InputFormat {
    BAM(ReadsFormat.BAM),
    SAM(ReadsFormat.SAM),
    CRAM(ReadsFormat.CRAM),
    BED(ReadsFormat.BED),
    BIGWIG(null);

    val internalReadsFormat: ReadsFormat?

    constructor(internalReadsFormat: ReadsFormat?) {
        this.internalReadsFormat = internalReadsFormat
    }

    fun check(path: Path) {
        when (this) {
            BIGWIG -> {
                if (path.extension !in listOf("bw", "bigwig")) {
                    LOG.warn("Unexpected file extension for BIGWIG file: $path, expected bw or bigwig")
                }
            }

            else -> internalReadsFormat!!.check(path)
        }
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(ReadsFormat::class.java)


        fun guess(path: Path): InputFormat? =
            if (path.extension in listOf("bw", "bigwig")) {
                BIGWIG
            }
            else {
                val guess = ReadsFormat.guess(path)
                if (guess != null) InputFormat.valueOf(guess.name) else null
            }
    }

}
