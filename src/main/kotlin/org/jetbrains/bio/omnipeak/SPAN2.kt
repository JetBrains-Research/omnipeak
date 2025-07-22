package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM

object SPAN2 {

    // Keep backward compatibility with SPAN2
    fun toOmnipeak(text: String): String = text
        .replace(
            "org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation",
            OmnipeakAnalyzeFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM",
            NB2ZHMM::class.qualifiedName!!
        )

}