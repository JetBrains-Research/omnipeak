package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.omnipeak.fit.OmnipeakCompareFitInformation
import org.jetbrains.bio.omnipeak.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.omnipeak.statistics.hmm.FreeNBZHMM
import org.jetbrains.bio.omnipeak.statistics.hmm.NB2ZHMM

object SPAN2 {

    // Keep backward compatibility with SPAN and SPAN2
    fun toOmnipeak(text: String): String = text
        .replace(
            "org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation",
            OmnipeakAnalyzeFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.span.fit.AnalyzeFitInformation",
            OmnipeakAnalyzeFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.span.AnalyzeFitInformation",
            OmnipeakAnalyzeFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.experiments.fit.SpanAnalyzeFitInformation",
            OmnipeakAnalyzeFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.experiments.fit.AnalyzeFitInformation",
            OmnipeakAnalyzeFitInformation::class.qualifiedName!!
        ).replace(
            "\"AnalyzeFitInformation\"",
            "\"${OmnipeakAnalyzeFitInformation::class.qualifiedName!!}\""
        ).replace(
            "\"SpanAnalyzeFitInformation\"",
            "\"${OmnipeakAnalyzeFitInformation::class.qualifiedName!!}\""
        ).replace(
            "org.jetbrains.bio.span.fit.SpanCompareFitInformation",
            OmnipeakCompareFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.span.fit.CompareFitInformation",
            OmnipeakCompareFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.span.CompareFitInformation",
            OmnipeakCompareFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.experiments.fit.SpanCompareFitInformation",
            OmnipeakCompareFitInformation::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.experiments.fit.CompareFitInformation",
            OmnipeakCompareFitInformation::class.qualifiedName!!
        ).replace(
            "\"CompareFitInformation\"",
            "\"${OmnipeakCompareFitInformation::class.qualifiedName!!}\""
        ).replace(
            "\"SpanCompareFitInformation\"",
            "\"${OmnipeakCompareFitInformation::class.qualifiedName!!}\""
        ).replace(
            "org.jetbrains.bio.span.statistics.hmm.NB2ZHMM",
            NB2ZHMM::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.span.statistics.hmm.ConstrainedNBZHMM",
            ConstrainedNBZHMM::class.qualifiedName!!
        ).replace(
            "org.jetbrains.bio.span.statistics.hmm.FreeNBZHMM",
            FreeNBZHMM::class.qualifiedName!!
        )

}