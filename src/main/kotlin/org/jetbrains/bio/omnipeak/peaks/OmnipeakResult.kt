package org.jetbrains.bio.omnipeak.peaks

data class OmnipeakResult(
    val fdr: Double,
    val sensitivity: Double,
    val gap: Int,
    val peaks: List<Peak>
)