package org.jetbrains.bio.omnipeak.peaks

@Suppress("ArrayInDataClass")
data class OmnipeakResult(
    val fdr: Double,
    val sensitivities: DoubleArray,
    val gap: Int,
    val peaks: List<Peak>
) {
    val sensitivity: Double get() = sensitivities.first()

    constructor(
        fdr: Double,
        sensitivity: Double,
        gap: Int,
        peaks: List<Peak>
    ) : this(fdr, doubleArrayOf(sensitivity), gap, peaks)
}
