package org.jetbrains.bio.omnipeak.fit

import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat

interface AbstractOmnipeakAnalyzeFitInformation : OmnipeakFitInformation {
    val paths: List<OmnipeakDataPaths>
    val explicitFormat: ReadsFormat?
    val fragment: Fragment
    val unique: Boolean
}
