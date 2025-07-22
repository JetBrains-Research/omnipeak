package org.jetbrains.bio.omnipeak.fit

import java.nio.file.Path

/**
 * Data class to represent treatment and control pair.
 * Later it is used to produce normalized coverage.
 */
data class OmnipeakDataPaths(val treatment: Path, val control: Path?) {
    override fun toString(): String {
        return "experiment ${treatment.toAbsolutePath()}" +
                if (control != null) " with control ${control.toAbsolutePath()}" else ""
    }
}

