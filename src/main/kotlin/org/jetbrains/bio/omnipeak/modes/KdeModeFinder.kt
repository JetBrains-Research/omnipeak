package org.jetbrains.bio.omnipeak.modes

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.omnipeak.fit.OmnipeakConstants
import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.exp

/** KDE smoothing (Gaussian) + local-maximum mode detection for 1-D signals.  */
object KdeModeFinder {

    fun detectModes(
        signal: DoubleArray,
        bandwidth: Int,
        start: Int = 0,
        end: Int = signal.size,
        workArray: DoubleArray? = null,
        minLx: Double = OmnipeakConstants.OMNIPEAK_SUMMITS_MIN_LENGTH,
        minDx: Double = OmnipeakConstants.OMNIPEAK_SUMMITS_MIN_DISTANCE,
        minDy: Double = 0.1
    ): List<Range> {
        if (end <= start) return emptyList()

        val size = end - start

        // 1) Build Gaussian kernel
        val (rad, kernel) = buildKernel(bandwidth)

        // 2) Apply Gaussian smoothing
        val smooth = applyKernel(rad, kernel, signal, start, end, workArray)

        // 3) Find sample-local maxima
        val localMaxs = BitList(size)
        (1 until size - 1).filter {
            smooth[it - 1] <= smooth[it] && smooth[it] >= smooth[it + 1]
        }.forEach { localMaxs.set(it) }
        val seeds = localMaxs.aggregate()
        if (seeds.isEmpty()) return emptyList()

        // 4) Fetch modes in smoothed signal
        val modes = ArrayList<Range>()
        for ((i, s) in seeds.withIndex()) {
            val leftBoundary = if (modes.isEmpty()) 0 else modes.last().endOffset + 1
            val rightBoundary = if (i < seeds.size - 1) seeds[i + 1].startOffset - rad + 1 else size - 1

            var left = s.startOffset
            var leftSlopeFound = false
            while (left > leftBoundary && smooth[left - 1] <= smooth[left]) {
                if (abs((smooth[left - 1] - smooth[left]) / smooth[left]) > minDy) {
                    leftSlopeFound = true
                } else {
                    if (leftSlopeFound) break
                }
                left--
            }

            var right = s.endOffset
            var rightSlopeFound = false
            while (right < rightBoundary && smooth[right + 1] <= smooth[right]) {
                if (abs(smooth[right + 1] - smooth[right] / smooth[right]) > minDy) {
                    rightSlopeFound = true
                } else {
                    if (rightSlopeFound) break
                }
                right++
            }

            if (left < right) {
                modes.add(Range(left, right))
            }
        }
        // Ensure enough space between modes and min length
        val minL = minLx * bandwidth
        val minD = minDx * bandwidth
        val result = arrayListOf<Range>()
        for (i in 0 until modes.size) {
            val current = modes[i]
            if (current.endOffset - current.startOffset < minL) {
                continue
            }
            if (result.isEmpty()) {
                result.add(current)
                continue
            }
            val prevMode = result.last()
            val d = current.startOffset - prevMode.endOffset
            if (d >= minD) {
                result.add(current)
                continue
            }
            val de = ceil((minD - d) * 0.5).toInt()
            if (prevMode.endOffset - prevMode.startOffset - de >= minL &&
                current.endOffset - current.startOffset - de >= minL
            ) {
                result[result.lastIndex] = Range(prevMode.startOffset, prevMode.endOffset - de)
                result.add(Range(current.startOffset + de, current.endOffset))
            } else {
                result[result.lastIndex] = Range(prevMode.startOffset, current.endOffset)
            }
        }
        return result
    }

    fun buildKernel(
        bandwidth: Int,
    ): Pair<Int, DoubleArray> {
        val sigma = bandwidth * 0.5
        val rad = ceil(bandwidth * 0.5).toInt()
        val kernel = DoubleArray(2 * rad + 1)
        var sum = 0.0
        val var2 = 2.0 * sigma * sigma
        for (k in -rad..rad) {
            val w = exp(-(k * k) / var2)
            kernel[k + rad] = w
            sum += w
        }
        for (i in kernel.indices) kernel[i] /= sum // normalize to area 1
        return rad to kernel
    }

    fun applyKernel(
        rad: Int,
        kernel: DoubleArray,
        signal: DoubleArray,
        start: Int = 0,
        end: Int = signal.size,
        workArray: DoubleArray? = null
    ): DoubleArray {
        val smooth = if (workArray != null && workArray.size >= end - start) {
            workArray.fill(0.0)
            workArray
        } else
            DoubleArray(end - start)
        for (i in 0 until end - start) {
            var acc = 0.0
            for (k in -rad..rad) {
                val j = start + i + k
                if (j in start until end) {
                    acc += kernel[k + rad] * signal[j]
                }
            }
            smooth[i] = acc
        }
        return smooth
    }
}
