package org.jetbrains.bio.omnipeak.modes

import java.util.*
import kotlin.math.max
import kotlin.math.min


/**
 * Hartigan & Hartigan (1985) "The Dip Test of Unimodality"
 * See https://www.jstor.org/stable/2241144
 *
 * @author Oleg Shpynov
 * @since 10/07/2015
*/
object HartiganDip {

    data class DipResult(
        val dip: Double, val pValue: Double
    )

    private val R = Random(1234567L)

    /** Convenience: dip only (no p-value).  */
    fun dipStatistic(x: DoubleArray): Double {
        val n = x.size
        if (n < 3) return 0.0

        // Sort array
        Arrays.sort(x)

        // Empirical CDF at order stats is i/n for i=1..n. We'll work with indices 0..n-1 and use (i+1)/n.
        // Build GCM (greatest convex minorant) of the empirical CDF graph (x_i, i/n).
        val gKnots = buildGCM(x) // increasing slopes
        val lKnots = buildLCM(x) // decreasing slopes

        // Evaluate the vertical gap between LCM and GCM across the support and take half of its maximum.
        // Because each is piecewise linear in x, the max difference occurs at a knot or data point.
        val maxGap = maxVerticalGap(x, gKnots, lKnots)
        return 0.5 * maxGap
    }

    /** Dip + bootstrap p-value with simulations under U[0,1].  */
    fun dipTest(x: DoubleArray, bootstraps: Int): DipResult {
        val dip = dipStatistic(x)
        if (dip.isNaN() || bootstraps <= 0) return DipResult(dip, Double.NaN)

        // Null: i.i.d. unimodal with least-favorable = Uniform[0,1]
        // (Hartigan & Hartigan show uniform yields largest dips under the null).
        if (x.size < 3) return DipResult(0.0, 1.0)

        var ge = 0
        val sim = DoubleArray(x.size)
        repeat(bootstraps) {
            for (i in x.indices) sim[i] = R.nextDouble()
            val d0 = dipStatistic(sim)
            if (d0 >= dip) ge++
        }
        val p = (ge + 1.0) / (bootstraps + 1.0) // small-sample smoothing
        return DipResult(dip, p)
    }

    // Build greatest convex minorant of points (x[i], y[i]=(i+1)/n).
    // Return indices of knot points in ascending order.
    private fun buildGCM(x: DoubleArray): IntArray {
        val n = x.size
        val stack = IntArray(n)
        var top = -1
        stack[++top] = 0
        stack[++top] = 1
        while (stack[top] < n - 1) {
            stack[++top] = stack[top - 1] + 1 // propose next point
            // enforce nondecreasing slopes
            while (top >= 2 && slope(x, stack[top - 2], stack[top - 1], n) >=
                slope(x, stack[top - 1], stack[top], n)
            ) {
                // pool adjacent segments
                stack[top - 1] = stack[top]
                top--
            }
        }
        return stack.copyOf(top + 1)
    }

    // Build least concave majorant of points (x[i], y[i]=(i+1)/n).
    // Return indices of knot points in ascending order.
    private fun buildLCM(x: DoubleArray): IntArray {
        val n = x.size
        val stack = IntArray(n)
        var top = -1
        stack[++top] = n - 1
        stack[++top] = n - 2
        while (stack[top] > 0) {
            stack[++top] = stack[top - 1] - 1 // propose previous point
            // enforce nonincreasing slopes
            while (top >= 2 && slope(x, stack[top - 2], stack[top - 1], n) <=
                slope(x, stack[top - 1], stack[top], n)
            ) {
                // pool adjacent segments
                stack[top - 1] = stack[top]
                top--
            }
        }
        // reverse to ascending
        val asc = IntArray(top + 1)
        for (i in 0..top) asc[i] = stack[top - i]
        return asc
    }

    // Slope of the secant of the empirical CDF between i and j (i<j): ( (j+1)/n - (i+1)/n ) / (x[j]-x[i])
    private fun slope(x: DoubleArray, i: Int, j: Int, n: Int): Double {
        if (x[j] == x[i]) return Double.POSITIVE_INFINITY // vertical in x -> infinite slope
        return ((j - i) * 1.0 / n) / (x[j] - x[i])
    }

    // Evaluate piecewise-linear GCM and LCM at each x[i], take maximum vertical gap.
    private fun maxVerticalGap(x: DoubleArray, gKnots: IntArray, lKnots: IntArray): Double {
        val n = x.size
        // Precompute envelopes y_gcm[i], y_lcm[i]
        val g = DoubleArray(n)
        val l = DoubleArray(n)

        // Fill GCM between knots with linear interpolation of y=(i+1)/n
        for (s in 0 until gKnots.size - 1) {
            val i0 = gKnots[s]
            val i1 = gKnots[s + 1]
            val x0 = x[i0]
            val x1 = x[i1]
            val y0 = (i0 + 1) / n.toDouble()
            val y1 = (i1 + 1) / n.toDouble()
            if (x1 == x0) {
                for (i in i0..i1) g[i] = max(g[i], max(y0, y1))
            } else {
                val a = (y1 - y0) / (x1 - x0)
                val b = y0 - a * x0
                for (i in i0..i1) g[i] = a * x[i] + b
            }
        }

        // Fill LCM between knots with linear interpolation
        for (s in 0 until lKnots.size - 1) {
            val i0 = lKnots[s]
            val i1 = lKnots[s + 1]
            val x0 = x[i0]
            val x1 = x[i1]
            val y0 = (i0 + 1) / n.toDouble()
            val y1 = (i1 + 1) / n.toDouble()
            if (x1 == x0) {
                for (i in i0..i1) l[i] =
                    min(if (l[i] == 0.0) Double.POSITIVE_INFINITY else l[i], min(y0, y1))
            } else {
                val a = (y1 - y0) / (x1 - x0)
                val b = y0 - a * x0
                for (i in i0..i1) l[i] = a * x[i] + b
            }
        }

        // The empirical CDF is stepwise i/n; the unimodal “pinch” lies between GCM and LCM.
        // Dip is half of the sup gap between them.
        var maxGap = 0.0
        for (i in 0 until n) {
            // ensure envelopes sandwich the empirical CDF: numerical guards
            val gi = min(g[i], (i + 1) / n.toDouble())
            val li = max(l[i], (i + 1) / n.toDouble())
            val gap = li - gi
            if (gap > maxGap) maxGap = gap
        }
        return maxGap
    }
}
