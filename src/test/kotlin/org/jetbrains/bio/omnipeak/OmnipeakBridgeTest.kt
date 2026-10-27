package org.jetbrains.bio.omnipeak

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.omnipeak.fit.OmnipeakAnalyzeFitInformation
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
import kotlin.test.assertTrue

class OmnipeakBridgeTest {

    @Test
    fun testBridgeInstantiationViaReflection() {
        val clazz = Class.forName(OmnipeakBridge.BRIDGE_CLASS_NAME)
        val instance = clazz.getDeclaredConstructor().newInstance()
        assertTrue(instance is OmnipeakBridge, "Instance must implement OmnipeakBridge")
    }

    @Test
    fun testBridgeMethods() {
        val bridge: OmnipeakBridge = OmnipeakBridgeImpl()
        assertTrue(bridge.recommendedMemoryGb() > 0)

        // Gap calculation through bridge
        assertEquals(0, bridge.estimateGapFromCandidateNs(intArrayOf()))
        assertEquals(0, bridge.estimateGapFromCandidateNs(intArrayOf(100)))

        val genomeQuery = GenomeQuery(Genome["to1"], "chr1", "chr2")
        val fitInfo = OmnipeakAnalyzeFitInformation(
            genomeQuery, emptyList(), emptyList(),
            null, AutoFragment, true, 200, true
        )
        val sensitivity = -2.0
        val logNullMembershipsMap = genomeMap(genomeQuery) { chr ->
            val len = (chr.length / 200) + 1
            F64Array(len) { (it * 37 % 100).toDouble() / -20.0 }
        }
        val bitList2reuseMap = genomeMap(genomeQuery) { chr ->
            BitList(logNullMembershipsMap[chr].length)
        }

        val candidateGapNs = bridge.estimateCandidateGapNs(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
            sensitivity, maxGap = 10, parallel = false
        )
        assertEquals(11, candidateGapNs.size)

        val (gap, candidateNs) = bridge.estimateGap(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
            sensitivity, maxGap = 10, parallel = false
        )
        assertEquals(candidateGapNs.size, candidateNs.size)
        assertTrue(gap >= 0)
    }
}
