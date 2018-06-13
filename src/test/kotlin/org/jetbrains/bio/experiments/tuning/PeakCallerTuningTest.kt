package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.big.ExtendedBedEntry
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.junit.Test
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Created by Aleksei Dievskii on 16.02.2018.
 */
class PeakCallerTuningTest {

    val RANDOM = Random(1234)

    @Test fun testErrorRateDegenerate() {
        val errorRate = ErrorRate()
        assertTrue(errorRate.rate().isNaN(), "Empty error rate is not NaN.")
    }

    @Test fun testErrorRateBasic() {
        val errorRate1 = ErrorRate()
        arrayOf(true, true, false, true, false).forEach(errorRate1::observe)
        assertEquals(0.6, errorRate1.rate())
        val errorRate2 = ErrorRate()
        arrayOf(true, true, true).forEach(errorRate2::observe)
        assertEquals(1.0, errorRate2.rate())
        errorRate1.combine(errorRate2)
        assertEquals(0.75, errorRate1.rate())
    }

    @Test fun testErrorRateCommutative() {
        val outcomes = (0 until 100).map { RANDOM.nextBoolean() }
        val errorRate = ErrorRate()
        outcomes.forEach(errorRate::observe)
        val outcomesShuffled = outcomes.shuffled(RANDOM)
        val errorRateShuffled = ErrorRate()
        outcomesShuffled.forEach(errorRateShuffled::observe)
        assertEquals(errorRate.rate(), errorRateShuffled.rate(), "The error rate is not commutative.")
    }

    @Test fun testErrorRateAssociative() {
        val outcomes1 = (0 until 100).map { RANDOM.nextBoolean() }
        val outcomes2 = (0 until 100).map { RANDOM.nextBoolean() }
        val errorRate1 = ErrorRate()
        val errorRate2 = ErrorRate()
        val errorRate3 = ErrorRate() // same as errorRate1, introduced since `coverageDataFrame` is destructive
        outcomes1.forEach(errorRate1::observe)
        outcomes2.forEach(errorRate2::observe)
        outcomes1.forEach(errorRate3::observe)
        errorRate1.combine(errorRate2)
        errorRate2.combine(errorRate3)
        assertEquals(errorRate1.rate(), errorRate2.rate(), "The error rate is not associative.")
    }

    @Test fun testLabelErrorsBasic() {
        val chr1 = Chromosome("to1", "chr1")
        val labelErrors = LabelErrors()
        val label1 = PeakAnnotation(Location(10, 20, chr1), PeakAnnotationType.NO_PEAKS)
        val label2 = PeakAnnotation(Location(30, 40, chr1), PeakAnnotationType.PEAKS)
        arrayOf(true, true, false, true, false).forEach { labelErrors.observe(label1, it) }
        arrayOf(true, true, true).forEach { labelErrors.observe(label2, it) }
        assertEquals(0.4, labelErrors.error(PeakAnnotationType.NO_PEAKS))
        assertEquals(0.0, labelErrors.error(PeakAnnotationType.PEAKS))
        assertEquals(0.25, labelErrors.error())
    }

    @Test fun testLabelErrorsCommutative() {
        val chr1 = Chromosome("to1", "chr1")
        val labelErrors = LabelErrors()
        val labelErrorsShuffled = LabelErrors()
        val label1 = PeakAnnotation(Location(10, 20, chr1), PeakAnnotationType.NO_PEAKS)
        val label2 = PeakAnnotation(Location(30, 40, chr1), PeakAnnotationType.PEAKS)
        val label3 = PeakAnnotation(Location(50, 60, chr1), PeakAnnotationType.NO_PEAKS)
        val labels = arrayListOf(label1, label2, label3)
        val outcomes = (0 until 100).map { RANDOM.nextBoolean() }
        val labelIndices = (0 until 100).map { RANDOM.nextInt(3) }
        (0 until 100).forEach { labelErrors.observe(labels[labelIndices[it]], outcomes[it]) }
        val shuffle = (0 until 100).shuffled(RANDOM)
        shuffle.forEach { labelErrorsShuffled.observe(labels[labelIndices[it]], outcomes[it]) }
        assertEquals(labelErrors.map, labelErrorsShuffled.map, "Label errors are not commutative.")
    }

    @Test fun testLabelErrorsAssociative() {
        val chr1 = Chromosome("to1", "chr1")
        val labelErrors1 = LabelErrors()
        val labelErrors2 = LabelErrors()
        val labelErrors3 = LabelErrors()
        val label1 = PeakAnnotation(Location(10, 20, chr1), PeakAnnotationType.NO_PEAKS)
        val label2 = PeakAnnotation(Location(30, 40, chr1), PeakAnnotationType.PEAKS)
        val label3 = PeakAnnotation(Location(50, 60, chr1), PeakAnnotationType.NO_PEAKS)
        val labels = arrayListOf(label1, label2, label3)
        val outcomes1 = (0 until 100).map { RANDOM.nextBoolean() }
        val outcomes2 = (0 until 100).map { RANDOM.nextBoolean() }
        val labelIndices1 = (0 until 100).map { RANDOM.nextInt(3) }
        val labelIndices2 = (0 until 100).map { RANDOM.nextInt(3) }
        (0 until 100).forEach {
            labelErrors1.observe(labels[labelIndices1[it]], outcomes1[it])
            labelErrors2.observe(labels[labelIndices2[it]], outcomes2[it])
            labelErrors3.observe(labels[labelIndices1[it]], outcomes1[it]) // 'coverageDataFrame' is destructive, we need a copy
        }
        labelErrors1.combine(labelErrors2)
        labelErrors2.combine(labelErrors3)
        assertEquals(labelErrors1.map, labelErrors2.map, "Label errors are not associative.")
    }

    @Test fun testLabelBasic() {
        val chr1 = Chromosome("to1", "chr1")
        val chrX = Chromosome("to1", "chrX")
        val genomeQuery = GenomeQuery("to1")

        val labelN = PeakAnnotation(Location(10, 20, chr1), PeakAnnotationType.NO_PEAKS)
        val labelP = PeakAnnotation(Location(30, 40, chr1), PeakAnnotationType.PEAKS)
        val labelS = PeakAnnotation(Location(50, 60, chr1), PeakAnnotationType.PEAK_START)
        val labelE = PeakAnnotation(Location(70, 80, chr1), PeakAnnotationType.PEAK_END)
        val labels = listOf(labelN, labelP, labelS, labelE)

        val extremePeaks = LocationsMergingList.create(genomeQuery,
                listOf(chr1.range.on(chr1, Strand.PLUS)))
        assertEquals(listOf(false, true, false, false), labels.map { it.check(extremePeaks) },
                "Wrong answer for test case 'extremePeaks'")

        val largePeaks = LocationsMergingList.create(genomeQuery,
                listOf(Location(35, 85, chr1)))
        assertEquals(listOf(true, true, false, false), labels.map { it.check(largePeaks) },
                "Wrong answer for test case 'largePeaks'")

        val smallPeaks = LocationsMergingList.create(genomeQuery,
                listOf(Location(31, 39, chr1), Location(55, 75, chr1)))
        assertEquals(listOf(true, true, true, true), labels.map { it.check(smallPeaks) },
                "Wrong answer for test case 'smallPeaks'")

        val noPeaks = LocationsMergingList.create(genomeQuery, emptyList())
        assertEquals(listOf(true, false, false, false), labels.map { it.check(noPeaks) },
                "Wrong answer for test case 'noPeaks'")

        val wrongChrPeaks = LocationsMergingList.create(genomeQuery,
                listOf(Location(35, 85, chrX)))
        assertEquals(listOf(true, false, false, false), labels.map { it.check(wrongChrPeaks) },
                "Wrong answer for test case 'wrongChrPeaks'")
    }

    @Test fun testLabelToBedEntry() {
        val chr1 = Chromosome("to1", "chr1")
        val annotationType = PeakAnnotationType.NO_PEAKS

        val labelN = PeakAnnotation(Location(10, 20, chr1), annotationType)
        assertEquals(ExtendedBedEntry("chr1", 10, 20, annotationType.id,
                                      itemRgb = annotationType.color.rgb),
                     labelN.asBedEntry())
    }

}