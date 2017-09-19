package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller.CalledCopyRatioSegment.Call.AMPLIFICATION_CALL;
import static org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller.CalledCopyRatioSegment.Call.DELETION_CALL;
import static org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller.CalledCopyRatioSegment.Call.NEUTRAL_CALL;

public final class ReCapSegCallerUnitTest extends BaseTest {
    @Test
    public void testMakeCalls() {

        final String sampleName = "Sample";
        final List<SimpleInterval> intervals = new ArrayList<>();
        final List<Double> testData = new ArrayList<>();

        //add amplification targets
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 101 + i, 101 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(2.0));
        }
        //add deletion targets
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 201 + i, 201 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(0.5));
        }
        //add obviously neutral targets with some small spread
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 301 + i, 301 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(0.01 * (i - 5) + 1));
        }
        //add spread-out targets to a neutral segment (mean near zero)
        for (int i = 0; i < 10; i++) {
            final SimpleInterval interval = new SimpleInterval("chr", 401 + i, 401 + i);
            intervals.add(interval);
            testData.add(ParamUtils.log2(0.1 * (i - 5) + 1));
        }

        final RealMatrix denoisedCopyRatioValues = new Array2DRowRealMatrix(1, intervals.size());
        denoisedCopyRatioValues.setRow(0, testData.stream().mapToDouble(x -> x).toArray());
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(sampleName, intervals, denoisedCopyRatioValues);

        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(sampleName,
                Arrays.asList(
                        new CopyRatioSegment(new SimpleInterval("chr", 101, 110), 10, ParamUtils.log2(2.0)),   //amplification
                        new CopyRatioSegment(new SimpleInterval("chr", 201, 210), 10, ParamUtils.log2(0.5)),   //deletion
                        new CopyRatioSegment(new SimpleInterval("chr", 301, 310), 10, ParamUtils.log2(1)),     //neutral
                        new CopyRatioSegment(new SimpleInterval("chr", 401, 410), 10, ParamUtils.log2(1))));   //neutral

        final CalledCopyRatioSegmentCollection calledCopyRatioSegments = new ReCapSegCaller(denoisedCopyRatios, copyRatioSegments).makeCalls();

        Assert.assertEquals(copyRatioSegments.getSampleName(), calledCopyRatioSegments.getSampleName());
        Assert.assertEquals(copyRatioSegments.getIntervals(), calledCopyRatioSegments.getIntervals());
        Assert.assertEquals(
                copyRatioSegments.getSegments().stream().map(CopyRatioSegment::getMeanLog2CopyRatio).collect(Collectors.toList()),
                calledCopyRatioSegments.getSegments().stream().map(CopyRatioSegment::getMeanLog2CopyRatio).collect(Collectors.toList()));
        Assert.assertEquals(
                calledCopyRatioSegments.getSegments().stream().map(CalledCopyRatioSegment::getCall).collect(Collectors.toList()),
                Arrays.asList(AMPLIFICATION_CALL, DELETION_CALL, NEUTRAL_CALL, NEUTRAL_CALL));
    }
}