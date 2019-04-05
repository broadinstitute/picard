package picard.sam;

import htsjdk.samtools.metrics.MetricBase;

public class SamComparisonMetric extends MetricBase {
    public String leftFile;
    public String rightFile;
    public int mappingsMatch;
    public int unmappedBoth;
    public int unmappedLeft;
    public int unmappedRight;
    public int mappingsDiffer;
    public int missingLeft;
    public int missingRight;
    public int duplicateMarkingsDiffer;
    public boolean areEqual;
}
