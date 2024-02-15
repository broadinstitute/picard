package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;

public class AsIsCollectDuplicateMetricsTest extends AsIsMarkDuplicatesTest{
    @Override
    AbstractMarkDuplicatesCommandLineProgramTester getTester(SAMFileHeader.SortOrder so) {
        return new CollectDuplicateMetricsTester(so);
    }
}
