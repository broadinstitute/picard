package picard.sam;

import htsjdk.samtools.*;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Unit tests for FilterFlowCellEdgeReads using TestNG.
 */
public class FilterFlowCellEdgeReadsTest {

    // Temporary files for input and output.
    private File inputSam;
    private File outputSam;

    /**
     * Helper method to create a temporary SAM file with one record per provided read name.
     * Each record is given minimal required fields.
     *
     * @param readNames an array of read names to include.
     * @return the temporary SAM file.
     * @throws IOException if an I/O error occurs.
     */
    private File createSamFile(String[] readNames) throws IOException {
        File tmpSam = File.createTempFile("FilterFlowCellEdgeReadsTest_input", ".sam");
        tmpSam.deleteOnExit();

        // Create a minimal SAM file header.
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        // Add one sequence record so that records have a reference.
        header.addSequence(new SAMSequenceRecord("chr1", 1000000));

        // Use SAMFileWriterFactory to write a SAM file.
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, false, tmpSam, null)) {
            // Create one record for each read name.
            for (String readName : readNames) {
                SAMRecord rec = new SAMRecord(header);
                rec.setReadName(readName);
                rec.setReferenceName("chr1");
                rec.setAlignmentStart(1);
                rec.setCigarString("50M");
                // Set dummy bases and qualities.
                rec.setReadString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                rec.setBaseQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
                writer.addAlignment(rec);
            }
        }
        return tmpSam;
    }

    /**
     * Helper method to count the number of SAM records in a given SAM file.
     *
     * @param samFile the SAM file.
     * @return the number of records.
     * @throws IOException if an I/O error occurs.
     */
    private int countRecords(File samFile) throws IOException {
        int count = 0;
        try (SamReader reader = SamReaderFactory.makeDefault().open(samFile)) {
            for (SAMRecord rec : reader) {
                count++;
            }
        }
        return count;
    }

    @AfterMethod
    public void tearDown() {
        if (inputSam != null && inputSam.exists()) {
            inputSam.delete();
        }
        if (outputSam != null && outputSam.exists()) {
            outputSam.delete();
        }
    }

    /**
     * Test with a mixed input:
     * – One read with a name that matches the default coordinates ("1000:1000") and should be filtered out.
     * – One read with non-matching coordinates ("2000:2000") that should be retained.
     */
    @Test
    public void testMixedReads() throws IOException {
        String[] readNames = new String[]{
                "EAS139:136:FC706VJ:2:1000:1000", // should be filtered out (matches default X_COORD and Y_COORD)
                "EAS139:136:FC706VJ:2:2000:2000"  // should be retained
        };
        inputSam = createSamFile(readNames);
        outputSam = File.createTempFile("FilterFlowCellEdgeReadsTest_output", ".sam");
        outputSam.deleteOnExit();

        FilterFlowCellEdgeReads tool = new FilterFlowCellEdgeReads();
        tool.INPUT = inputSam.getAbsolutePath();
        tool.OUTPUT = outputSam.getAbsolutePath();
        // Use default X_COORD=1000, Y_COORD=1000

        int ret = tool.doWork();
        Assert.assertEquals(ret, 0, "doWork() should return 0");

        // Only the record that does not match the filter should be written.
        int recordCount = countRecords(outputSam);
        Assert.assertEquals(recordCount, 1, "Only one record should be written");
    }

    /**
     * Test with a read whose name does not contain colon-delimited coordinates.
     * The method hasFlowcellCoordinates should catch the exception and return false,
     * so the record should be retained.
     */
    @Test
    public void testNonConformingReadName() throws IOException {
        String[] readNames = new String[]{
                "nonconforming_read"  // no colon-separated parts → not filtered
        };
        inputSam = createSamFile(readNames);
        outputSam = File.createTempFile("FilterFlowCellEdgeReadsTest_output", ".sam");
        outputSam.deleteOnExit();

        FilterFlowCellEdgeReads tool = new FilterFlowCellEdgeReads();
        tool.INPUT = inputSam.getAbsolutePath();
        tool.OUTPUT = outputSam.getAbsolutePath();
        // Defaults are used.

        int ret = tool.doWork();
        Assert.assertEquals(ret, 0);

        // The read should be retained.
        int recordCount = countRecords(outputSam);
        Assert.assertEquals(recordCount, 1, "The nonconforming read should be kept");
    }

    /**
     * Test with an input that has only a read with coordinates matching the filter.
     * In this case, the tool should filter out the only record and write an empty output.
     */
    @Test
    public void testAllReadsFiltered() throws IOException {
        String[] readNames = new String[]{
                "EAS139:136:FC706VJ:2:1000:1000"  // matches filter → filtered out
        };
        inputSam = createSamFile(readNames);
        outputSam = File.createTempFile("FilterFlowCellEdgeReadsTest_output", ".sam");
        outputSam.deleteOnExit();

        FilterFlowCellEdgeReads tool = new FilterFlowCellEdgeReads();
        tool.INPUT = inputSam.getAbsolutePath();
        tool.OUTPUT = outputSam.getAbsolutePath();
        // Defaults: X_COORD=1000, Y_COORD=1000

        int ret = tool.doWork();
        Assert.assertEquals(ret, 0);

        // Expect zero records in the output.
        int recordCount = countRecords(outputSam);
        Assert.assertEquals(recordCount, 0, "No records should be written");
    }
}
