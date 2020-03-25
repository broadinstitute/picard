package picard.arrays.illumina;

import htsjdk.samtools.util.IOUtil;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


public class IlluminaAdpcFileWriterTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_EXPECTED_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestVcfToAdpc.adpc.bin");

    @Test
    public void testWriteIlluminaAdpcFile() throws Exception {
        final File output = File.createTempFile("testIlluminaAdpcFileWriter.", ".adpc.bin");
        output.deleteOnExit();

        try (final IlluminaAdpcFileWriter adpcFileWriter = new IlluminaAdpcFileWriter(output)) {
            final List<IlluminaAdpcFileWriter.Record> adpcRecordList = new ArrayList<>();
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record(11352, 405, 1.444f, 0.088f, 0.705f, IlluminaGenotype.AA));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record(458, 2743, 0.043f, 0.852f, 0.818f, IlluminaGenotype.BB));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record(7548, 303, 1.072f, 0.076f, 0.0f, IlluminaGenotype.NN));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record(7414, 2158, 0.805f, 0.597f, 0.881f, IlluminaGenotype.AB));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record(222, 215, 0.0f, 0.0f, 0.91f, IlluminaGenotype.NN));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record(232, 246, null, null, 0.926f, IlluminaGenotype.NN));
            adpcFileWriter.write(adpcRecordList);
        }
        IOUtil.assertFilesEqual(TEST_EXPECTED_ADPC_BIN_FILE, output);
    }
}
