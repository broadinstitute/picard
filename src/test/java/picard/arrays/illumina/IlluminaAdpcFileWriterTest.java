package picard.arrays.illumina;

import htsjdk.samtools.util.IOUtil;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class IlluminaAdpcFileWriterTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_EXPECTED_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestIlluminaAdpcFileWriter.adpc.bin");

    @Test
    public void testWriteIlluminaAdpcFile() throws Exception {
        final File output = File.createTempFile("testIlluminaAdpcFileWriter.", ".adpc.bin");
        output.deleteOnExit();

        try (final IlluminaAdpcFileWriter adpcFileWriter = new IlluminaAdpcFileWriter(output)) {
            final List<IlluminaAdpcFileWriter.Record> adpcRecordList = new ArrayList<>();
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record((short) 11352, (short) 405, 1.444f, 0.088f, 0.705f, IlluminaGenotype.AA));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record((short) 458, (short) 2743, 0.043f, 0.852f, 0.818f, IlluminaGenotype.BB));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record((short) 7548, (short) 303, 1.072f, 0.076f, 0.0f, IlluminaGenotype.NN));
            adpcRecordList.add(new IlluminaAdpcFileWriter.Record((short) 7414, (short) 2158, 0.805f, 0.597f, 0.881f, IlluminaGenotype.AB));
            adpcFileWriter.write(adpcRecordList);
        }
        IOUtil.assertFilesEqual(TEST_EXPECTED_ADPC_BIN_FILE, output);
    }
}
