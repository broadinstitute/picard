package picard.sam;

import htsjdk.samtools.SAMFormatException;
import org.apache.commons.io.FileUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;

public class SortSamTest {
    @Test
    public void bugTest() throws Exception {
        File input = File.createTempFile("testIn",".bam");
        File output = File.createTempFile("testOut",".bam");
        FileUtils.write(input, "not valid sam input");
        try {
            new SortSam().instanceMain(new String[]{
                    "I=" + input.getPath(),
                    "O=" + output.getPath(),
                    "SORT_ORDER=coordinate"});
        } catch (SAMFormatException ex) {}
        Files.delete(output.toPath());
    }
}