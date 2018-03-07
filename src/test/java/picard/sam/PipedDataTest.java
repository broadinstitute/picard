package picard.sam;

import htsjdk.samtools.*;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

public class PipedDataTest extends CommandLineProgramTest {


    @Override
    public String getCommandLineProgramName() {
        return SortSam.class.getSimpleName();
    }

    private static final File BAM = new File("testdata/picard/sam/test.bam");


    @Test
    public void testPipedData() throws IOException, InterruptedException, ExecutionException {

        ExecutorService exec = Executors.newCachedThreadPool();

        Path tmpFile = Files.createTempFile("sortSamPipeTest", "tmp");
        Path tmpOutputFile = Files.createTempFile("sortSamPipeTest", "output");

        // Read the sam file and do some simple operation
        Callable<Integer> readSam = () -> {
            SamReader reader = SamReaderFactory.makeDefault().open(BAM);
            final SAMFileWriter outputSam = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(),
                    true, tmpFile);
            for (SAMRecord record : reader) {
                record.setReadName(record.getReadName().toUpperCase());
                outputSam.addAlignment(record);
            }
            outputSam.close();
            return 0;
        };

        // Read the sam file output by htsjdk tool.
        Callable<Integer> sortSam = () -> {
            final String[] args = new String[]{
                    "INPUT=" + tmpFile.toAbsolutePath().toString(),
                    "OUTPUT=" + tmpOutputFile.toAbsolutePath().toString(),
                    "SORT_ORDER=queryname"
            };
            return runPicardCommandLine(args);
        };

        // Run them at the same time to simulate piping and waiting for EOF
        List<Future<Integer>> futures = exec.invokeAll(Arrays.asList(readSam, sortSam));

        Assert.assertEquals(futures.get(1).get(), new Integer(0));
    }
}
