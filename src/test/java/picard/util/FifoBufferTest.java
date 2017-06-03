package picard.util;

import htsjdk.samtools.util.Md5CalculatingInputStream;
import htsjdk.samtools.util.Md5CalculatingOutputStream;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;

/**
 * Tests for the FifoBuffer class. Generates (deterministically!) random files of certain sizes
 * and then pipes them through the FifoBuffer class and checks that input and output MD5s match.
 */
public class FifoBufferTest {

    /* Simply invokes the real test method many times with different inputs. */
    @Test public void testFifoBuffer() throws IOException {
        test(1);
        test(2);
        test(3);
        test(4);
        test(5);
        test(6);
        test(7);
        test(8);
        test(9);
        test(10);
        test(10.1345);
        test(150);
    }

    /**
     * Writes a file of the given size in megabytes fill with random bytes, initialized using (int) megabytes as the
     * seed to the random number generator. Then pipes through FifoBuffer and confirms that output == input via md5.
     */
    public void test(final double megabytes) throws IOException {
        final File inputFile  = File.createTempFile("fifo_input.", ".foo");
        inputFile.deleteOnExit();
        // Generate a file with a set number of megabytes of random data
        final int nBytes = (int) (megabytes * 1024 * 1024);
        {
            final Random random = new Random(nBytes);
            final BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(inputFile));
            final int bytesPerWrite = 127;
            final byte[] bytes = new byte[bytesPerWrite];
            for (int i=0; i<nBytes; i+=bytesPerWrite) {
                random.nextBytes(bytes);
                out.write(bytes);
            }
            out.close();
        }

        // Run the input file through the FifoBuffer and back into another file
        final Md5CalculatingInputStream in   = new Md5CalculatingInputStream(new FileInputStream(inputFile), (File)null);
        final Md5CalculatingOutputStream out = new Md5CalculatingOutputStream(new FileOutputStream("/dev/null"), (File)null);
        final PrintStream outStream = new PrintStream(out);
        final FifoBuffer buffer = new FifoBuffer(in, outStream);
        buffer.BUFFER_SIZE = 128 * 1024 * 1024;
        buffer.doWork();

        in.close();
        out.close();
        final String inputMd5  = in.md5();
        final String outputMd5 = out.md5();
        Assert.assertEquals(outputMd5, inputMd5, "MD5s do not match between input and output.");
        inputFile.delete();
    }
}
