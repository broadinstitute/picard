package picard.vcf;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;

public class SamTestUtils {
    /**
     * Useful test method.  Creates a (temporary) indexed BAM so that we can store a sam in the testdata set.
     *
     * @param samFile the sam file to convert to bam and index
     * @return File a (temporary) bam file (index file is created in same path).
     */
    public static File createIndexedBam(final File samFile, final File tempFilePrefix) throws IOException {
        final File output = File.createTempFile(tempFilePrefix.getAbsolutePath(), ".bam");
        output.deleteOnExit();
        final File indexFile = new File(output.getAbsolutePath() + ".bai");
        indexFile.deleteOnExit();

        final SamReader in = SamReaderFactory.makeDefault().open(samFile);
        SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(in.getFileHeader(), true, output);

        in.iterator().stream().forEach(out::addAlignment);
        out.close();
        in.close();

        return output;
    }
}
