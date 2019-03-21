package picard.vcf;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import picard.PicardException;

import java.io.File;
import java.io.IOException;

public class SamTestUtils {
    /**
     * Useful test method.  Creates a (temporary) indexed BAM or CRAM so that we can store a sam in the testdata set.
     *
     * @param samFile the sam file to convert to bam and index
     * @return File a (temporary) bam or cram file (index file is created in same path).
     */
    public static File createIndexedBamOrCram(final File samFile, final File tempFilePrefix, final SamReader.Type type) throws IOException {
        return createIndexedBamOrCram(samFile, tempFilePrefix, type, null);
    }

    public static File createIndexedBamOrCram(final File samFile, final File tempFilePrefix, final SamReader.Type type, final File referenceFasta) throws IOException {
        if (type != SamReader.Type.BAM_TYPE && type != SamReader.Type.CRAM_TYPE) {
            throw new PicardException("Cannot create indexed file of type "+type);
        }

        if (type == SamReader.Type.CRAM_TYPE && referenceFasta == null) {
            throw new PicardException("Cannot create indexed CRAM file without reference");
        }

        final File output = File.createTempFile(tempFilePrefix.getAbsolutePath(), "."+type.fileExtension());
        output.deleteOnExit();
        final File indexFile = new File(output.getAbsolutePath() + "."+type.indexExtension());
        indexFile.deleteOnExit();

        final SamReader in = SamReaderFactory.makeDefault().open(samFile);
        SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(in.getFileHeader(), true, output,referenceFasta);

        in.iterator().stream().forEach(out::addAlignment);
        out.close();
        in.close();

        return output;
    }
}
