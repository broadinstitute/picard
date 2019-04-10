package picard.vcf;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

public class SamTestUtils {
    /**
     * Useful test method.  Creates a (temporary) indexed BAM or CRAM so that we can store a sam in the testdata set.
     *
     * @param samFile the sam file to convert to bam and index
     * @param type the type of file to be created
     * @param referenceFasta reference to be used if creating CRAM
     * @return File a (temporary) bam or cram file (index file is created in same path).
     */

    public static File createIndexedBamOrCram(final File samFile, final File tempFilePrefix, final SamReader.Type type, final File referenceFasta) throws IOException {
        if (type != SamReader.Type.BAM_TYPE && type != SamReader.Type.CRAM_TYPE) {
            throw new PicardException("Cannot create indexed file of type " + type);
        }

        if (type == SamReader.Type.CRAM_TYPE && referenceFasta == null) {
            throw new PicardException("Cannot create indexed CRAM file without reference");
        }

        final File output = File.createTempFile(tempFilePrefix.getAbsolutePath(), "." + type.fileExtension());
        output.deleteOnExit();
        final File indexFile; //SamReaderFactor uses diffefent conventions for naming Bam indexes and Cram indexes, need to handle each case separately
        if (type == SamReader.Type.BAM_TYPE) {
            final String fileBase = output.getAbsolutePath().substring(0, output.getAbsolutePath().lastIndexOf('.'));
            indexFile = new File(fileBase + BAMIndex.BAI_INDEX_SUFFIX);
            indexFile.deleteOnExit();

        } else if (type == SamReader.Type.CRAM_TYPE) {
            indexFile = IOUtil.addExtension(Paths.get(output.getAbsolutePath()), BAMIndex.BAI_INDEX_SUFFIX).toFile();
            indexFile.deleteOnExit();

        }

        try (final SamReader in = SamReaderFactory.makeDefault().open(samFile);
             final SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(true).makeWriter(in.getFileHeader(), true, output, referenceFasta);) {
            in.iterator().stream().forEach(out::addAlignment);
        }

        return output;
    }

    public static File createIndexedBamOrCram(final File samFile, final File tempFilePrefix, final SamReader.Type type) throws IOException {
        return createIndexedBamOrCram(samFile, tempFilePrefix, type, null);
    }
}
