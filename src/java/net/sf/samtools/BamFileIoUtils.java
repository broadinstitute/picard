package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedFilePointerUtil;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.BlockCompressedStreamConstants;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.IOUtil;
import net.sf.samtools.util.Md5CalculatingOutputStream;
import net.sf.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

public class BamFileIoUtils {
    public static final String BAM_FILE_EXTENSION = ".bam";

    public static boolean isBamFile(final File file) {
        return ((file != null) && file.getName().endsWith(BAM_FILE_EXTENSION));
    }

    public static void reheaderBamFile(final SAMFileHeader samFileHeader, final File inputFile, final File outputFile) {
        reheaderBamFile(samFileHeader, inputFile, outputFile, true, true);
    }

    /**
     * Copy a BAM file but replacing the header
     * @param samFileHeader The header to use in the new file
     * @param inputFile The BAM file to copy, sans header
     * @param outputFile The new BAM file, constructed with the new header and the content from inputFile
     * @param createMd5 Whether or not to create an MD5 file for the new BAM
     * @param createIndex Whether or not to create an index file for the new BAM
     */
    public static void reheaderBamFile(final SAMFileHeader samFileHeader, final File inputFile, final File outputFile, final boolean createMd5, final boolean createIndex) {
        // TODO: In a future world where IoUtil and IOUtil are merged, de-comment these
//        IoUtil.assertFileIsReadable(inputFile);
//        IoUtil.assertFileIsWritable(outputFile);

        try {
            BlockCompressedInputStream.assertNonDefectiveFile(inputFile);
            assertSortOrdersAreEqual(samFileHeader, inputFile);

            final OutputStream outputStream = buildOutputStream(outputFile, createMd5, createIndex);

            BAMFileWriter.writeHeader(outputStream, samFileHeader);
            blockCopyBamFile(inputFile, outputStream, true, false);

            CloserUtil.close(inputFile);
            outputStream.close();
        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    /**
     * Copy data from a BAM file to an OutputStream by directly copying the gzip blocks
     *
     * @param inputFile The file to be copied
     * @param outputStream The stream to write the copied data to
     * @param skipHeader If true, the header of the input file will not be copied to the output stream
     * @param skipTerminator If true, the terminator block of the input file will not be written to the output stream
     */
    public static void blockCopyBamFile(final File inputFile, final OutputStream outputStream, final boolean skipHeader, final boolean skipTerminator) {
        FileInputStream in = null;
        try {
            in = new FileInputStream(inputFile);

            // a) It's good to check that the end of the file is valid and b) we need to know if there's a terminator block and not copy it if skipTerminator is true
            final BlockCompressedInputStream.FileTermination term = BlockCompressedInputStream.checkTermination(inputFile);
            if (term == BlockCompressedInputStream.FileTermination.DEFECTIVE)
                throw new SAMException(inputFile.getAbsolutePath() + " does not have a valid GZIP block at the end of the file.");

            if (skipHeader) {
                final long vOffsetOfFirstRecord = SAMUtils.findVirtualOffsetOfFirstRecordInBam(inputFile);
                final BlockCompressedInputStream blockIn = new BlockCompressedInputStream(inputFile);
                blockIn.seek(vOffsetOfFirstRecord);
                final long remainingInBlock = blockIn.available();

                // If we found the end of the header then write the remainder of this block out as a
                // new gzip block and then break out of the while loop
                if (remainingInBlock >= 0) {
                    final BlockCompressedOutputStream blockOut = new BlockCompressedOutputStream(outputStream, null);
                    IOUtil.transferByStream(blockIn, blockOut, remainingInBlock);
                    blockOut.flush();
                    // Don't close blockOut because closing underlying stream would break everything
                }

                long pos = BlockCompressedFilePointerUtil.getBlockAddress(blockIn.getFilePointer());
                blockIn.close();
                while (pos > 0) {
                    pos -= in.skip(pos);
                }
            }

            // Copy remainder of input stream into output stream
            final long currentPos = in.getChannel().position();
            final long length = inputFile.length();
            final long skipLast = ((term == BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK) && skipTerminator) ?
                    BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length : 0;
            final long bytesToWrite = length - skipLast - currentPos;

            IOUtil.transferByStream(in, outputStream, bytesToWrite);
        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        } finally {
            CloserUtil.close(in);
        }
    }

    /**
     * Assumes that all inputs and outputs are block compressed VCF files and copies them without decompressing and parsing
     * most of the gzip blocks. Will decompress and parse blocks up to the one containing the end of the header in each file
     * (often the first block) and re-compress any data remaining in that block into a new block in the output file. Subsequent
     * blocks (excluding a terminator block if present) are copied directly from input to output.
     */
    public static void gatherWithBlockCopying(final List<File> bams, final File output, final boolean createIndex, final boolean createMd5) {
        try {
            OutputStream out = new FileOutputStream(output);
            if (createMd5) out   = new Md5CalculatingOutputStream(out, new File(output.getAbsolutePath() + ".md5"));
            if (createIndex) out = new StreamInflatingIndexingOutputStream(out, new File(output.getParentFile(), IOUtil.basename(output) + BAMIndex.BAMIndexSuffix));
            boolean isFirstFile = true;

            for (final File f : bams) {
                blockCopyBamFile(f, out, !isFirstFile, true);
                isFirstFile = false;
            }

            // And lastly add the Terminator block and close up
            out.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
            out.close();
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    private static OutputStream buildOutputStream(final File outputFile, final boolean createMd5, final boolean createIndex) throws IOException {
        OutputStream outputStream = new FileOutputStream(outputFile);
        if (createMd5) outputStream   = new Md5CalculatingOutputStream(outputStream, new File(outputFile.getAbsolutePath() + ".md5"));
        if (createIndex) outputStream = new StreamInflatingIndexingOutputStream(outputStream, new File(outputFile.getParentFile(), IOUtil.basename(outputFile) + BAMIndex.BAMIndexSuffix));
        return outputStream;
    }

    private static void assertSortOrdersAreEqual(final SAMFileHeader newHeader, final File inputFile) throws IOException {
        final SAMFileHeader origHeader = new SAMFileReader(inputFile).getFileHeader();
        final SAMFileHeader.SortOrder newSortOrder = newHeader.getSortOrder();
        if (newSortOrder != SAMFileHeader.SortOrder.unsorted && newSortOrder != origHeader.getSortOrder()) {
            throw new SAMException("Sort order of new header does not match the original file, needs to be " + origHeader.getSortOrder());
        }
    }
}
