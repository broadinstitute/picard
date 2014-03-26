package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.util.BlockCompressedFilePointerUtil;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.BlockCompressedStreamConstants;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.Md5CalculatingOutputStream;
import net.sf.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.List;

/**
 * Program to perform a rapid "gather" operation on BAM files after a scatter operations where
 * the same process has been performed on different regions of a BAM file creating many smaller
 * BAM files that now need to be concatenated back together.
 *
 * @author Tim Fennell
 */
public class GatherBamFiles extends CommandLineProgram {
    @Usage public final String USAGE = "Concatenates one or more BAM files together as efficiently as possible. Assumes that the " +
            "list of BAM files provided as INPUT are in the order that they should be concatenated and simply concatenates the bodies " +
            "of the BAM files while retaining the header from the first file.  Operates via copying of the gzip blocks directly for speed " +
            "but also supports generation of an MD5 on the output and indexing of the output BAM file. Only support BAM files, does not " +
            "support SAM files.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc="One or more BAM files or text files containing lists of BAM files one per line.")
    public List<File> INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output BAM file to write.")
    public File OUTPUT;

    @Option(doc="Whether to use low-level gzip block copying for performance. The block-copying method is usually 3-5 times faster and " +
            "is recommended. Non block-copying is supported primarily for testing.")
    public boolean BLOCK_COPY = true;

    private static final Log log = Log.getInstance(GatherBamFiles.class);

    // Stock main method.
    public static void main(final String[] args) {
        final GatherBamFiles gatherer = new GatherBamFiles();
        gatherer.CREATE_INDEX = true;
        gatherer.instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        final List<File> inputs = IoUtil.unrollFiles(INPUT, ".bam", ".sam");
        for (final File f: inputs) IoUtil.assertFileIsReadable(f);
        IoUtil.assertFileIsWritable(OUTPUT);

        if (BLOCK_COPY) {
            gatherWithBlockCopying(inputs, OUTPUT, CREATE_INDEX, CREATE_MD5_FILE);
        }
        else {
            gatherNormally(inputs, OUTPUT, CREATE_INDEX, CREATE_MD5_FILE);
        }

        return 0;
    }

    /**
     * Simple implementation of a gather operations that uses SAMFileReaders and Writers in order to concatenate
     * multiple BAM files.
     */
    private static void gatherNormally(final List<File> inputs, final File output, final boolean createIndex, final boolean createMd5) {
        final SAMFileHeader header;
        {
            final SAMFileReader tmp = new SAMFileReader(inputs.get(0));
            header = tmp.getFileHeader();
            tmp.close();
        }

        final SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(createIndex).setCreateMd5File(createMd5).makeSAMOrBAMWriter(header, true, output);

        for (final File f : inputs) {
            log.info("Gathering " + f.getAbsolutePath());
            final SAMFileReader in = new SAMFileReader(f);
            for (final SAMRecord rec : in) out.addAlignment(rec);
            CloserUtil.close(in);
        }

        out.close();
    }

    /**
     * Assumes that all inputs and outputs are block compressed VCF files and copies them without decompressing and parsing
     * most of the gzip blocks. Will decompress and parse blocks up to the one containing the end of the header in each file
     * (often the first block) and re-compress any data remaining in that block into a new block in the output file. Subsequent
     * blocks (excluding a terminator block if present) are copied directly from input to output.
     */
    private static void gatherWithBlockCopying(final List<File> bams, final File output, final boolean createIndex, final boolean createMd5) {
        try {
            OutputStream out = new FileOutputStream(output);
            if (createMd5) out   = new Md5CalculatingOutputStream(out, new File(output.getAbsolutePath() + ".md5"));
            if (createIndex) out = new IndexingOutputStream(out, new File(output.getParentFile(), IoUtil.basename(output) + ".bai"));
            boolean isFirstFile = true;

            for (final File f : bams) {
                log.info("Gathering " + f.getAbsolutePath());
                final FileInputStream in = new FileInputStream(f);

                // a) It's good to check that the end of the file is valid and b) we need to know if there's a terminator block and not copy it
                final BlockCompressedInputStream.FileTermination term = BlockCompressedInputStream.checkTermination(f);
                if (term == BlockCompressedInputStream.FileTermination.DEFECTIVE) throw new PicardException(f.getAbsolutePath() + " does not have a valid GZIP block at the end of the file.");

                if (!isFirstFile) {
                    final long vOffsetOfFirstRecord = SAMUtils.findVirtualOffsetOfFirstRecordInBam(f);
                    final BlockCompressedInputStream blockIn = new BlockCompressedInputStream(f);
                    blockIn.seek(vOffsetOfFirstRecord);
                    final long remainingInBlock = blockIn.available();

                    // If we found the end of the header then write the remainder of this block out as a
                    // new gzip block and then break out of the while loop
                    if (remainingInBlock >= 0) {
                        final BlockCompressedOutputStream blockOut = new BlockCompressedOutputStream(out, null);
                        IoUtil.transferByStream(blockIn, blockOut, remainingInBlock);
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
                final long length     = f.length();
                final long skipLast   = (term == BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK) ?
                        BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length : 0;
                final long bytesToWrite = length - skipLast - currentPos;

                IoUtil.transferByStream(in, out, bytesToWrite);
                in.close();
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
}

/**
 * A little class that takes an InputStream from which it reads a BAM file, generates
 * a BAMIndex and then writes the index to the File provided.  All operations are designed
 * to be carried out in a separate thread.
 */
class Indexer implements Runnable {
    private final File index;
    private final InputStream stream;

    /** Constructs an indexer that reads from the stream provided and writes an index to the File provided. */
    Indexer(final File index, final InputStream stream) {
        this.index = index;
        this.stream = stream;
    }

    /** Runnable implementation that reads the entire stream and writes the index. */
    @Override
    public void run() {
        final SAMFileReader in = new SAMFileReader(this.stream);
        in.enableFileSource(true);
        in.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        in.enableCrcChecking(false);
        final BAMIndexer indexer = new BAMIndexer(this.index, in.getFileHeader());
        for (final SAMRecord rec : in) {
            indexer.processAlignment(rec);
        }

        indexer.finish();
        in.close();
    }
}

/**
 * OutputStream implementation that writes output to an underlying output stream while also copying the
 * same bytes to a PipedOutputStream that routes the data back into an Indexer to generate a BAMIndex
 * by inflating and decoding the stream and feeding the SAMRecords to a BAMIndexer.
 */
class IndexingOutputStream extends OutputStream {
    private final OutputStream s1;
    private final PipedOutputStream s2;
    private final Thread thread;

    IndexingOutputStream(final OutputStream s1, final File indexFile) {
        try {
            this.s1 = s1;
            this.s2 = new PipedOutputStream();
            final PipedInputStream pin = new PipedInputStream(this.s2, Defaults.NON_ZERO_BUFFER_SIZE);
            this.thread = new Thread(new Indexer(indexFile, pin), "BamIndexingThread");
            this.thread.start();
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    @Override public final void write(final int b) throws IOException { this.s1.write(b); this.s2.write(b); }
    @Override public final void write(final byte[] b) throws IOException { this.s1.write(b); this.s2.write(b); }
    @Override public final void write(final byte[] b, final int off, final int len) throws IOException { this.s1.write(b, off, len); this.s2.write(b, off, len); }
    @Override public final void flush() throws IOException { this.s1.flush(); this.s2.flush(); }

    @Override public final void close() throws IOException {
        this.s1.close();
        this.s2.close();

        try { this.thread.join(); }
        catch (final InterruptedException ie) { throw new RuntimeException(ie); }
    }
}