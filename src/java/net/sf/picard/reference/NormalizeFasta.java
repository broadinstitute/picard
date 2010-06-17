package net.sf.picard.reference;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.CloserUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 * Little program to "normalize" a fasta file to ensure that all line of sequence are the
 * same length, and are a reasonable length!
 */
public class NormalizeFasta extends CommandLineProgram {
    @Usage public final String USAGE = "Takes any file that conforms to the fasta format and " +
            "normalizes it so that all lines of sequence except the last line per named sequence " +
            "are of the same length.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input fasta file to normalize.")
    public File INPUT;

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output fasta file to write.")
    public File OUTPUT;

    @Option(doc="The line length to be used for the output fasta file.")
    public int LINE_LENGTH=100;

    @Option(doc="Truncate sequence names at first whitespace.")
    public boolean TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=false;

    public static void main(final String[] args) {
        new NormalizeFasta().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        if (INPUT.getAbsoluteFile().equals(OUTPUT.getAbsoluteFile())) {
            throw new IllegalArgumentException("Input and output cannot be the same file.");
        }

        final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(INPUT, TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE);
        final BufferedWriter out = IoUtil.openFileForBufferedWriting(OUTPUT);

        ReferenceSequence seq = null;
        while ((seq = ref.nextSequence()) != null) {
            final String name  = seq.getName();
            final byte[] bases = seq.getBases();

            try {
                out.write(">");
                out.write(name);
                out.newLine();

                for (int i=0; i<bases.length; ++i) {
                    if (i > 0 && i % LINE_LENGTH == 0) out.write("\n");
                    out.write(bases[i]);
                }

                out.write("\n");
            }
            catch (IOException ioe) {
                throw new PicardException("Error writing to file " + OUTPUT.getAbsolutePath(), ioe);

            }
        }

        CloserUtil.close(out);
        return 0;
    }
}
