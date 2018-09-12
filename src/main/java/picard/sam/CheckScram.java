package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import static htsjdk.samtools.SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS;

@CommandLineProgramProperties(
        summary = "Compare SCRAM to its source bam",
        oneLineSummary = "compare SCRAM to its sorce bam",
        programGroup = DiagnosticsAndQCProgramGroup.class)
@DocumentedFeature
public class CheckScram extends CommandLineProgram {
    @Argument(shortName = "S", doc = "The scram file.")
    public File SCRAM;

    @Argument(shortName = "B", doc = "The bam file.")
    public File BAM;

    String safeToString(Object o) {
        if(o == null) {
            return "null";
        } else {
            return o.toString();
        }
    }

    class ScramMismatchException extends RuntimeException {
        public ScramMismatchException(SAMRecord scram, SAMRecord bam, String why, Object scramval, Object bamval) {
            super(  scram.toString() + " doesn't match " + bam.toString() + " because " + why + ":\n" +
                    safeToString(scramval) + " -- scram\n" +
                    safeToString(bamval) + " -- bam\n\n" +
                    scram.getSAMString() + " -- scram\n" +
                    bam.getSAMString() + " -- bam\n");
        }
    }

    protected void compare(SAMRecord scramRec, SAMRecord bamRec) {
        if( !scramRec.getContig().equals(bamRec.getContig()) ) {
            throw new ScramMismatchException(scramRec, bamRec, "contigs mismatch", scramRec.getContig(), bamRec.getContig() );
        } else if ( scramRec.getAlignmentStart() != bamRec.getAlignmentStart() ) {
            throw new ScramMismatchException(scramRec, bamRec, "alignment start mismatch", scramRec.getAlignmentStart(), bamRec.getAlignmentStart() );
        } else if ( scramRec.getAlignmentEnd() != bamRec.getAlignmentEnd() ) {
            throw new ScramMismatchException(scramRec, bamRec, "alignment end mismatch", scramRec.getAlignmentEnd(), bamRec.getAlignmentEnd() );
        } else if ( !scramRec.getReadString().equals(bamRec.getReadString()) ) {
            throw new ScramMismatchException(scramRec, bamRec, "bases mismatch", scramRec.getReadString(), bamRec.getReadString() );
            //NOTE: no quals, because they're deoscillated.
        } else if ( scramRec.getFlags() != bamRec.getFlags()) {
            //FIXME: delta supplementary shenanigans
            throw new ScramMismatchException(scramRec, bamRec, "flags mismatch", scramRec.getFlags(), bamRec.getFlags() );
        } else if ( scramRec.getMappingQuality() != bamRec.getMappingQuality() ) {
            throw new ScramMismatchException(scramRec, bamRec, "MAPQ mismatch", scramRec.getMappingQuality(), bamRec.getMappingQuality());
        } else if ( !scramRec.getCigarString().equals(bamRec.getCigarString()) ) {
            //TODO: M-ify both cigars
            throw new ScramMismatchException(scramRec, bamRec, "cigar mismatch", scramRec.getCigarString(), bamRec.getCigarString());
        }
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(SCRAM);
        IOUtil.assertFileIsReadable(BAM);

        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().setOption(INCLUDE_SOURCE_IN_RECORDS, true).validationStringency(ValidationStringency.SILENT);

        try (final SamReader scramReader = readerFactory.open(SamInputResource.of(SCRAM)) ) {
            try (final SamReader bamReader = readerFactory.open(SamInputResource.of(BAM)) ) {
                Iterator<SAMRecord> scramit = scramReader.iterator();
                Iterator<SAMRecord> bamit = bamReader.iterator();
                while(true) {
                    boolean scramNext = scramit.hasNext();
                    boolean bamNext = bamit.hasNext();

                    if(scramNext != bamNext) {
                        throw new RuntimeException("Iterators fell out of sync! This means there is a different number of records.");
                    }

                    if(!scramNext) {
                        break;
                    }

                    SAMRecord scramRec = scramit.next();
                    SAMRecord bamRec = bamit.next();

                    compare(scramRec, bamRec);
                }
            }
        } catch (IOException e) {
            System.out.println( "A problem occurred: " + e.getMessage());
            return 1;
        }
        return 0;
    }
}
