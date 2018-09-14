package picard.sam;

import com.sun.tools.javac.util.List;
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
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

    protected CigarOperator unscramCigarOp(CigarOperator op) {
        switch(op) {
            case EQ:
            case X:
                return CigarOperator.M;
            default:
                return op;
        }
    }

    //Turns EQ and X back into M.
    protected Cigar unscramCigar(Cigar c) {
        return Cigar.fromCigarOperators(
                c.getCigarElements().stream()
                    .flatMap( elem -> Stream.generate(() -> unscramCigarOp(elem.getOperator())).limit(elem.getLength()) )
                    .collect(Collectors.toList()));
    }

    protected boolean properCigarCompare(SAMRecord scramRec, SAMRecord bamRec) {
        return unscramCigar(scramRec.getCigar()).equals(unscramCigar(bamRec.getCigar()));
    }

    protected int set(int value, int flag) {
        return value | flag;
    }

    protected int unset(int value, int flag) {
        return value & (~flag);
    }

    protected boolean properFlagCompare(SAMRecord scramRec, SAMRecord bamRec) {
        //if mapQ = 0, it's okay to say that the scram is a secondary if the bam says supplementary.
        //we'll edit the flags on the bam to reflect this.
        int bamFlags = bamRec.getFlags();
        if( bamRec.getMappingQuality() == 0 && bamRec.getSupplementaryAlignmentFlag() ) {
            bamFlags = unset(bamFlags, SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue());
            bamFlags = set(bamFlags, SAMFlag.SECONDARY_ALIGNMENT.intValue());
        }

        return scramRec.getFlags() == bamFlags;
    }

    protected void compare(SAMRecord scramRec, SAMRecord bamRec) {
        if( !Objects.equals(scramRec.getContig(), bamRec.getContig()) ) {
            throw new ScramMismatchException(scramRec, bamRec, "contigs mismatch", scramRec.getContig(), bamRec.getContig() );
        } else if ( scramRec.getAlignmentStart() != bamRec.getAlignmentStart() ) {
            throw new ScramMismatchException(scramRec, bamRec, "alignment start mismatch", scramRec.getAlignmentStart(), bamRec.getAlignmentStart() );
        } else if ( scramRec.getAlignmentEnd() != bamRec.getAlignmentEnd() ) {
            throw new ScramMismatchException(scramRec, bamRec, "alignment end mismatch", scramRec.getAlignmentEnd(), bamRec.getAlignmentEnd() );
        } else if ( !Objects.equals(scramRec.getReadString(), bamRec.getReadString()) ) {
            throw new ScramMismatchException(scramRec, bamRec, "bases mismatch", scramRec.getReadString(), bamRec.getReadString() );
            //NOTE: no quals, because they're deoscillated.
        } else if ( scramRec.getMappingQuality() != bamRec.getMappingQuality() ) {
            throw new ScramMismatchException(scramRec, bamRec, "MAPQ mismatch", scramRec.getMappingQuality(), bamRec.getMappingQuality());
        } else if ( scramRec.getFlags() != bamRec.getFlags()
                && !properFlagCompare(scramRec, bamRec)) {
            //FIXME: delta supplementary shenanigans
            throw new ScramMismatchException(scramRec, bamRec, "flags mismatch", scramRec.getFlags(), bamRec.getFlags() );
        } else if ( !Objects.equals(scramRec.getCigarString(), bamRec.getCigarString())
                    && !properCigarCompare(scramRec, bamRec) ) {
            throw new ScramMismatchException(scramRec, bamRec, "cigar mismatch", scramRec.getCigarString(), bamRec.getCigarString());
        }

        //TODO: mate information
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
