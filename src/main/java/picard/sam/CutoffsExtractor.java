package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static htsjdk.samtools.SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = ReadDataManipulationProgramGroup.class)
public class CutoffsExtractor extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to extract.")
    public File INPUT;

    @Argument(shortName = "OC", doc = "Cutoffs output file.")
    public File OUTPUT_CUTOFFS;

    @Argument(shortName = "OQ", doc = "Quals output file.")
    public File OUTPUT_QUALS;

    @Argument(shortName = "C", doc = "Cutoff threshold.")
    public Integer CUTOFF = 10;

    //return the index of the first item in the list that is above the cutoff
    protected static int getIndexAboveCutoff(List<Integer> ls, Integer cutoff) {
        for(int i = 0; i < ls.size(); ++i) {
            if( ls.get(i) > cutoff ) {
                return i;
            }
        }
        return ls.size();
    }

    private int[] byte2intQuals(byte[] bquals) {
        int[] iquals = new int[bquals.length];
        for (int i = 0; i < bquals.length; ++i) {
            iquals[i] = Byte.toUnsignedInt(bquals[i]);
        }
        return iquals;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT_CUTOFFS);
        IOUtil.assertFileIsWritable(OUTPUT_QUALS);

        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().setOption(INCLUDE_SOURCE_IN_RECORDS, true).validationStringency(ValidationStringency.SILENT);

        //TODO: writers

        try (final SamReader samReader = readerFactory.open(SamInputResource.of(INPUT)) ) {
            try (BufferedWriter cutoffsW = new BufferedWriter(new FileWriter(OUTPUT_CUTOFFS))) {
                try (BufferedWriter qualsW = new BufferedWriter(new FileWriter(OUTPUT_QUALS))) {
                    for (SAMRecord record : samReader) {

                        int[] iquals = byte2intQuals(record.getBaseQualities());
                        List<Integer> lquals = Arrays.stream(iquals).boxed().collect(Collectors.toList());

                        //flip so the bad-quals are at the beginning
                        if(!record.getReadNegativeStrandFlag()) {
                            Collections.reverse(lquals);
                        }

                        int idx = getIndexAboveCutoff(lquals, CUTOFF);
                        int avgQual = (int) Math.round(lquals.subList(idx, iquals.length).stream().mapToInt(i->i).average().orElse(0));

                        //write how many quals are GOOD
                        cutoffsW.write(iquals.length - idx);
                        qualsW.write(avgQual);
                    }
                }
            }
        } catch (IOException e) {
            System.out.println( "A problem occurred: " + e.getMessage());
            return 1;
        }

        return 0;
    }
}
