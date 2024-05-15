package picard.fingerprint;

import htsjdk.samtools.util.*;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.nio.PicardHtsPath;
import picard.util.IntervalListTools;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Convert Attempts to modify the representative Snp in the HaplotypeDatabase to overlap with" +
                "an input VCF",
        oneLineSummary = "Aligns a Haplotype database file with vcf",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class AlignHaplotypeDatabase extends CommandLineProgram {

    @Argument(doc = "Haplotype database to be adjusted.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Where to write converted haplotype database VCF.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "VCF with which to align representative SNP", shortName = "V")
    public File VCF;

    @Override
    protected boolean requiresReference() {
        return false;
    }


    @Override
    protected int doWork() {

        final HaplotypeMap haplotypeMap = new HaplotypeMap(INPUT);

        final IntervalList intervals = VCFFileReader.toIntervalList(VCF.toPath(), false);
        final OverlapDetector<Interval> detector = new OverlapDetector<>(0, 0);
        detector.addAll(intervals.getIntervals(), intervals.getIntervals());

        final HaplotypeMap newHaplotypeMap = new HaplotypeMap(haplotypeMap.getHeader());

        haplotypeMap.getHaplotypes().forEach(block-> {
            Optional<Snp> overlappingSnp = block.getSnps().stream().filter(detector::overlapsAny).findAny();
            if (overlappingSnp.isEmpty()){
                newHaplotypeMap.addHaplotype(block);
            }else{
                final HaplotypeBlock newBlock = new HaplotypeBlock(block.getMaf());
                newBlock.
            }
        );
        return 0;
    }

}