/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.vcf.MendelianViolations;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;
import picard.pedigree.PedFile;
import picard.vcf.processor.VariantProcessor;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import static htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY;

/**
 * Finds mendelian violations of all types within a VCF.  Takes in VCF or BCF and a pedigree file and looks
 * for high confidence calls where the genotype of the offspring is incompatible with the genotypes of the
 * parents.  Key features:
 * 1) Checks for regular MVs in diploid regions and invalid transmissions in haploid regions
 * 2) Outputs metrics about the different kinds of MVs found
 * 3) Can output a per-trio VCF with violations; INFO field will contain a MV= attribute with the type of violation
 *
 * This CLP ignores variants that are:
 * - Not SNPs
 * - Filtered
 * - Multiallelic
 * - Monomorphic
 * - Within the SKIP_CHROMS contigs
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Finds mendelian violations of all types within a VCF. " +
           "Takes in VCF or BCF and a pedigree file and looks for high confidence calls " +
           "where the genotype of the offspring is incompatible with the genotypes of the parents. " +
                "Assumes the existence of format fields AD, DP, GT, GQ, and PL fields. " +
                "\n" +
                "Take note that the implementation assumes that reads from the PAR will be mapped to the female chromosome" +
                "rather than the male. This requires that the PAR in the male chromosome be masked so that the aligner " +
                "has a single coting to map to. This is normally done for the public releases of the human reference."+
                "\n" +
                "Usage example: java -jar picard.jar FindMendelianViolations I=input.vcf \\\n" +
                "                 TRIO=family.ped \\\n" +
                "                 OUTPUT=mendelian.txt \\\n" +
                "                 MIN_DP=20 \n" +
                "\n"
        ,
        oneLineSummary = "Finds mendelian violations of all types within a VCF",
        programGroup = VcfOrBcf.class
)
public class FindMendelianViolations extends CommandLineProgram {
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input VCF or BCF with genotypes.")
    public File INPUT;

    @Argument(shortName = "PED", doc = "File of Trio information in PED format (with no genotype columns).")
    public File TRIOS;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Argument(shortName = "GQ", doc = "Minimum genotyping quality (or non-ref likelihood) to perform tests.")
    public int MIN_GQ = 30;

    @Argument(shortName = "DP", doc="Minimum depth in each sample to consider possible mendelian violations.")
    public int MIN_DP = 0;

    @Argument(shortName = "MINHET", doc = "Minimum allele balance at sites that are heterozygous in the offspring.")
    public double MIN_HET_FRACTION = 0.3;

    @Argument(optional = true, doc = "If provided, output per-family VCFs of mendelian violations into this directory.")
    public File VCF_DIR;

    @Argument(doc = "List of chromosome names to skip entirely.")
    public Set<String> SKIP_CHROMS = CollectionUtil.makeSet("MT", "chrM");

    @Argument(doc = "List of possible names for male sex chromosome(s)")
    public Set<String> MALE_CHROMS = CollectionUtil.makeSet("Y", "chrY");

    @Argument(doc = "List of possible names for female sex chromosome(s)")
    public Set<String> FEMALE_CHROMS = CollectionUtil.makeSet("X", "chrX");

    @Argument(doc = "List of chr:start-end for pseudo-autosomal regions on the female sex chromosome. Defaults to HG19/b37 & HG38 coordinates.")
    public Set<String> PSEUDO_AUTOSOMAL_REGIONS = CollectionUtil.makeSet("X:10001-2649520", "X:59034050-59373566", "chrX:10000-2781479", "chrX:155701382-156030895");

    @Argument(doc = "The number of threads that will be used to collect the metrics. ")
    public int THREAD_COUNT = 1;

    @Argument(doc = "If true then fields need to be delimited by a single tab. If false the delimiter is one or more whitespace characters." +
            " Note that tab mode does not strictly follow the PED spec")
    public boolean TAB_MODE = false;

    private final Log LOG = Log.getInstance(FindMendelianViolations.class);
    private final ProgressLogger progressLogger = new ProgressLogger(LOG, 10000, "variants analyzed");

    // The following two members are "lazy" since they need to wait until the commandline program populates
    // the files from the arguments before they can be realized (while keeping them final)
    private final Lazy<VCFHeader> inputHeader = new Lazy<>(
            () -> {
                try (final VCFFileReader in = new VCFFileReader(INPUT)) {
                    return in.getFileHeader();
                }
            });

    private final Lazy<PedFile> pedFile = new Lazy<>(()-> PedFile.fromFile(TRIOS, TAB_MODE));

    private final Lazy<Set<Interval>> parIntervals = new Lazy<>(() -> Collections.unmodifiableSet(parseIntervalLists(PSEUDO_AUTOSOMAL_REGIONS)));

    private Set<Interval> parseIntervalLists(final Set<String> intervalLists){
        // Parse the PARs
        final Set<Interval> intervals = new HashSet<>(intervalLists.size());
        for (final String par : intervalLists) {
            final String[] splits1 = par.split(":");
            final String[] splits2 = splits1[1].split("-");
            intervals.add(new Interval(splits1[0], Integer.parseInt(splits2[0]), Integer.parseInt(splits2[1])));
        }
        return intervals;
    }
     /**
     * Validates that the sex chromosomes don't overlap and parses the pseudo-autosomal regions into usable
     * objects to ensure their parsability
     *
     * @return
     */
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();

        // Check that the sex chromosomes are not overlapping
        final Set<String> sexChroms = new HashSet<>();
        sexChroms.addAll(FEMALE_CHROMS);
        sexChroms.retainAll(MALE_CHROMS);
        if (!sexChroms.isEmpty()) errors.add("The following chromosomes were listed as both male and female sex chromosomes: " + sexChroms);

        if (errors.isEmpty()) return null;
        else return errors.toArray(new String[errors.size()]);
    }

    @Override
    protected int doWork() {
        ///////////////////////////////////////////////////////////////////////
        // Test and then load up the inputs
        ///////////////////////////////////////////////////////////////////////
        final boolean outputVcfs = VCF_DIR != null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(TRIOS);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (outputVcfs) IOUtil.assertDirectoryIsWritable(VCF_DIR);

        LOG.info("Loading and filtering trios.");

        final MendelianViolationDetector.Result result =
                VariantProcessor.Builder
                        .generatingAccumulatorsBy(this::buildDetector)
                        .withInput(INPUT)
                        .combiningResultsBy(MendelianViolationDetector.Result::merge)
                        .multithreadingBy(THREAD_COUNT)
                        .build()
                        .process();

        // Write out the metrics!
        final MetricsFile<MendelianViolationMetrics, ?> metricsFile = getMetricsFile();
        for (final MendelianViolationMetrics m : result.metrics()) {
            m.calculateDerivedFields();
            metricsFile.addMetric(m);
        }

        metricsFile.write(OUTPUT);
        writeAllViolations(result);

        return 0;
    }

    private void writeAllViolations(final MendelianViolationDetector.Result result) {
        if (VCF_DIR != null) {
            LOG.info(String.format("Writing family violation VCFs to %s/", VCF_DIR.getAbsolutePath()));

            final VariantContextComparator vcComparator = new VariantContextComparator(inputHeader.get().getContigLines());
            final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputHeader.get().getMetaDataInInputOrder());

            headerLines.add(new VCFInfoHeaderLine(MendelianViolationDetector.MENDELIAN_VIOLATION_KEY, 1, VCFHeaderLineType.String, "Type of mendelian violation."));
            headerLines.add(new VCFInfoHeaderLine(MendelianViolationDetector.ORIGINAL_AC, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Original AC"));
            headerLines.add(new VCFInfoHeaderLine(MendelianViolationDetector.ORIGINAL_AF, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Original AF"));
            headerLines.add(new VCFInfoHeaderLine(MendelianViolationDetector.ORIGINAL_AN, 1, VCFHeaderLineType.Integer, "Original AN"));

            for (final PedFile.PedTrio trio : pedFile.get().values()) {
                final File outputFile = new File(VCF_DIR, IOUtil.makeFileNameSafe(trio.getFamilyId() + ".vcf"));
                LOG.info(String.format("Writing %s violation VCF to %s", trio.getFamilyId(), outputFile.getAbsolutePath()));

                final VariantContextWriter out = new VariantContextWriterBuilder()
                        .setOutputFile(outputFile)
                        .unsetOption(INDEX_ON_THE_FLY)
                        .build();

                final VCFHeader newHeader = new VCFHeader(headerLines, CollectionUtil.makeList(trio.getMaternalId(), trio.getPaternalId(), trio.getIndividualId()));
                final TreeSet<VariantContext> orderedViolations = new TreeSet<>(vcComparator);

                orderedViolations.addAll(result.violations().get(trio.getFamilyId()));
                out.writeHeader(newHeader);
                orderedViolations.forEach(out::add);

                out.close();
            }
        }
    }

    private MendelianViolationDetector buildDetector() {
        return new MendelianViolationDetector(
                ImmutableSet.copyOf(SKIP_CHROMS),
                ImmutableSet.copyOf(MALE_CHROMS),
                ImmutableSet.copyOf(FEMALE_CHROMS),
                MIN_HET_FRACTION,
                MIN_GQ,
                MIN_DP,
                generateTrioMetricsBase(),
                ImmutableList.copyOf(parIntervals.get()),
                progressLogger
        );
    }

    /** Generates a {@link MendelianViolationMetrics} objects for each trio. */
    private List<MendelianViolationMetrics> generateTrioMetricsBase() {
        final List<MendelianViolationMetrics> metrics = new ArrayList<>();

        for (final PedFile.PedTrio trio : pedFile.get().values()) {
            final MendelianViolationMetrics m = new MendelianViolationMetrics();
            m.MOTHER = trio.getMaternalId();
            m.FATHER = trio.getPaternalId();
            m.OFFSPRING = trio.getIndividualId();
            m.FAMILY_ID = trio.getFamilyId();
            m.OFFSPRING_SEX = trio.getSex();
            metrics.add(m);
        }

        ///////////////////////////////////////////////////////////////////////
        // Filter out trios where one or more of the samples is not in the VCF
        ///////////////////////////////////////////////////////////////////////
        final Set<String> allSamples = new HashSet<>(inputHeader.get().getSampleNamesInOrder());
        final Iterator<MendelianViolationMetrics> trioIterator = metrics.iterator();

        while (trioIterator.hasNext()) {
            final MendelianViolationMetrics m = trioIterator.next();
            final Set<String> trio = CollectionUtil.makeSet(m.MOTHER, m.FATHER, m.OFFSPRING);
            trio.removeAll(allSamples);

            if (!trio.isEmpty()) {
                LOG.warn("Removing trio due to the following missing samples in VCF: " + trio);
                trioIterator.remove();
            }
        }
        return metrics;
    }
}
