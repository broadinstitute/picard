/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.variantcontext;

import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReaderUtil;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.*;
import org.testng.Assert;

import java.io.*;
import java.util.*;

/**
 * Routines for generating all sorts of VCs for testing
 *
 * @author Your Name
 * @since Date created
 */
public class VariantContextTestProvider {
    final private static boolean ENABLE_GENOTYPE_TESTS = true;
    final private static boolean ENABLE_A_AND_G_TESTS = true;
    final private static boolean ENABLE_VARARRAY_TESTS = true;
    final private static boolean ENABLE_PLOIDY_TESTS = true;
    final private static boolean ENABLE_PL_TESTS = true;
    final private static boolean ENABLE_SYMBOLIC_ALLELE_TESTS = true;
    final private static boolean ENABLE_SOURCE_VCF_TESTS = true;
    final private static boolean ENABLE_VARIABLE_LENGTH_GENOTYPE_STRING_TESTS = true;
    final private static List<Integer> TWENTY_INTS = Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20);

    private static VCFHeader syntheticHeader;
    final static List<VariantContextTestData> TEST_DATAs = new ArrayList<VariantContextTestData>();
    private static VariantContext ROOT;

    private final static List<File> testSourceVCFs = new ArrayList<File>();
    static {
        testSourceVCFs.add(new File(VariantBaseTest.variantTestDataRoot + "ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf"));
        testSourceVCFs.add(new File(VariantBaseTest.variantTestDataRoot + "ex2.vcf"));
        testSourceVCFs.add(new File(VariantBaseTest.variantTestDataRoot + "dbsnp_135.b37.1000.vcf"));
        if ( ENABLE_SYMBOLIC_ALLELE_TESTS ) {
            testSourceVCFs.add(new File(VariantBaseTest.variantTestDataRoot + "diagnosis_targets_testfile.vcf"));
            testSourceVCFs.add(new File(VariantBaseTest.variantTestDataRoot + "VQSR.mixedTest.recal"));
        }
    }

    public static class VariantContextContainer {
        private VCFHeader header;
        private Iterable<VariantContext> vcs;

        public VariantContextContainer( VCFHeader header, Iterable<VariantContext> vcs ) {
            this.header = header;
            this.vcs = vcs;
        }

        public VCFHeader getHeader() {
            return header;
        }

        public Iterable<VariantContext> getVCs() {
            return vcs;
        }
    }

    public abstract static class VariantContextIOTest<CODECTYPE> {
        public String toString() {
            return "VariantContextIOTest:" + getExtension();
        }
        public abstract String getExtension();
        public abstract CODECTYPE makeCodec();
        public abstract VariantContextWriter makeWriter(final File outputFile, final EnumSet<Options> baseOptions);

        public abstract VariantContextContainer readAllVCs(final File input) throws IOException;
        
        public List<VariantContext> preprocess(final VCFHeader header, List<VariantContext> vcsBeforeIO) {
            return vcsBeforeIO;
        }

        public List<VariantContext> postprocess(final VCFHeader header, List<VariantContext> vcsAfterIO) {
            return vcsAfterIO;
        }
    }

    public static class VariantContextTestData {
        public final VCFHeader header;
        public List<VariantContext> vcs;

        public VariantContextTestData(final VCFHeader header, final VariantContextBuilder builder) {
            this(header, Collections.singletonList(builder.fullyDecoded(true).make()));
        }

        public VariantContextTestData(final VCFHeader header, final List<VariantContext> vcs) {
            final Set<String> samples = new HashSet<String>();
            for ( final VariantContext vc : vcs )
                if ( vc.hasGenotypes() )
                    samples.addAll(vc.getSampleNames());
            this.header = samples.isEmpty() ? header : new VCFHeader(header.getMetaDataInSortedOrder(), samples);
            this.vcs = vcs;
        }

        public boolean hasGenotypes() {
            return vcs.get(0).hasGenotypes();
        }

        public String toString() {
            StringBuilder b = new StringBuilder();
            b.append("VariantContextTestData: [");
            final VariantContext vc = vcs.get(0);
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            builder.noGenotypes();
            b.append(builder.make().toString());
            if ( vc.getNSamples() < 5 ) {
                for ( final Genotype g : vc.getGenotypes() )
                    b.append(g.toString());
            } else {
                b.append(" nGenotypes = ").append(vc.getNSamples());
            }

            if ( vcs.size() > 1 ) b.append(" ----- with another ").append(vcs.size() - 1).append(" VariantContext records");
            b.append("]");
            return b.toString();
        }
    }

    private final static VariantContextBuilder builder() {
        return new VariantContextBuilder(ROOT);
    }

    private final static void add(VariantContextBuilder builder) {
        TEST_DATAs.add(new VariantContextTestData(syntheticHeader, builder));
    }

    public static void initializeTests() throws IOException {
        createSyntheticHeader();
        makeSyntheticTests();
        makeEmpiricalTests();
    }

    private static void makeEmpiricalTests() throws IOException {
        if ( ENABLE_SOURCE_VCF_TESTS ) {
            for ( final File file : testSourceVCFs ) {
                VCFCodec codec = new VCFCodec();
                VariantContextContainer x = readAllVCs( file, codec );
                List<VariantContext> fullyDecoded = new ArrayList<VariantContext>();

                for ( final VariantContext raw : x.getVCs() ) {
                    if ( raw != null )
                        fullyDecoded.add(raw.fullyDecode(x.getHeader(), false));
                }

                TEST_DATAs.add(new VariantContextTestData(x.getHeader(), fullyDecoded));
            }
        }
    }

    private final static void addHeaderLine(final Set<VCFHeaderLine> metaData, final String id, final int count, final VCFHeaderLineType type) {
        metaData.add(new VCFInfoHeaderLine(id, count, type, "x"));
        if ( type != VCFHeaderLineType.Flag )
            metaData.add(new VCFFormatHeaderLine(id, count, type, "x"));
    }

    private final static void addHeaderLine(final Set<VCFHeaderLine> metaData, final String id, final VCFHeaderLineCount count, final VCFHeaderLineType type) {
        metaData.add(new VCFInfoHeaderLine(id, count, type, "x"));
        if ( type != VCFHeaderLineType.Flag )
            metaData.add(new VCFFormatHeaderLine(id, count, type, "x"));
    }

    private static void createSyntheticHeader() {
        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();

        addHeaderLine(metaData, "STRING1", 1, VCFHeaderLineType.String);
        addHeaderLine(metaData, "END", 1, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "STRING3", 3, VCFHeaderLineType.String);
        addHeaderLine(metaData, "STRING20", 20, VCFHeaderLineType.String);
        addHeaderLine(metaData, "VAR.INFO.STRING", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String);

        addHeaderLine(metaData, "GT", 1, VCFHeaderLineType.String);
        addHeaderLine(metaData, "GQ", 1, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "ADA", VCFHeaderLineCount.A, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "PL", VCFHeaderLineCount.G, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "GS", 2, VCFHeaderLineType.String);
        addHeaderLine(metaData, "GV", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String);
        addHeaderLine(metaData, "FT", 1, VCFHeaderLineType.String);

        // prep the header
        metaData.add(new VCFContigHeaderLine(Collections.singletonMap("ID", "1"), 0));

        metaData.add(new VCFFilterHeaderLine("FILTER1"));
        metaData.add(new VCFFilterHeaderLine("FILTER2"));

        addHeaderLine(metaData, "INT1", 1, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "INT3", 3, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "INT20", 20, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "INT.VAR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "FLOAT1", 1, VCFHeaderLineType.Float);
        addHeaderLine(metaData, "FLOAT3", 3, VCFHeaderLineType.Float);
        addHeaderLine(metaData, "FLAG", 0, VCFHeaderLineType.Flag);

        syntheticHeader = new VCFHeader(metaData);
    }


    private static void makeSyntheticTests() {
        VariantContextBuilder rootBuilder = new VariantContextBuilder();
        rootBuilder.source("test");
        rootBuilder.loc("1", 10, 10);
        rootBuilder.alleles("A", "C");
        rootBuilder.unfiltered();
        ROOT = rootBuilder.make();

        add(builder());
        add(builder().alleles("A"));
        add(builder().alleles("A", "C", "T"));
        add(builder().alleles("A", "AC"));
        add(builder().alleles("A", "ACAGT"));
        add(builder().loc("1", 10, 11).alleles("AC", "A"));
        add(builder().loc("1", 10, 13).alleles("ACGT", "A"));

        // make sure filters work
        add(builder().unfiltered());
        add(builder().passFilters());
        add(builder().filters("FILTER1"));
        add(builder().filters("FILTER1", "FILTER2"));

        add(builder().log10PError(VariantContext.NO_LOG10_PERROR));
        add(builder().log10PError(-1));
        add(builder().log10PError(-1.234e6));

        add(builder().noID());
        add(builder().id("rsID12345"));


        add(builder().attribute("INT1", 1));
        add(builder().attribute("INT1", 100));
        add(builder().attribute("INT1", 1000));
        add(builder().attribute("INT1", 100000));
        add(builder().attribute("INT1", null));
        add(builder().attribute("INT3", Arrays.asList(1, 2, 3)));
        add(builder().attribute("INT3", Arrays.asList(1000, 2000, 3000)));
        add(builder().attribute("INT3", Arrays.asList(100000, 200000, 300000)));
        add(builder().attribute("INT3", null));
        add(builder().attribute("INT20", TWENTY_INTS));

        add(builder().attribute("FLOAT1", 1.0));
        add(builder().attribute("FLOAT1", 100.0));
        add(builder().attribute("FLOAT1", 1000.0));
        add(builder().attribute("FLOAT1", 100000.0));
        add(builder().attribute("FLOAT1", null));
        add(builder().attribute("FLOAT3", Arrays.asList(1.0, 2.0, 3.0)));
        add(builder().attribute("FLOAT3", Arrays.asList(1000.0, 2000.0, 3000.0)));
        add(builder().attribute("FLOAT3", Arrays.asList(100000.0, 200000.0, 300000.0)));
        add(builder().attribute("FLOAT3", null));

        add(builder().attribute("FLAG", true));
        //add(builder().attribute("FLAG", false)); // NOTE -- VCF doesn't allow false flags

        add(builder().attribute("STRING1", "s1"));
        add(builder().attribute("STRING1", null));
        add(builder().attribute("STRING3", Arrays.asList("s1", "s2", "s3")));
        add(builder().attribute("STRING3", null));
        add(builder().attribute("STRING20", Arrays.asList("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20")));

        add(builder().attribute("VAR.INFO.STRING", "s1"));
        add(builder().attribute("VAR.INFO.STRING", Arrays.asList("s1", "s2")));
        add(builder().attribute("VAR.INFO.STRING", Arrays.asList("s1", "s2", "s3")));
        add(builder().attribute("VAR.INFO.STRING", null));

        if ( ENABLE_GENOTYPE_TESTS ) {
            addGenotypesToTestData();
            addComplexGenotypesTest();
        }

        if ( ENABLE_A_AND_G_TESTS )
            addGenotypesAndGTests();

        if ( ENABLE_SYMBOLIC_ALLELE_TESTS )
            addSymbolicAlleleTests();
    }

    private static void addSymbolicAlleleTests() {
        // two tests to ensure that the end is computed correctly when there's (and not) an END field present
        add(builder().alleles("N", "<VQSR>").start(10).stop(11).attribute("END", 11));
        add(builder().alleles("N", "<VQSR>").start(10).stop(10));
    }

    private static void addGenotypesToTestData() {
        final ArrayList<VariantContext> sites = new ArrayList<VariantContext>();

        sites.add(builder().alleles("A").make());
        sites.add(builder().alleles("A", "C", "T").make());
        sites.add(builder().alleles("A", "AC").make());
        sites.add(builder().alleles("A", "ACAGT").make());

        for ( VariantContext site : sites ) {
            addGenotypes(site);
        }
    }

    private static void addGenotypeTests( final VariantContext site, Genotype ... genotypes ) {
        // for each sites VC, we are going to add create two root genotypes.
        // The first is the primary, and will be added to each new test
        // The second is variable.  In some tests it's absent (testing 1 genotype), in others it is duplicated
        // 1 once, 10, 100, or 1000 times to test scaling

        final VariantContextBuilder builder = new VariantContextBuilder(site);

        // add a single context
        builder.genotypes(genotypes[0]);
        add(builder);

        if ( genotypes.length > 1 ) {
            // add all
            add(builder.genotypes(Arrays.asList(genotypes)));

            // add all with the last replicated 10x and 100x times
            for ( int nCopiesOfLast : Arrays.asList(10, 100, 1000) ) {
                final GenotypesContext gc = new GenotypesContext();
                final Genotype last = genotypes[genotypes.length-1];
                for ( int i = 0; i < genotypes.length - 1; i++ )
                    gc.add(genotypes[i]);
                for ( int i = 0; i < nCopiesOfLast; i++ )
                    gc.add(new GenotypeBuilder(last).name("copy" + i).make());
                add(builder.genotypes(gc));
            }
        }
    }

    private static void addGenotypes( final VariantContext site) {
        // test ref/ref
        final Allele ref = site.getReference();
        final Allele alt1 = site.getNAlleles() > 1 ? site.getAlternateAllele(0) : null;
        final Genotype homRef = GenotypeBuilder.create("homRef", Arrays.asList(ref, ref));
        addGenotypeTests(site, homRef);

        if ( alt1 != null ) {
            final Genotype het = GenotypeBuilder.create("het", Arrays.asList(ref, alt1));
            final Genotype homVar = GenotypeBuilder.create("homVar", Arrays.asList(alt1, alt1));
            addGenotypeTests(site, homRef, het);
            addGenotypeTests(site, homRef, het, homVar);

            // test no GT at all
            addGenotypeTests(site, new GenotypeBuilder("noGT", new ArrayList<Allele>(0)).attribute("INT1", 10).make());

            final List<Allele> noCall = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

            // ploidy
            if ( ENABLE_PLOIDY_TESTS ) {
                addGenotypeTests(site,
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall",  noCall),
                        GenotypeBuilder.create("noCall2", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall", noCall),
                        GenotypeBuilder.create("noCall2", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("nocall", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));
            }


            //
            //
            // TESTING PHASE
            //
            //
            final Genotype gUnphased = new GenotypeBuilder("gUnphased", Arrays.asList(ref, alt1)).make();
            final Genotype gPhased = new GenotypeBuilder("gPhased", Arrays.asList(ref, alt1)).phased(true).make();
            final Genotype gPhased2 = new GenotypeBuilder("gPhased2", Arrays.asList(alt1, alt1)).phased(true).make();
            final Genotype gPhased3 = new GenotypeBuilder("gPhased3", Arrays.asList(ref, ref)).phased(true).make();
            final Genotype haploidNoPhase = new GenotypeBuilder("haploidNoPhase", Arrays.asList(ref)).make();
            addGenotypeTests(site, gUnphased, gPhased);
            addGenotypeTests(site, gUnphased, gPhased2);
            addGenotypeTests(site, gUnphased, gPhased3);
            addGenotypeTests(site, gPhased, gPhased2);
            addGenotypeTests(site, gPhased, gPhased3);
            addGenotypeTests(site, gPhased2, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased);
            addGenotypeTests(site, haploidNoPhase, gPhased2);
            addGenotypeTests(site, haploidNoPhase, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased, gPhased2);
            addGenotypeTests(site, haploidNoPhase, gPhased, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased2, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased, gPhased2, gPhased3);

            final Genotype gUnphasedTet = new GenotypeBuilder("gUnphasedTet", Arrays.asList(ref, alt1, ref, alt1)).make();
            final Genotype gPhasedTet = new GenotypeBuilder("gPhasedTet", Arrays.asList(ref, alt1, alt1, alt1)).phased(true).make();
            addGenotypeTests(site, gUnphasedTet, gPhasedTet);
        }

        if ( ENABLE_PL_TESTS ) {
            if ( site.getNAlleles() == 2 ) {
                // testing PLs
                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{0, -1, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2, -3}));

                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{-1, 0, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2, -3}));

                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{-1, 0, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2000, -1000}));

                addGenotypeTests(site, // missing PLs
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{-1, 0, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref)));
            }
            else if ( site.getNAlleles() == 3 ) {
                // testing PLs
                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{0, -1, -2, -3, -4, -5}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2, -3, -4, -5, -6}));
            }
        }

        // test attributes
        addGenotypeTests(site,
                attr("g1", ref, "INT1", 1),
                attr("g2", ref, "INT1", 2));
        addGenotypeTests(site,
                attr("g1", ref, "INT1", 1),
                attr("g2", ref, "INT1"));
        addGenotypeTests(site,
                attr("g1", ref, "INT3", 1, 2, 3),
                attr("g2", ref, "INT3", 4, 5, 6));
        addGenotypeTests(site,
                attr("g1", ref, "INT3", 1, 2, 3),
                attr("g2", ref, "INT3"));

        addGenotypeTests(site,
                attr("g1", ref, "INT20", TWENTY_INTS),
                attr("g2", ref, "INT20", TWENTY_INTS));


        if (ENABLE_VARARRAY_TESTS) {
            addGenotypeTests(site,
                    attr("g1", ref, "INT.VAR", 1, 2, 3),
                    attr("g2", ref, "INT.VAR", 4, 5),
                    attr("g3", ref, "INT.VAR", 6));
            addGenotypeTests(site,
                    attr("g1", ref, "INT.VAR", 1, 2, 3),
                    attr("g2", ref, "INT.VAR"),
                    attr("g3", ref, "INT.VAR", 5));
        }

        addGenotypeTests(site,
                attr("g1", ref, "FLOAT1", 1.0),
                attr("g2", ref, "FLOAT1", 2.0));
        addGenotypeTests(site,
                attr("g1", ref, "FLOAT1", 1.0),
                attr("g2", ref, "FLOAT1"));
        addGenotypeTests(site,
                attr("g1", ref, "FLOAT3", 1.0, 2.0, 3.0),
                attr("g2", ref, "FLOAT3", 4.0, 5.0, 6.0));
        addGenotypeTests(site,
                attr("g1", ref, "FLOAT3", 1.0, 2.0, 3.0),
                attr("g2", ref, "FLOAT3"));

        if (ENABLE_VARIABLE_LENGTH_GENOTYPE_STRING_TESTS) {
            //
            //
            // TESTING MULTIPLE SIZED LISTS IN THE GENOTYPE FIELD
            //
            //
            addGenotypeTests(site,
                    attr("g1", ref, "GS", Arrays.asList("S1", "S2")),
                    attr("g2", ref, "GS", Arrays.asList("S3", "S4")));

            addGenotypeTests(site, // g1 is missing the string, and g2 is missing FLOAT1
                    attr("g1", ref, "FLOAT1", 1.0),
                    attr("g2", ref, "GS", Arrays.asList("S3", "S4")));

            // variable sized lists
            addGenotypeTests(site,
                    attr("g1", ref, "GV", "S1"),
                    attr("g2", ref, "GV", Arrays.asList("S3", "S4")));

            addGenotypeTests(site,
                    attr("g1", ref, "GV", Arrays.asList("S1", "S2")),
                    attr("g2", ref, "GV", Arrays.asList("S3", "S4", "S5")));

            addGenotypeTests(site, // missing value in varlist of string
                    attr("g1", ref, "FLOAT1", 1.0),
                    attr("g2", ref, "GV", Arrays.asList("S3", "S4", "S5")));
        }

        //
        //
        // TESTING GENOTYPE FILTERS
        //
        //
        addGenotypeTests(site,
                new GenotypeBuilder("g1-x", Arrays.asList(ref, ref)).filters("X").make(),
                new GenotypeBuilder("g2-x", Arrays.asList(ref, ref)).filters("X").make());
        addGenotypeTests(site,
                new GenotypeBuilder("g1-unft", Arrays.asList(ref, ref)).unfiltered().make(),
                new GenotypeBuilder("g2-x", Arrays.asList(ref, ref)).filters("X").make());
        addGenotypeTests(site,
                new GenotypeBuilder("g1-unft", Arrays.asList(ref, ref)).unfiltered().make(),
                new GenotypeBuilder("g2-xy", Arrays.asList(ref, ref)).filters("X", "Y").make());
        addGenotypeTests(site,
                new GenotypeBuilder("g1-unft", Arrays.asList(ref, ref)).unfiltered().make(),
                new GenotypeBuilder("g2-x", Arrays.asList(ref, ref)).filters("X").make(),
                new GenotypeBuilder("g3-xy", Arrays.asList(ref, ref)).filters("X", "Y").make());
    }

    private static void addGenotypesAndGTests() {
//        for ( final int ploidy : Arrays.asList(2)) {
        for ( final int ploidy : Arrays.asList(1, 2, 3, 4, 5)) {
            final List<List<String>> alleleCombinations =
                    Arrays.asList(
                            Arrays.asList("A"),
                            Arrays.asList("A", "C"),
                            Arrays.asList("A", "C", "G"),
                            Arrays.asList("A", "C", "G", "T"));

            for ( final List<String> alleles : alleleCombinations ) {
                final VariantContextBuilder vcb = builder().alleles(alleles);
                final VariantContext site = vcb.make();
                final int nAlleles = site.getNAlleles();
                final Allele ref = site.getReference();

                // base genotype is ref/.../ref up to ploidy
                final List<Allele> baseGenotype = new ArrayList<Allele>(ploidy);
                for ( int i = 0; i < ploidy; i++) baseGenotype.add(ref);
                final int nPLs = GenotypeLikelihoods.numLikelihoods(nAlleles, ploidy);

                // ada is 0, 1, ..., nAlleles - 1
                final List<Integer> ada = new ArrayList<Integer>(nAlleles);
                for ( int i = 0; i < nAlleles - 1; i++ ) ada.add(i);

                // pl is 0, 1, ..., up to nPLs (complex calc of nAlleles and ploidy)
                final int[] pl = new int[nPLs];
                for ( int i = 0; i < pl.length; i++ ) pl[i] = i;

                final GenotypeBuilder gb = new GenotypeBuilder("ADA_PL_SAMPLE");
                gb.alleles(baseGenotype);
                gb.PL(pl);
                gb.attribute("ADA", nAlleles == 2 ? ada.get(0) : ada);
                vcb.genotypes(gb.make());

                add(vcb);
            }
        }
    }

    private static Genotype attr(final String name, final Allele ref, final String key, final Object ... value) {
        if ( value.length == 0 )
            return GenotypeBuilder.create(name, Arrays.asList(ref, ref));
        else {
            final Object toAdd = value.length == 1 ? value[0] : Arrays.asList(value);
            return new GenotypeBuilder(name, Arrays.asList(ref, ref)).attribute(key, toAdd).make();
        }
    }

    public static List<VariantContextTestData> generateSiteTests() {
        return TEST_DATAs;
    }

    public static void testReaderWriterWithMissingGenotypes(final VariantContextIOTest tester, final VariantContextTestData data) throws IOException {
        final int nSamples = data.header.getNGenotypeSamples();
        if ( nSamples > 2 ) {
            for ( final VariantContext vc : data.vcs )
                if ( vc.isSymbolic() )
                    // cannot handle symbolic alleles because they may be weird non-call VCFs
                    return;

            final File tmpFile = File.createTempFile("testReaderWriter", tester.getExtension());
            tmpFile.deleteOnExit();

            // write expected to disk
            final EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
            final VariantContextWriter writer = tester.makeWriter(tmpFile, options);

            final Set<String> samplesInVCF = new HashSet<String>(data.header.getGenotypeSamples());
            final List<String> missingSamples = Arrays.asList("MISSING1", "MISSING2");
            final List<String> allSamples = new ArrayList<String>(missingSamples);
            allSamples.addAll(samplesInVCF);

            final VCFHeader header = new VCFHeader(data.header.getMetaDataInInputOrder(), allSamples);
            writeVCsToFile(writer, header, data.vcs);

            // ensure writing of expected == actual
            final VariantContextContainer p = tester.readAllVCs(tmpFile);
            final Iterable<VariantContext> actual = p.getVCs();

            int i = 0;
            for ( final VariantContext readVC : actual ) {
                if ( readVC == null ) continue; // sometimes we read null records...
                final VariantContext expected = data.vcs.get(i++);
                for ( final Genotype g : readVC.getGenotypes() ) {
                    Assert.assertTrue(allSamples.contains(g.getSampleName()));
                    if ( samplesInVCF.contains(g.getSampleName()) ) {
                        assertEquals(g, expected.getGenotype(g.getSampleName()));
                    } else {
                        // missing
                        Assert.assertTrue(g.isNoCall());
                    }
                }
            }

        }
    }

    public static void testReaderWriter(final VariantContextIOTest tester, final VariantContextTestData data) throws IOException {
        testReaderWriter(tester, data.header, data.vcs, data.vcs, true);
    }

    public static void testReaderWriter(final VariantContextIOTest tester,
                                        final VCFHeader header,
                                        final List<VariantContext> expected,
                                        final Iterable<VariantContext> vcs,
                                        final boolean recurse) throws IOException {
        final File tmpFile = File.createTempFile("testReaderWriter", tester.getExtension());
        tmpFile.deleteOnExit();

        // write expected to disk
        final EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter writer = tester.makeWriter(tmpFile, options);
        writeVCsToFile(writer, header, vcs);

        // ensure writing of expected == actual
        final VariantContextContainer p = tester.readAllVCs(tmpFile);
        final Iterable<VariantContext> actual = p.getVCs();
        assertEquals(actual, expected);

        if ( recurse ) {
            // if we are doing a recursive test, grab a fresh iterator over the written values
            final Iterable<VariantContext> read = tester.readAllVCs(tmpFile).getVCs();
            testReaderWriter(tester, p.getHeader(), expected, read, false);
        }
    }

    private static void writeVCsToFile(final VariantContextWriter writer, final VCFHeader header, final Iterable<VariantContext> vcs) {
        // write
        writer.writeHeader(header);
        for ( VariantContext vc : vcs )
            if (vc != null)
                writer.add(vc);
        writer.close();
    }

    public static abstract class VCIterable<SOURCE> implements Iterable<VariantContext>, Iterator<VariantContext> {
        final FeatureCodec<VariantContext, SOURCE> codec;
        final VCFHeader header;

        public VCIterable(final FeatureCodec<VariantContext, SOURCE> codec, final VCFHeader header) {
            this.codec = codec;
            this.header = header;
        }

        @Override
        public Iterator<VariantContext> iterator() {
            return this;
        }

        @Override
        public abstract boolean hasNext();

        public abstract SOURCE nextSource();
        
        @Override
        public VariantContext next() {
            try {
                final VariantContext vc = codec.decode(nextSource());
                return vc == null ? null : vc.fullyDecode(header, false);
            } catch ( IOException e ) {
                throw new RuntimeException(e);
            }
        }

        @Override
        public void remove() { }
    }

    public static VariantContextContainer readAllVCs(final File input, final BCF2Codec codec) throws IOException {
        PositionalBufferedStream headerPbs = new PositionalBufferedStream(new FileInputStream(input));
        FeatureCodecHeader header = codec.readHeader(headerPbs);
        headerPbs.close();

        final PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(input));
        pbs.skip(header.getHeaderEnd());

        final VCFHeader vcfHeader = (VCFHeader)header.getHeaderValue();
        return new VariantContextTestProvider.VariantContextContainer(vcfHeader, new VariantContextTestProvider.VCIterable(codec, vcfHeader) {
            @Override
            public boolean hasNext() {
                try {
                    return !pbs.isDone();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }

            @Override
            public Object nextSource() {
                return pbs;
            }
        });
    }

    public static VariantContextContainer readAllVCs(final File input, final VCFCodec codec) throws FileNotFoundException {
        final LineIterator lineIterator = new LineIteratorImpl(LineReaderUtil.fromBufferedStream(new BufferedInputStream(new FileInputStream(input))));
        final VCFHeader vcfHeader = (VCFHeader) codec.readActualHeader(lineIterator);
        return new VariantContextTestProvider.VariantContextContainer(vcfHeader, new VariantContextTestProvider.VCIterable<LineIterator>(codec, vcfHeader) {
            @Override
            public boolean hasNext() {
                return lineIterator.hasNext();
            }

            @Override
            public LineIterator nextSource() {
                return lineIterator;
            }
        });
    }
    
    public static void assertVCFandBCFFilesAreTheSame(final File vcfFile, final File bcfFile) throws IOException {
        final VariantContextContainer vcfData = readAllVCs(vcfFile, new VCFCodec());
        final VariantContextContainer bcfData = readAllVCs(bcfFile, new BCF2Codec());
        assertEquals(bcfData.getHeader(), vcfData.getHeader());
        assertEquals(bcfData.getVCs(), vcfData.getVCs());
    }

    public static void assertEquals(final Iterable<VariantContext> actual, final Iterable<VariantContext> expected) {
        final Iterator<VariantContext> actualIT = actual.iterator();
        final Iterator<VariantContext> expectedIT = expected.iterator();

        while ( expectedIT.hasNext() ) {
            final VariantContext expectedVC = expectedIT.next();
            if ( expectedVC == null )
                continue;

            VariantContext actualVC;
            do {
                Assert.assertTrue(actualIT.hasNext(), "Too few records found in actual");
                actualVC = actualIT.next();
            } while ( actualIT.hasNext() && actualVC == null );

            if ( actualVC == null )
                Assert.fail("Too few records in actual");

            assertEquals(actualVC, expectedVC);
        }
        Assert.assertTrue(! actualIT.hasNext(), "Too many records found in actual");
    }

    /**
     * Assert that two variant contexts are actually equal
     * @param actual
     * @param expected
     */
    public static void assertEquals( final VariantContext actual, final VariantContext expected ) {
        Assert.assertNotNull(actual, "VariantContext expected not null");
        Assert.assertEquals(actual.getChr(), expected.getChr(), "chr");
        Assert.assertEquals(actual.getStart(), expected.getStart(), "start");
        Assert.assertEquals(actual.getEnd(), expected.getEnd(), "end");
        Assert.assertEquals(actual.getID(), expected.getID(), "id");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "alleles for " + expected + " vs " + actual);

        assertAttributesEquals(actual.getAttributes(), expected.getAttributes());
        Assert.assertEquals(actual.filtersWereApplied(), expected.filtersWereApplied(), "filtersWereApplied");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "isFiltered");
        VariantBaseTest.assertEqualsSet(actual.getFilters(), expected.getFilters(), "filters");
        VariantBaseTest.assertEqualsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual());

        Assert.assertEquals(actual.hasGenotypes(), expected.hasGenotypes(), "hasGenotypes");
        if ( expected.hasGenotypes() ) {
            VariantBaseTest.assertEqualsSet(actual.getSampleNames(), expected.getSampleNames(), "sample names set");
            Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "sample names");
            final Set<String> samples = expected.getSampleNames();
            for ( final String sample : samples ) {
                assertEquals(actual.getGenotype(sample), expected.getGenotype(sample));
            }
        }
    }

    public static void assertEquals(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype names");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "Genotype alleles");
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString(), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");

        // filters are the same
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "Genotype fields");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "Genotype isFiltered");

        // inline attributes
        Assert.assertEquals(actual.getDP(), expected.getDP(), "Genotype dp");
        Assert.assertTrue(Arrays.equals(actual.getAD(), expected.getAD()));
        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype gq");
        Assert.assertEquals(actual.hasPL(), expected.hasPL(), "Genotype hasPL");
        Assert.assertEquals(actual.hasAD(), expected.hasAD(), "Genotype hasAD");
        Assert.assertEquals(actual.hasGQ(), expected.hasGQ(), "Genotype hasGQ");
        Assert.assertEquals(actual.hasDP(), expected.hasDP(), "Genotype hasDP");

        Assert.assertEquals(actual.hasLikelihoods(), expected.hasLikelihoods(), "Genotype haslikelihoods");
        Assert.assertEquals(actual.getLikelihoodsString(), expected.getLikelihoodsString(), "Genotype getlikelihoodsString");
        Assert.assertEquals(actual.getLikelihoods(), expected.getLikelihoods(), "Genotype getLikelihoods");
        Assert.assertTrue(Arrays.equals(actual.getPL(), expected.getPL()));

        Assert.assertEquals(actual.getPhredScaledQual(), expected.getPhredScaledQual(), "Genotype phredScaledQual");
        assertAttributesEquals(actual.getExtendedAttributes(), expected.getExtendedAttributes());
        Assert.assertEquals(actual.isPhased(), expected.isPhased(), "Genotype isPhased");
        Assert.assertEquals(actual.getPloidy(), expected.getPloidy(), "Genotype getPloidy");
    }

    private static void assertAttributesEquals(final Map<String, Object> actual, Map<String, Object> expected) {
        final Set<String> expectedKeys = new HashSet<String>(expected.keySet());

        for ( final Map.Entry<String, Object> act : actual.entrySet() ) {
            final Object actualValue = act.getValue();
            if ( expected.containsKey(act.getKey()) && expected.get(act.getKey()) != null ) {
                final Object expectedValue = expected.get(act.getKey());
                if ( expectedValue instanceof List ) {
                    final List<Object> expectedList = (List<Object>)expectedValue;
                    Assert.assertTrue(actualValue instanceof List, act.getKey() + " should be a list but isn't");
                    final List<Object> actualList = (List<Object>)actualValue;
                    Assert.assertEquals(actualList.size(), expectedList.size(), act.getKey() + " size");
                    for ( int i = 0; i < expectedList.size(); i++ )
                        assertAttributeEquals(act.getKey(), actualList.get(i), expectedList.get(i));
                } else
                    assertAttributeEquals(act.getKey(), actualValue, expectedValue);
            } else {
                // it's ok to have a binding in x -> null that's absent in y
                Assert.assertNull(actualValue, act.getKey() + " present in one but not in the other");
            }
            expectedKeys.remove(act.getKey());
        }

        // now expectedKeys contains only the keys found in expected but not in actual,
        // and they must all be null
        for ( final String missingExpected : expectedKeys ) {
            final Object value = expected.get(missingExpected);
            Assert.assertTrue(isMissing(value), "Attribute " + missingExpected + " missing in one but not in other" );
        }
    }

    private static final boolean isMissing(final Object value) {
        if ( value == null ) return true;
        else if ( value.equals(VCFConstants.MISSING_VALUE_v4) ) return true;
        else if ( value instanceof List ) {
            // handles the case where all elements are null or the list is empty
            for ( final Object elt : (List)value)
                if ( elt != null )
                    return false;
            return true;
        } else
            return false;
    }

    private static void assertAttributeEquals(final String key, final Object actual, final Object expected) {
        if ( expected instanceof Double ) {
            // must be very tolerant because doubles are being rounded to 2 sig figs
            VariantBaseTest.assertEqualsDoubleSmart(actual, (Double)expected, 1e-2);
        } else
            Assert.assertEquals(actual, expected, "Attribute " + key);
    }

    public static void addComplexGenotypesTest() {
        final List<Allele> allAlleles = Arrays.asList(
                Allele.create("A", true),
                Allele.create("C", false),
                Allele.create("G", false));

        for ( int nAlleles : Arrays.asList(2, 3) ) {
            for ( int highestPloidy : Arrays.asList(1, 2, 3) ) {
                // site alleles
                final List<Allele> siteAlleles = allAlleles.subList(0, nAlleles);

                // possible alleles for genotypes
                final List<Allele> possibleGenotypeAlleles = new ArrayList<Allele>(siteAlleles);
                possibleGenotypeAlleles.add(Allele.NO_CALL);

                // there are n^ploidy possible genotypes
                final List<List<Allele>> possibleGenotypes = makeAllGenotypes(possibleGenotypeAlleles, highestPloidy);
                final int nPossibleGenotypes = possibleGenotypes.size();

                VariantContextBuilder vb = new VariantContextBuilder("unittest", "1", 1, 1, siteAlleles);

                // first test -- create n copies of each genotype
                for ( int i = 0; i < nPossibleGenotypes; i++ ) {
                    final List<Genotype> samples = new ArrayList<Genotype>(3);
                    samples.add(GenotypeBuilder.create("sample" + i, possibleGenotypes.get(i)));
                    add(vb.genotypes(samples));
                }

                // second test -- create one sample with each genotype
                {
                    final List<Genotype> samples = new ArrayList<Genotype>(nPossibleGenotypes);
                    for ( int i = 0; i < nPossibleGenotypes; i++ ) {
                        samples.add(GenotypeBuilder.create("sample" + i, possibleGenotypes.get(i)));
                    }
                    add(vb.genotypes(samples));
                }

                // test mixed ploidy
                for ( int i = 0; i < nPossibleGenotypes; i++ ) {
                    for ( int ploidy = 1; ploidy < highestPloidy; ploidy++ ) {
                        final List<Genotype> samples = new ArrayList<Genotype>(highestPloidy);
                        final List<Allele> genotype = possibleGenotypes.get(i).subList(0, ploidy);
                        samples.add(GenotypeBuilder.create("sample" + i, genotype));
                        add(vb.genotypes(samples));
                    }
                }
            }
        }
    }

    private static List<List<Allele>> makeAllGenotypes(final List<Allele> alleles, final int highestPloidy) {
        return GeneralUtils.makePermutations(alleles, highestPloidy, true);
    }

    public static void assertEquals(final VCFHeader actual, final VCFHeader expected) {
        Assert.assertEquals(actual.getMetaDataInSortedOrder().size(), expected.getMetaDataInSortedOrder().size(), "No VCF header lines");

        // for some reason set.equals() is returning false but all paired elements are .equals().  Perhaps compare to is busted?
        //Assert.assertEquals(actual.getMetaDataInInputOrder(), expected.getMetaDataInInputOrder());
        final List<VCFHeaderLine> actualLines = new ArrayList<VCFHeaderLine>(actual.getMetaDataInSortedOrder());
        final List<VCFHeaderLine> expectedLines = new ArrayList<VCFHeaderLine>(expected.getMetaDataInSortedOrder());
        for ( int i = 0; i < actualLines.size(); i++ ) {
            Assert.assertEquals(actualLines.get(i), expectedLines.get(i), "VCF header lines");
        }
    }

    public static void main( String argv[] ) {
        final File variants1 = new File(argv[0]);
        final File variants2 = new File(argv[1]);
        try {
            VariantContextTestProvider.assertVCFandBCFFilesAreTheSame(variants1, variants2);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }
}