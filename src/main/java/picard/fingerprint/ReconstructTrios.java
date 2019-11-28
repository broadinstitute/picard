/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.fingerprint;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.pedigree.PedFile;
import picard.pedigree.Sex;
import picard.util.MathUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static java.lang.Math.log10;

/**
 * Attempts to use a set of genome data samples to determine trio relationships within the cohort provided.
 * Roughly speaking it:
 * Determines genders
 * Performs pairwise tests between samples for 1st vs. 2nd degree relatedness
 * Tests all plausible trios from 1st degree related samples of appropriate genders
 *
 * @author Tim Fennell
 * @author Jonathan Barlev
 * @author Yossi Farjoun
 */
public abstract class ReconstructTrios extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output pedigree file name (PED format).")
    public File OUTPUT;

    @Argument(doc = "Output sex-only file name. If provided this file will be populated with the sex inference done internally.", optional = true)
    public File OUTPUT_SEX = null;

    @Argument(doc = "Input sex PED file. If provided this file will be used instead of trying to figure out the sex of the samples.", optional = true)
    public File INPUT_SEX = null;

    @Argument(doc = "Don't use sex to limit search.", optional = true)
    public Boolean NO_SEX = null;

    @Argument(shortName = "H", doc = "Haplotype map to be used for fingerprinting. FORMAT of Halotype map is, standard SAM header (with one @HD and multiple @SQs) " +
            "followed by a tab-separated header line and lines describing SNPs and Haplotype blocks. For example :" +
            "#CHROMOSOME    POSITION    NAME        MAJOR_ALLELE    MINOR_ALLELE    MAF         ANCHOR_SNP  PANELS" +
            "1              29478419    rs2503005   C               T               0.483791    rs2230678   panel1" +
            "1              29497820    rs1994859   C               T               0.483791    rs2230678   panel2" +
            "1              29500995    rs2486204   T               C               0.484167    rs2230678         " +
            "1              29509603    rs2428556   C               T               0.728976    rs2230679   panel1" +
            "" +
            "NOTE: SNPs listed with the same ANCHOR_SNP will be in the same haplotype. In case of discrepancy between the MAFs within a block, " +
            "the MAF of the first (smallest genomic position) SNP in the block is considered the MAF of the block." +
            "NOTE: the PANEL field is optional (as a value, not in the header)")
    public File HAPLOTYPE_MAP;

    @Argument(doc = "Minimum LOD score to achieve in the putative parent test to treat samples as possible parents for trio testing.")
    public double PARENT_ALONE_LOD = 2;

    @Argument(doc = "Minimum LOD difference between unrelated and the other models to stop the calculation short and call the relationship unrelated.")
    public double UNRELATED_SHORTCUT_LOD = 10;

    @Argument(doc = "Minimum LOD for the trio test in order to report the trio in the output ped file.")
    public double TRIO_LOD = 8;

    @Argument(shortName = "MC", doc = "List of possible names for male sex chromosome(s)")
    public Set<String> MALE_CHROMOSOMES = CollectionUtil.makeSet("Y", "chrY");

    @Argument(shortName = "FC", doc = "List of possible names for female sex chromosome(s)")
    public Set<String> FEMALE_CHROMOSOMES = CollectionUtil.makeSet("X", "chrX");

    @Argument(doc = "If given, will only try to reconstruct trios involving these samples", optional = true)
    public List<String> SAMPLE = null;

    @Argument(doc = "Same as SAMPLE but file containing samples to use", optional = true)
    public File SAMPLE_FILE = null;

    @Argument(doc = "Probability of denovo mutation. Includes the probability of an incorrect genotyping due to large structural variation and other incorrect modeling assumptions.")
    public double P_DENOVO = 0.000_01;

    @Argument(doc = "Emit all putative trios (that have putative mother and father) regardless of LOD and rank (for debugging purposes.)")
    public boolean EMIT_ALL_TRIOS = false;

    @Argument(doc = "Emit all samples into pedigree file (so that founding members get a line on their own).")
    public boolean EMIT_ALL_SAMPLES = false;


    protected Set<String> sexChromosomes;
    protected final Log log = Log.getInstance(ReconstructTrios.class);
    protected HaplotypeMap hapmap;

    protected abstract SexInferenceEngine getSexInferencer();

    private final int IDENTITY_CHECK_MIN_LOD = -3;

    @Override
    protected int doWork() {
        initializeMatrices();

        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (SAMPLE_FILE != null) IOUtil.assertFileIsReadable(SAMPLE_FILE);
        if (OUTPUT_SEX != null) IOUtil.assertFileIsWritable(OUTPUT_SEX);
        if (INPUT_SEX != null) IOUtil.assertFileIsReadable(INPUT_SEX);

        // Create the map of haploptypes without any on the sex chromosomes
        sexChromosomes = new HashSet<>();
        sexChromosomes.addAll(MALE_CHROMOSOMES);
        sexChromosomes.addAll(FEMALE_CHROMOSOMES);

        hapmap = new HaplotypeMap(HAPLOTYPE_MAP).withoutChromosomes(sexChromosomes);

        // Read all the fingerprints into per-sample FPs in memory
        log.info("Fingerprinting files");
        final Map<String, Fingerprint> fingerprints;

        if (SAMPLE_FILE != null) {
            try {
                SAMPLE = IOUtil.slurpLines(SAMPLE_FILE);
            } catch (FileNotFoundException e) {
                log.error(e);
                return -1;
            }
        }

        final Function<FingerprintIdDetails, String> idDetailsToSample = Fingerprint.getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE);

        final Map<String, Fingerprint> fingerprintsAll = Fingerprint.mergeFingerprintsBy(getFingerprints(), idDetailsToSample)
                .entrySet().stream()
                .collect(Collectors.toMap(e->e.getKey().getSample(), Map.Entry::getValue));

        if (SAMPLE != null && !SAMPLE.isEmpty()) {
            log.info("subsetting to " + SAMPLE.size() + " samples.");
            fingerprints = fingerprintsAll.entrySet()
                    .stream()
                    .filter(fpEntry -> SAMPLE.contains(fpEntry.getKey()))
                    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

            log.info("resulting map contains " + fingerprints.size() + " samples.");

        } else {
            fingerprints = fingerprintsAll;
        }

        if (fingerprints.size() < 3) {
            log.warn("Not enough samples to perform trio testing: " + fingerprints.size());
            new PedFile(true).write(OUTPUT); // Write out empty ped file.
            return 0;
        }

        // Determine gender for everyone that we can assign mother vs. father
        final Map<String, Sex> sampleSexes;
        if( NO_SEX ) {
            log.info("Not using sex. Will try both roles for every sample. Resulting genders in PED file are arbitrary!");
            sampleSexes = new CollectionUtil.DefaultingMap<>(Sex.Unknown);
        } else if (INPUT_SEX != null) {
            log.info("Reading sample sexes");
            final PedFile sampleSex = PedFile.fromFile(INPUT_SEX, false);
            sampleSexes = sampleSex.entrySet().stream()
                    .collect(Collectors.toMap(e -> e.getValue().getIndividualId(), e -> e.getValue().getSex()));

            SAMPLE.stream().filter(s -> !sampleSex.containsKey(s)).findAny().ifPresent(s -> {
                final String error = "Couldn't find sex for individual " + s + " in provided file " + INPUT_SEX;
                log.error(error);
                throw new RuntimeException(error);
            });

        } else {
            log.info("Determining sample sexes");
            sampleSexes = getSexInferencer().determineSexes();
        }

        if (OUTPUT_SEX != null) {
            log.info("writing sample sexes to output file.");
            PedFile.fromSexMap(sampleSexes).write(OUTPUT_SEX);
        }

        final ListMap<String, String> sampleToPotentialMothers = new ListMap<>();
        final ListMap<String, String> sampleToPotentialFathers = new ListMap<>();
        final ListMap<String, String> sampleToUngenderedParents = new ListMap<>();


        final NumberFormat fmt = new DecimalFormat("#.00");
        log.info("Looking for putative parents.");
        // Do a first pass and pick out samples that appear related to one another
        for (final Fingerprint offspring : fingerprints.values()) {
            for (final Fingerprint parent : fingerprints.values()) {
                if (offspring == parent) continue;

                final double lodParent = testParent(parent, offspring);
                final Sex parentSex = sampleSexes.get(parent.getSample());

                //otherwise debug messages get unruly
                if (lodParent > -30.0) {
                    log.debug(String.format("%s is %s of %s. LOD: %g.", parent.getSample(), (parentSex == Sex.Male ? "father" : "mother"), offspring.getSample(), lodParent));
                }
                //If putative parent seems more like a parent than a grandparent (and
                if (lodParent >= PARENT_ALONE_LOD) {
                    switch(parentSex ) {
                        case Female:
                            sampleToPotentialMothers.add(offspring.getSample(), parent.getSample());
                            break;
                        case Male:
                            sampleToPotentialFathers.add(offspring.getSample(), parent.getSample());
                            break;
                        case Unknown:
                            sampleToUngenderedParents.add(offspring.getSample(), parent.getSample());
                            break;
                        default:
                            throw new RuntimeException("Unpossible!");
                    }
                }
            }
            // Once all the putative parents have been collected, see if they all are identical (within gender sets).
            // if not, emit a warning since this seems a little odd. This entire loop only checks for warning, so will
            // be bypassed if warnings are not emitted.

            // TODO: fix this...
            // TODO: 1. not a problem since parent and child will show up here.
            // TODO: 2. need to account for ungendered parents.
            if (Log.isEnabled(Log.LogLevel.WARNING)) {
                for (final List<String> sampleToParent : CollectionUtil.makeList(
                        sampleToPotentialFathers.get(offspring.getSample()),
                        sampleToPotentialMothers.get(offspring.getSample()))) {

                    if (sampleToParent == null) {
                        continue;
                    } // no parents here, move along (makeList returns null in this case, not an empty list...)

                    for (int i = 0; i < sampleToParent.size(); i++) {
                        final String parent1 = sampleToParent.get(i);
                        final Fingerprint parent1FP = fingerprints.get(parent1);

                        for (int j = 0; j < i; j++) {
                            final String parent2 = sampleToParent.get(j);
                            final Fingerprint parent2FP = fingerprints.get(parent2);

                            final double lodIdentical = FingerprintChecker.calculateMatchResults(parent1FP, parent2FP).getLOD();
                            if (lodIdentical < IDENTITY_CHECK_MIN_LOD) { //  different individuals
                                log.warn(String.format("Found more than one distinct putative %s, for sample %s: %s, and, %s. LOD=%s", (sampleSexes.get(parent1) == Sex.Male ? "father" : "mother"), offspring.getSample(), parent1, parent2, fmt.format(lodIdentical)));
                            }
                        }
                    }
                }
            }
        }

        log.info("Looking for trios within possible families.");

        // Check for actual trios within families
        final PedFile pedfile = new PedFile(true);
        final Random rng = new Random();

        for (final Fingerprint offspring : fingerprints.values()) {
            final String offspringName = offspring.getSample();
            final Collection<String> potentialMothers = Optional.ofNullable(sampleToPotentialMothers.get(offspringName)).orElseGet(ArrayList::new);
            potentialMothers.addAll(Optional.ofNullable(sampleToUngenderedParents.get(offspringName)).orElseGet(ArrayList::new));

            final Collection<String> potentialFathers = Optional.ofNullable(sampleToPotentialFathers.get(offspringName)).orElseGet(ArrayList::new);
            potentialFathers.addAll(Optional.ofNullable(sampleToUngenderedParents.get(offspringName)).orElseGet(ArrayList::new));

            final SortedSet<TrioTestResults> results = new TreeSet<>();

            if (potentialMothers.isEmpty() || potentialFathers.isEmpty()) {
                log.debug(String.format("sample %s doesn't have a putative %s", offspringName, potentialMothers.isEmpty() ? "mother" : "father"));
                continue;
            }

            for (final String motherName : potentialMothers) {
                // once a sample has been considered as a mother, do not consider it as a father
                potentialFathers.remove(motherName);
                for (final String fatherName : potentialFathers) {

                    final Fingerprint mother = fingerprints.get(motherName);
                    final Fingerprint father = fingerprints.get(fatherName);

                    final TrioTestResults result = testTrios(mother, father, offspring);
                    log.debug(String.format("Trio test results for  %s + %s -> %s are llTrio=%g, llFatherPlusRandom=%g, llMotherPlusRandom=%g, llRandomParents=%g, llSibsAsParents=%g, llThreeGenChain=%g, lodTrio=%g",
                            mother.getSample(), father.getSample(), offspring.getSample(), result.llTrio, result.llFatherPlusRandom, result.llMotherPlusRandom, result.llRandomParents, result.llSibsAsParents, result.llThreeGenChain, result.lodTrio));

                    if (result.lodTrio > Math.min(TRIO_LOD, 0)) results.add(result);
                }
            }

            log.debug(results.size() + " trio matches for offspring" + offspringName);
            if (!results.isEmpty()) {

                final TrioTestResults best = results.first();
                if (best.lodTrio >= TRIO_LOD) {
                    log.debug(String.format("Best match for %s has LOD %s for Father=%s, and Mother=%s.", best.offspring, fmt.format(best.lodTrio), best.father, best.mother));

                    results.stream()
                            .limit(EMIT_ALL_TRIOS ? Integer.MAX_VALUE : 1)
                            .filter(t -> EMIT_ALL_TRIOS || t.lodTrio > TRIO_LOD)
                            .forEach(trio -> pedfile.put(trio.offspring + "_" + rng.nextInt(),pedfile.new PedTrio(trio.offspring, trio.offspring, trio.father, trio.mother,
                                    sampleSexes.get(trio.offspring), trio.lodTrio)));
                }
            }
        }

        if (EMIT_ALL_SAMPLES) {
            log.info("adding founding members");

            final Set<String> offsprings = pedfile.values().stream()
                    .map(PedFile.PedTrio::getIndividualId)
                    .collect(Collectors.toSet());

            fingerprints.keySet().stream()
                    .filter(s -> !offsprings.contains(s))
                    .map(sample -> pedfile.new PedTrio(sample, sample, null, null, sampleSexes.get(sample), -9))
                    .forEach(pedfile::add);
        }

        log.info("writing pedigree to file");
        pedfile.write(OUTPUT);

        return 0;
    }

    /**
     * Tests whether the three samples appear to constitute a mother/father/offspring trio. Calculates the
     * likelihood of this arrangement vs. the likelihoods of the three samples being unrelated random samples
     * and also of one of the parents being a true parent and one of the samples being a random unrelated sample.
     *
     * @param mother    the fingerprint of the putative mother
     * @param father    the fingerprint of the putative father
     * @param offspring the fingerprint of the offspring
     * @return a TrioTestResults object that has the likelihoods of each scenarios along with the names of the samples
     */
    private TrioTestResults testTrios(final Fingerprint mother, final Fingerprint father, final Fingerprint offspring) {
        double llTrio = 0, llMotherPlusRandom = 0, llFatherPlusRandom = 0, llRandomParents = 0, llSibsAsParents = 0;
        double llThreeGenerationChain = 0;


        for (final HaplotypeBlock hap : mother.keySet()) {
            final HaplotypeProbabilities mom = mother.get(hap);
            final HaplotypeProbabilities dad = father.get(hap);
            final HaplotypeProbabilities os = offspring.get(hap);

            // Skip any haplotypes where we don't have information for all three people
            if (mom == null || dad == null || os == null) continue;

            final double[] probsMom = mom.getPosteriorProbabilities();
            final double[] probsDad = dad.getPosteriorProbabilities();
            final double[] random = hap.getHaplotypeFrequencies();

            llTrio += llChildofKnownParents(probsMom, probsDad, os);
            llMotherPlusRandom += llChildofKnownParents(probsMom, random, os);
            llFatherPlusRandom += llChildofKnownParents(random, probsDad, os);
            llRandomParents += llChildofKnownParents(random, random, os);

            // here we are giving the putative parents as if they are the children, since that's what we
            // want to protect against.
            llSibsAsParents += llParentOfKnownChildren(probsDad, probsMom, os);

            // protect against a scenario where one of the parents is actually the offspring of the child
            llThreeGenerationChain += Math.max( llMiddleOfChain(probsDad, os, probsMom), llMiddleOfChain(probsMom, os, probsDad));

            final double lodTrio=llChildofKnownParents(probsMom, probsDad, os) - MathUtil.max(new double[] {
                    llChildofKnownParents(random, probsDad, os), llChildofKnownParents(probsMom, random, os),
                    llChildofKnownParents(random, random, os), llParentOfKnownChildren(probsDad, probsMom, os),
                    llMiddleOfChain(probsDad, os, probsMom), llMiddleOfChain(probsMom, os, probsDad)});

            log.debug(String.format("Trio test for  %s + %s -> %s  at Haplotype %s are llTrio=%g, llFatherPlusRandom=%g, llMotherPlusRandom=%g, llRandomParents=%g, llSibsAsParents=%g, llThreeGenChainDadisGramps=%g, llThreeGenChainDadisSon=%g lodTrio=%g",
                    mother.getSample(), father.getSample(), offspring.getSample(), hap.toString(), llChildofKnownParents(probsMom, probsDad, os),
                    llChildofKnownParents(random, probsDad, os), llChildofKnownParents(probsMom, random, os),
                    llChildofKnownParents(random, random, os), llParentOfKnownChildren(probsDad, probsMom, os),
                    llMiddleOfChain(probsDad, os, probsMom), llMiddleOfChain(probsMom, os, probsDad), lodTrio));
        }

        return new TrioTestResults(offspring.getSample(), father.getSample(), mother.getSample(),
                llTrio, llMotherPlusRandom, llFatherPlusRandom, llRandomParents, llSibsAsParents, llThreeGenerationChain);
    }

    /**
     * Little class to encapsulate the results of testing three samples that may or may not form a Trio.
     */
    private static class TrioTestResults implements Comparable<TrioTestResults> {
        public TrioTestResults(final String offspring, final String father, final String mother,
                               final double llTrio, final double llMotherPlusRandom,
                               final double llFatherPlusRandom, final double llRandomParents, double llSibsAsParents, double llThreeGenChain) {
            this.offspring = offspring;
            this.father = father;
            this.mother = mother;
            this.llTrio = llTrio;
            this.llMotherPlusRandom = llMotherPlusRandom;
            this.llFatherPlusRandom = llFatherPlusRandom;
            this.llRandomParents = llRandomParents;
            this.llSibsAsParents = llSibsAsParents;
            this.llThreeGenChain = llThreeGenChain;

            this.lodTrio = llTrio - MathUtil.max(new double[]{llMotherPlusRandom, llFatherPlusRandom, llRandomParents, llSibsAsParents, llThreeGenChain});
        }

        public final String offspring;

        public final String father;
        public final String mother;
        public final double llTrio;
        public final double llMotherPlusRandom;
        public final double llFatherPlusRandom;
        public final double llRandomParents;
        public final double llSibsAsParents;
        public final double llThreeGenChain;

        public final double lodTrio;

        /**
         * Comparison function that orders things with higher LODs first.
         */
        @Override
        public int compareTo(final TrioTestResults that) {
            // Multiply by 1000 so that the conversion to int only truncates after 3 significant digits
            return (int) ((that.lodTrio * 1000) - (this.lodTrio * 1000));
        }
    }

    /**
     * Attempts to determine if samples have 50% sharing of genetic material (and therefore could be parent/offspring)
     * or are more likely to be 25% related or unrelated.
     * <p/>
     * Works as follows:
     * - Assume one FP comes from the parent and ask what the likelihood of the two fingerprints is
     * given the model that one sample is a parent of the other sample.
     * - Assume that the putative parent is actually a grandparent and synthesize the likelihoods of
     * postulative - parent and perform the same test as above
     * - Lastly do the same but using the population frequencies as the likelihood of a "random" parent
     */
    private double testParent(final Fingerprint parentFp, final Fingerprint offspringFp) {
        double llUnrelated = 0, llIdentical = 0, ll50Related = 0, ll25Related = 0;

        for (final HaplotypeBlock hap : parentFp.keySet()) {
            final HaplotypeProbabilities parent = parentFp.get(hap);
            final HaplotypeProbabilities offspring = offspringFp.get(hap);

            if (parent == null || offspring == null) continue;
            final double[] pParentModel = parent.getPosteriorProbabilities();
            final double[] pPopulationModel = hap.getHaplotypeFrequencies();

            // Calculate a set of probabilities for a hypothetical parent who is themselves the offspring of the putative "parent"
            // passed in and a random other sample from the population
            final double[] pPseudoParentModel = getModelFromParents(pParentModel, pPopulationModel);

            // Now add to our likelihoods
            llUnrelated += offspring.shiftedLogEvidenceProbability();
            llIdentical += offspring.shiftedLogEvidenceProbabilityUsingGenotypeFrequencies(pParentModel);

            ll50Related += llChildofKnownParents(pParentModel, pPopulationModel, offspring);
            ll25Related += llChildofKnownParents(pPseudoParentModel, pPopulationModel, offspring);

            // if the unrelated model has grown so much as to overshadow the other models, call the situation unrelated and stop the calculation
            // This is the slowest part of the program, so looking for ways to speed it up.

            if (llUnrelated > MathUtil.max(new double[]{llIdentical, ll25Related, ll50Related}) + UNRELATED_SHORTCUT_LOD) {
                break;
            }
        }
        //The logic here is that the likelihood of the trio is compared to the most likely other option for the offspring and the ratio
        //of those two likelihoods is used for calculating the LOD
        final double LodParent = ll50Related - MathUtil.max(new double[]{llIdentical, ll25Related, llUnrelated});
        log.debug(String.format("Parent test results for  %s -> %s are ll50Related=%g, llIdentical=%g, ll25Related=%g, llUnrelated=%g. LOD=%g",
                parentFp.getSample(),offspringFp.getSample(),ll50Related, llIdentical, ll25Related, llUnrelated, LodParent));

        return LodParent;
    }

    private static double[][] getDenovoMatrix(final double delta){
        final double tau = 1 - delta;
        final double tau_delta = tau * delta;
        final double tau_2 = tau * tau;
        final double delta_2 = delta * delta;
        return new double[][] {
                {tau_2, tau_delta, delta_2},
                {2 * tau_delta, delta_2 + tau_2, 2 * tau_delta},
                {delta_2, tau_delta, tau_2}
        };
    }

    // Mendelian matrix. M[a][b][c] is the probability that
    // genotype a come from parents with genotype b and c
    private double [][][] MWithDenovo = null;
    private double [][][][] M2WithDenovo = null;


    private static double[][][] getMWithDenovoMutations(final double delta) {
        final double[][][] M = {
            {
                    {1, .5, 0}, //hom_ref child
                    {.5, .25, 0},
                    {0, 0, 0}
            },
            {
                    {0, .5, 1},//het child
                    {.5, .5, .5},
                    {1, .5, 0}
            },
            {
                    {0, 0, 0},//hom_var child
                    {0, .25, .5},
                    {0, .5, 1}
            }
        };

        final double [][][] retval = new double[3][][];
        final double [][] denovoMatrix = getDenovoMatrix(delta);

        for (final int child : genotypes) {
            retval[child] = new double[genotypes.length][];
            for (final int parent1 : genotypes) {
                retval[child][parent1] = new double[genotypes.length];
                for (final int parent2 : genotypes) {
                    for (final int altchild : genotypes) {
                        retval[child][parent1][parent2] += denovoMatrix[child][altchild] * M[altchild][parent1][parent2];
                    }
                }
            }
        }
        return retval;
    }

    private void initializeMatrices(){
        MWithDenovo = getMWithDenovoMutations(P_DENOVO);
        M2WithDenovo = getParentsAndTwoChildrenTensor(MWithDenovo);
    }

    static public int[] genotypes = {0, 1, 2};

    /**
     * Returns the probabilities of each of the three genotypes given a random mating between two parents with specified probabilities.
     */
    private double[] getModelFromParents(final double[] pModel1, final double[] pModel2) {
        final double[] pPriorOffspring = new double[genotypes.length];

        for (final int child : genotypes) {
            for (final int parent1 : genotypes) {
                for (final int parent2 : genotypes) {
                    pPriorOffspring[child] += pModel1[parent1] * pModel2[parent2] * MWithDenovo[child][parent1][parent2];
                }
            }
        }
        return pPriorOffspring;
    }

    /** M^2 the tensor for two parents and two children  */
    private static double[][][][] getParentsAndTwoChildrenTensor(double[][][] M) {
        final int n = genotypes.length;
        final double[][][][] parentsAndTwoChildrenTensor = new double[n][n][n][n];

        for (final int child1 : genotypes) {
            for (final int child2 : genotypes) {
                for (final int parent1 : genotypes) {
                    for (final int parent2 : genotypes) {
                        parentsAndTwoChildrenTensor[parent1][parent2][child1][child2] = M[child1][parent1][parent2] * M[child2][parent1][parent2];
                    }
                }
            }
        }
        return parentsAndTwoChildrenTensor;
    }


    /**
     * Calculates the likelihood of the offspring being the product of two parents by synthesizing the genotype probabilities
     * for an offspring given the parent's genotypes and then comparing that to the actual offspring's genotypes likelihoods.
     */
    private double llChildofKnownParents(final double[] parent1, final double[] parent2, final HaplotypeProbabilities offspring) {
        final double[] predictedOffspring = getModelFromParents(parent1, parent2);
        return offspring.shiftedLogEvidenceProbabilityUsingGenotypeFrequencies(predictedOffspring);
    }

    /**
     * Calculates the likelihood that the data is from the parent of both given genotypes
     */
    private double pParentOfKnownChildren(final double[] child1, final double[] child2, final HaplotypeProbabilities parentHP) {

        double probability = 0;
        final double [] parent = parentHP.getLikelihoods();
        final double [] prior = parentHP.getPriorProbablities();
        for (int child1G : genotypes) {
            for (int child2G : genotypes) {
                for (int parent1G : genotypes) {
                    for (int parent2G : genotypes) {
                        probability += M2WithDenovo[parent1G][parent2G][child1G][child2G] * child1[child1G] * child2[child2G] * parent[parent1G] * prior[parent2G];
                    }
                }
            }
        }
        return probability;
    }

    /**
     * Calculates the likelihood that the data is from the parent of both given genotypes
     */
    private double llParentOfKnownChildren(final double[] child1, final double[] child2, final HaplotypeProbabilities parentHP) {
        return log10(pParentOfKnownChildren(child1, child2, parentHP));
    }

    /**
     * calculates the likelihood that the data of child agrees with a model where the parent and offspring are given
     * @param parent probabilities of parent of subject
     * @param childHP probabilities of subject
     * @param grandchild probabilities of putative child of subject
     * @return likelihood that data would arise from such an arrangement
     */
    private double pMiddleOfChain(final double[] parent, final HaplotypeProbabilities childHP, final double[] grandchild) {
        double probability = 0;
        final double [] child = childHP.getLikelihoods();
        final double [] prior = childHP.getPriorProbablities();
        for (int parent1G : genotypes) {
            for (int unknownSpouce1G : genotypes) {
                for (int childG : genotypes) {
                    for (int unknownSpouce2G : genotypes) {
                        for (int grandChildG : genotypes) {
                            probability +=
                                    MWithDenovo[childG][parent1G][unknownSpouce1G] * prior[unknownSpouce1G] *
                                    MWithDenovo[grandChildG][childG][unknownSpouce2G] * prior[unknownSpouce2G] * parent[parent1G] * child[childG] * grandchild[grandChildG];
                        }
                    }
                }
            }
        }
        return probability;
     }

    /**
     * calculates the log likelihood that the data of child agrees with a model where the parent and offspring are given
     * @param parent probabilities of parent of subject
     * @param childHP probabilities of subject
     * @param grandchild probabilities of putative child of subject
     * @return probability that data would arise from such an arrangement
     */
    private double llMiddleOfChain(final double[] parent, final HaplotypeProbabilities childHP, final double[] grandchild) {
        return log10(pMiddleOfChain(parent, childHP, grandchild));
    }

    protected abstract Map<FingerprintIdDetails, Fingerprint> getFingerprints();
}

