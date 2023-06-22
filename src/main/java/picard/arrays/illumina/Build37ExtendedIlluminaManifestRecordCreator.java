package picard.arrays.illumina;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import picard.PicardException;

import java.io.File;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Build37ExtendedIlluminaManifestRecordCreator {

    private final String targetBuild;
    private final Map<String, ReferenceSequenceFile> referenceFilesMap;
    private final Map<String, File> chainFilesMap;
    private final Set<String> supportedBuilds;

    // Builds found in records that we don't support (and so are failed)
    private final Set<String> unsupportedBuilds;

    private boolean refStrandDefinedInManifest = true;

    private final Log log = Log.getInstance(Build37ExtendedIlluminaManifestRecordCreator.class);

    public static final String BUILD_37 = "37";

    // These are the IUPAC nucleotide codes as described here: https://www.bioinformatics.org/sms/iupac.html
    public static final String IUPAC_NUCLEOTIDE_CODES = "ACGTRYSWKMBDHVNacgtryswkmbdhvn";
    public static final String ACGT_CODES = "ACGTacgt";

    private static final String SRC_SEQ_REGEX = "([" + IUPAC_NUCLEOTIDE_CODES + "]*)\\[([" + ACGT_CODES + "-])\\/([" + ACGT_CODES + "]*)\\]([" + IUPAC_NUCLEOTIDE_CODES + "]*)";
    private static final Pattern pattern = Pattern.compile(SRC_SEQ_REGEX);

    private static final String ACGT_REGEX = "^[" + ACGT_CODES + "]+$";
    private static final Pattern ACGT_PATTERN = Pattern.compile(ACGT_REGEX);

    // Symbolics for the regex groups...
    public static final int FIVE_PRIME_SEQUENCE = 1;
    public static final int PRE_INDEL_SEQUENCE = 2;
    public static final int INDEL_SEQUENCE = 3;
    public static final int THREE_PRIME_SEQUENCE = 4;

    Build37ExtendedIlluminaManifestRecordCreator(final String targetBuild,
                                                 final Map<String, ReferenceSequenceFile> referenceFilesMap,
                                                 final Map<String, File> chainFilesMap) {
        this.targetBuild = targetBuild;
        this.referenceFilesMap = referenceFilesMap;
        this.chainFilesMap = chainFilesMap;
        this.supportedBuilds = new HashSet<>();
        this.unsupportedBuilds = new TreeSet<>();
        supportedBuilds.add(targetBuild);
        supportedBuilds.addAll(chainFilesMap.keySet());
    }

    public Build37ExtendedIlluminaManifestRecord createRecord(
              final IlluminaManifestRecord illuminaManifestRecord) {

        Build37ExtendedIlluminaManifestRecord newRecord = new Build37ExtendedIlluminaManifestRecord(illuminaManifestRecord,
                Build37ExtendedIlluminaManifestRecord.Flag.PASS,
                "",
                null,
                "",
                "",
                "",
                "");

        // Look for entries which Illumina has marked as invalid
        if ((illuminaManifestRecord.getChr().equals(IlluminaManifestRecord.ILLUMINA_FLAGGED_BAD_CHR)) ||
            (illuminaManifestRecord.getPosition() == 0)) {
            newRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.ILLUMINA_FLAGGED;
            return newRecord;
        }

        if (!this.supportedBuilds.contains(illuminaManifestRecord.getMajorGenomeBuild())) {
            this.unsupportedBuilds.add(illuminaManifestRecord.getMajorGenomeBuild());
            newRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.UNSUPPORTED_GENOME_BUILD;
            return newRecord;
        }

        if (illuminaManifestRecord.getMajorGenomeBuild().equals(targetBuild)) {
            // no liftover needed
            newRecord.b37Chr = illuminaManifestRecord.getChr();
            newRecord.b37Pos = illuminaManifestRecord.getPosition();
        } else {
            liftOverToTargetBuild(newRecord, illuminaManifestRecord);
            if (newRecord.isFail()) {
                return newRecord;
            }
        }

        final ReferenceSequenceFile refFile = referenceFilesMap.get(targetBuild);

        setReferenceStrand(newRecord, refFile);

        if (!newRecord.isFail()) {
            if (newRecord.isSnp()) {
                processSnp(newRecord, refFile);
            } else {
                processIndel(newRecord, refFile);
            }
        }

        return newRecord;
    }

    /**
     * Uses source sequence to determine snpAlleleA, snpAlleleB, and Allele.
     * <p>
     */
    private void processSnp(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                            final ReferenceSequenceFile refFile) {
        build37ExtendedIlluminaManifestRecord.snpAlleleA = build37ExtendedIlluminaManifestRecord.getSnp().substring(1, 2);
        build37ExtendedIlluminaManifestRecord.snpAlleleB = build37ExtendedIlluminaManifestRecord.getSnp().substring(3, 4);

        if (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.NEGATIVE) {
            build37ExtendedIlluminaManifestRecord.snpAlleleA = SequenceUtil.reverseComplement(build37ExtendedIlluminaManifestRecord.snpAlleleA);
            build37ExtendedIlluminaManifestRecord.snpAlleleB = SequenceUtil.reverseComplement(build37ExtendedIlluminaManifestRecord.snpAlleleB);
        }

        if (build37ExtendedIlluminaManifestRecord.isAmbiguous()) {
            if (build37ExtendedIlluminaManifestRecord.getAlleleBProbeSeq() != null) {
                final String probeAAllele = build37ExtendedIlluminaManifestRecord.getAlleleAProbeSeq().substring(build37ExtendedIlluminaManifestRecord.getAlleleAProbeSeq().length() - 1);
                validateThatSequenceOnlyContainsACGTCharacters("AlleleAProbeSeq for record: " + build37ExtendedIlluminaManifestRecord, probeAAllele);
                final String probeBAllele = build37ExtendedIlluminaManifestRecord.getAlleleBProbeSeq().substring(build37ExtendedIlluminaManifestRecord.getAlleleBProbeSeq().length() - 1);
                validateThatSequenceOnlyContainsACGTCharacters("AlleleBProbeSeq for record: " + build37ExtendedIlluminaManifestRecord, probeBAllele);
                if (!probeAAllele.equals(build37ExtendedIlluminaManifestRecord.snpAlleleA) &&
                    !probeBAllele.equals(build37ExtendedIlluminaManifestRecord.snpAlleleB) &&
                        (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.POSITIVE)) {
                    build37ExtendedIlluminaManifestRecord.snpAlleleA = probeAAllele;
                    build37ExtendedIlluminaManifestRecord.snpAlleleB = probeBAllele;
                }
            } else {
                // This manifest contains no Allele B Probe Sequence.
                build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.MISSING_ALLELE_B_PROBESEQ;
                log.warn("Error in processSnp.  Record:" + build37ExtendedIlluminaManifestRecord);
                log.warn("  Ambiguous probe without alleleBProbeSeq");
            }
        }

        build37ExtendedIlluminaManifestRecord.snpRefAllele = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos, build37ExtendedIlluminaManifestRecord.b37Pos);
    }

    /**
     * This method sets the reference strand
     *
     * If the refStrand is provided in the Illumina manifest(s) we will use that.
     * Otherwise, we will attempt to calculate the reference strand.
     *
     * @param refFile reference to use for finding the probe sequence
     */
    private void setReferenceStrand(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                                    final ReferenceSequenceFile refFile) {
        if (build37ExtendedIlluminaManifestRecord.referenceStrand != null) {
            return;
        }

        build37ExtendedIlluminaManifestRecord.referenceStrand = build37ExtendedIlluminaManifestRecord.getRefStrand();
        if (build37ExtendedIlluminaManifestRecord.referenceStrand == Strand.NONE) {
            // If here, the referenceStrand is not defined in the manifest.
            refStrandDefinedInManifest = false;

            if (build37ExtendedIlluminaManifestRecord.isSnp()) {

                String probeSeq = build37ExtendedIlluminaManifestRecord.getAlleleAProbeSeq();
                if (build37ExtendedIlluminaManifestRecord.isAmbiguous()) {
                    probeSeq = probeSeq.substring(0, probeSeq.length() - 1);        // Ambiguous SNPs contain the probed base, so truncate
                }
                validateThatSequenceOnlyContainsACGTCharacters("AlleleAProbeSeq for record: " + build37ExtendedIlluminaManifestRecord, probeSeq);

                final String reference = getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos - probeSeq.length(), build37ExtendedIlluminaManifestRecord.b37Pos - 1);
                final String reverseReference = SequenceUtil.reverseComplement(getSequenceAt(refFile, build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos + 1, build37ExtendedIlluminaManifestRecord.b37Pos + probeSeq.length()));

                if (reference.equals(probeSeq)) {
                    build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.POSITIVE;
                } else if (reverseReference.equals(probeSeq)) {
                    build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.NEGATIVE;
                } else {
                    build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.PROBE_SEQUENCE_MISMATCH;
                    log.warn("Error in getStrand.  Record:" + build37ExtendedIlluminaManifestRecord);
                    log.warn("  Couldn't find alleleAProbeSeq in reference");
                    log.debug("  AlleleAProbeSeq: " + build37ExtendedIlluminaManifestRecord.getAlleleAProbeSeq());
                    log.debug("  Reference:       " + reference);
                    log.debug("  Reverse Ref:     " + reverseReference);
                }
            } else {
                // Some Illumina manifests do not have refStrand defined.  We will use the illumina strand instead.
                if (build37ExtendedIlluminaManifestRecord.getIlmnStrand() == IlluminaManifestRecord.IlluminaStrand.PLUS)  {
                    build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.POSITIVE;
                }
                else if (build37ExtendedIlluminaManifestRecord.getIlmnStrand() == IlluminaManifestRecord.IlluminaStrand.MINUS) {
                    build37ExtendedIlluminaManifestRecord.referenceStrand = Strand.NEGATIVE;
                }
                else {
                    throw new PicardException("Unexpected value for Illumina Strand: " + build37ExtendedIlluminaManifestRecord.getIlmnStrand());
                }
            }
        }
    }


    private void processIndel(final Build37ExtendedIlluminaManifestRecord extendedIlluminaManifestRecord,
                              final ReferenceSequenceFile refFile) {
        if (extendedIlluminaManifestRecord.isSnp()) {
            throw new PicardException("This shouldn't happen");
        }

        // Validate the source sequence
        final Matcher matcher = parseSourceSeq(extendedIlluminaManifestRecord.getSourceSeq());
        if (!matcher.group(PRE_INDEL_SEQUENCE).equals("-")) {       // In indels it's always of the form: [-/GCA]
            throw new PicardException("Unexpected allele '-' Record: " + extendedIlluminaManifestRecord);
        }

        String fivePrimeSeq = matcher.group(FIVE_PRIME_SEQUENCE).toUpperCase();
        String indelSeq = matcher.group(INDEL_SEQUENCE).toUpperCase();
        String threePrimeSeq = matcher.group(THREE_PRIME_SEQUENCE).toUpperCase();

        validateThatSequenceOnlyContainsACGTCharacters("Indel sequence for record: " + extendedIlluminaManifestRecord, indelSeq);

        boolean isSourceOnDesignStrand = extendedIlluminaManifestRecord.getSourceStrand() == extendedIlluminaManifestRecord.getIlmnStrand();
        if (isSourceOnDesignStrand != (extendedIlluminaManifestRecord.referenceStrand == Strand.POSITIVE)) {
            final String temp = threePrimeSeq;
            threePrimeSeq = SequenceUtil.reverseComplement(fivePrimeSeq);
            indelSeq = SequenceUtil.reverseComplement(indelSeq);
            fivePrimeSeq = SequenceUtil.reverseComplement(temp);
        }

        ImmutablePair<Boolean, Boolean> illuminaIsIndel = calculateIsInsertionOrDeletion(extendedIlluminaManifestRecord, refFile, fivePrimeSeq, indelSeq, threePrimeSeq);
        boolean isInsertion = illuminaIsIndel.left;
        boolean isDeletion = illuminaIsIndel.right;

        if (!isInsertion && !isDeletion) {
            extendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.INDEL_NOT_FOUND;
            log.warn("Error in processIndel.  Record: " + extendedIlluminaManifestRecord);
            log.warn("  Couldn't find source sequence with or without variant in reference");
            return;
        }

        if (isInsertion && isDeletion) {
            extendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.INDEL_CONFLICT;
            log.warn("Error in processIndel.  Record: " + extendedIlluminaManifestRecord);
            log.warn("  Conflict.  Both source sequence with and without variation found in reference");
            return;
        }

        if (isDeletion) {
            // A deletion.  Position in VCF is before the deletion
            extendedIlluminaManifestRecord.b37Pos--;
        }
        final String refAllele = getSequenceAt(refFile, extendedIlluminaManifestRecord.getB37Chr(), extendedIlluminaManifestRecord.getB37Pos(), extendedIlluminaManifestRecord.getB37Pos());
        if (isDeletion) {
            extendedIlluminaManifestRecord.snpRefAllele = refAllele + indelSeq;
        } else {
            extendedIlluminaManifestRecord.snpRefAllele = refAllele;
        }
        if (extendedIlluminaManifestRecord.getSnp().equals("[I/D]")) {
            extendedIlluminaManifestRecord.snpAlleleA = refAllele + indelSeq;
            extendedIlluminaManifestRecord.snpAlleleB = refAllele;
        } else {
            extendedIlluminaManifestRecord.snpAlleleA = refAllele;
            extendedIlluminaManifestRecord.snpAlleleB = refAllele + indelSeq;
        }
    }

    private ImmutablePair<Boolean, Boolean> calculateIsInsertionOrDeletion(final Build37ExtendedIlluminaManifestRecord extendedIlluminaManifestRecord,
                                                                          final ReferenceSequenceFile refFile,
                                                                          final String fivePrimeSeq, final String indelSeq, final String threePrimeSeq) {

        // Do Illumina's strange left shift http://github.com/Illumina/GTCtoVCF/blob/develop/BPMRecord.py
        ImmutablePair<String, String> leftShiftedSourceSeq = illuminaLeftShift(fivePrimeSeq, indelSeq, threePrimeSeq);
        final String leftShiftedFivePrimeSeq = leftShiftedSourceSeq.left;
        final String leftShiftedThreePrimeSeq = leftShiftedSourceSeq.right;

        int startIndex = extendedIlluminaManifestRecord.b37Pos;
        String indelSequenceFromReference = getSequenceAt(refFile, extendedIlluminaManifestRecord.b37Chr, startIndex, startIndex + indelSeq.length() - 1);

        boolean indelSequenceMatch = indelSeq.equals(indelSequenceFromReference);

        // CONTAINS the indel
        String genomicDeletionFivePrime = getSequenceAt(refFile, extendedIlluminaManifestRecord.b37Chr, startIndex - leftShiftedFivePrimeSeq.length(), startIndex - 1);
        String genomicDeletionThreePrime = getSequenceAt(refFile, extendedIlluminaManifestRecord.b37Chr, startIndex + indelSeq.length(), startIndex + indelSeq.length() + leftShiftedThreePrimeSeq.length() - 1);

        // Calculate the deletionContextScore:
        ImmutablePair<Double, Integer> deletionContextInfo = calculateIsDeletion(leftShiftedSourceSeq.left, leftShiftedSourceSeq.right,
                indelSequenceMatch, genomicDeletionFivePrime, indelSeq, genomicDeletionThreePrime);

        String genomicInsertionFivePrime = getSequenceAt(refFile, extendedIlluminaManifestRecord.b37Chr, startIndex - leftShiftedFivePrimeSeq.length() + 1, startIndex);
        String genomicInsertionThreePrime = getSequenceAt(refFile, extendedIlluminaManifestRecord.b37Chr, startIndex + 1, startIndex + leftShiftedThreePrimeSeq.length());

        // Calculate the insertionContextScore:
        ImmutablePair<Double, Integer> insertionContextInfo = calculateIsInsertion(leftShiftedSourceSeq.left, leftShiftedSourceSeq.right,
                genomicInsertionFivePrime, indelSeq, genomicInsertionThreePrime);

        boolean illuminaIsDeletion = indelSequenceMatch & (deletionContextInfo.left > insertionContextInfo.left) & (deletionContextInfo.right >= 1);

        boolean illuminaIsInsertion = insertionContextInfo.left > deletionContextInfo.left & (insertionContextInfo.right >= 1);

        return new ImmutablePair<>(illuminaIsInsertion, illuminaIsDeletion);
    }

    /**
     * Calculates the deletion context score and minimum deletion context, used in determining whether an indel, as encoded
     * in the Illumina manifest represents a genomic deletion
     * (that is, the reference contains the indel, and this assay is testing for a deletion)
     *
     * @param fivePrimeSeq The Five Prime component of the Source Sequence, as pulled from the Illumina Manifest (and left aligned)
     * @param threePrimeSeq The Three Prime component of the Source Sequence, as pulled from the Illumina Manifest (and left aligned)
     * @param indelSequenceMatch A boolean indicating whether the indel sequence (as pulled from the Illumina Manifest) is found in the reference sequence
     * @param genomicDeletionFivePrime The Five Prime component of the source sequence as pulled from reference WITH the indel
     * @param genomicDeletionThreePrime The Three Prime component of the source sequence as pulled from reference WITH the indel
     */
    static ImmutablePair<Double, Integer> calculateIsDeletion(final String fivePrimeSeq, final String threePrimeSeq, final boolean indelSequenceMatch,
                                        final String genomicDeletionFivePrime, final String indelSeq, final String genomicDeletionThreePrime) {

        ImmutablePair<String, String> leftShiftedGenomicDeletionSeqs = illuminaLeftShift(genomicDeletionFivePrime, indelSeq, genomicDeletionThreePrime);

        // Calculate the deletionContextScore:
        int deletionContextSuffixMatchLength = commonSuffixLength(leftShiftedGenomicDeletionSeqs.left, fivePrimeSeq);
        int deletionContextPrefixMatchLength = commonPrefixLength(leftShiftedGenomicDeletionSeqs.right, threePrimeSeq);
        final int minDeletionContextMaxLength = Math.min(deletionContextPrefixMatchLength, deletionContextSuffixMatchLength);

        int maxDeletionContext = Math.min(genomicDeletionFivePrime.length(), fivePrimeSeq.length()) +
                                 Math.min(genomicDeletionThreePrime.length(), threePrimeSeq.length()) +
                                 indelSeq.length();

        double deletionContextScore = 0.0;
        if (indelSequenceMatch) {
            deletionContextScore = ((double) (deletionContextPrefixMatchLength + deletionContextSuffixMatchLength + indelSeq.length())) / maxDeletionContext;
        }
        return new ImmutablePair<>(deletionContextScore, minDeletionContextMaxLength);
    }

    /**
     *
     * Calculates the insertion context score and minimum insertion context, used in determining whether an indel, as encoded
     * in the Illumina manifest represents a genomic insertion
     * (that is, the reference DOES NOT contain the indel, and this assay is testing for an insertion)
     *
     * @param fivePrimeSeq The Five Prime component of the Source Sequence, as pulled from the Illumina Manifest (and left aligned)
     * @param threePrimeSeq The Three Prime component of the Source Sequence, as pulled from the Illumina Manifest (and left aligned)
     * @param genomicInsertionFivePrime The Five Prime component of the source sequence as pulled from reference WITH the indel
     * @param indelSeq The sequence of the indel
     * @param genomicInsertionThreePrime The Three Prime component of the source sequence as pulled from reference WITH the indel
     */
    static ImmutablePair<Double, Integer> calculateIsInsertion(final String fivePrimeSeq, final String threePrimeSeq,
                                                               final String genomicInsertionFivePrime, final String indelSeq, final String genomicInsertionThreePrime) {
        ImmutablePair<String, String> leftShiftedGenomicInsertionSeqs = illuminaLeftShift(genomicInsertionFivePrime, indelSeq, genomicInsertionThreePrime);

        // Calculate the insertionContextScore:
        int insertionContextSuffixMatchLength = commonSuffixLength(leftShiftedGenomicInsertionSeqs.left, fivePrimeSeq);
        int insertionContextPrefixMatchLength = commonPrefixLength(leftShiftedGenomicInsertionSeqs.right, threePrimeSeq);
        final int minInsertionContextMaxLength = Math.min(insertionContextPrefixMatchLength, insertionContextSuffixMatchLength);

        int maxInsertionContext = Math.min(genomicInsertionFivePrime.length(), fivePrimeSeq.length()) + Math.min(genomicInsertionThreePrime.length(), threePrimeSeq.length());
        double insertionContextScore =  ((double) (insertionContextSuffixMatchLength + insertionContextPrefixMatchLength)) / maxInsertionContext;

        return new ImmutablePair<>(insertionContextScore, minInsertionContextMaxLength);
    }


    /**
     * Adjust 5' and 3' context of indel such that indel is fully shifted to 5'
     *
     * See: http://github.com/Illumina/GTCtoVCF/blob/develop/BPMRecord.py
     *
     * @param fivePrimeSeq 5' sequence
     * @param indelSeq indel sequence
     * @param threePrimeSeq 3' sequence

     * @return ImmutablePair containing the 5' and 3' sequences, left shifted
     *
     */
    static ImmutablePair<String, String> illuminaLeftShift(final String fivePrimeSeq, final String indelSeq, final String threePrimeSeq) {
        String modifiedFivePrimeSeq = fivePrimeSeq;
        StringBuilder modifiedThreePrimeSeq = new StringBuilder(threePrimeSeq);

        final int indelLength = indelSeq.length();

        while (modifiedFivePrimeSeq.endsWith(indelSeq)) {
            modifiedFivePrimeSeq = modifiedFivePrimeSeq.substring(0, modifiedFivePrimeSeq.length() - indelLength);
            modifiedThreePrimeSeq.insert(0, indelSeq);
        }
        // May have not fully shifted homopolymer
        while (StringUtils.repeat(modifiedFivePrimeSeq.charAt(modifiedFivePrimeSeq.length() - 1), indelLength).equals(indelSeq)) {
            modifiedThreePrimeSeq.insert(0, modifiedFivePrimeSeq.substring(modifiedFivePrimeSeq.length() - 1));
            modifiedFivePrimeSeq = modifiedFivePrimeSeq.substring(0, modifiedFivePrimeSeq.length() - 1);
        }
        return new ImmutablePair<>(modifiedFivePrimeSeq, modifiedThreePrimeSeq.toString());
    }

    /**
     * Find the number of common characters in the prefixes of two strings
     *
     * @param seq1 sequence
     * @param seq2 sequence

     * @return The number of matching characters
     */
    static int commonPrefixLength(final String seq1, final String seq2) {
        if (StringUtils.isEmpty(seq1) || (StringUtils.isEmpty(seq2))) {
            return 0;
        }
        int commonLength = 0;
        while ((commonLength < seq1.length()) && (commonLength < seq2.length()) && (seq1.charAt(commonLength) == seq2.charAt(commonLength))) {
            commonLength++;
        }
        return commonLength;
    }


    /**
     * Find the number of common characters in the suffixes of two strings
     *
     * @param seq1 sequence
     * @param seq2 sequence

     * @return The number of matching characters
     */
    static int commonSuffixLength(final String seq1, final String seq2) {
        if (StringUtils.isEmpty(seq1) || (StringUtils.isEmpty(seq2))) {
            return 0;
        }
        int index1 = seq1.length() - 1;
        int index2 = seq2.length() - 1;
        int commonLength = 0;
        while ((index1 >= 0) && (index2 >=0) && (seq1.charAt(index1--) == seq2.charAt(index2--))) {
            commonLength++;
        }
        return commonLength;
    }


    /**
     * Use regex to capture the insertion sequence and the sequence after the indel.
     * <p>
     * The source sequence is of the from V[W/X]Y, where V,W,X,Y are sequences.
     * <p>
     * - A SNP example looks like:    AGGGAGTC[A/G]GGTTGCGA
     * V     W X    Y
     * <p>
     * - A InDel example looks like:  AGCCTCGA[-/CGAA]TCACC
     * V     W   X   Y
     */
    private static Matcher parseSourceSeq(final String sourceSeq) {
        final Matcher matcher = pattern.matcher(sourceSeq);
        if (matcher.find()) {
            return matcher;
        } else {
            throw new PicardException("Could not find the pattern V[W/X]Y in the SourceSeq: " + sourceSeq);
        }
    }

    private static void validateThatSequenceOnlyContainsACGTCharacters(final String sequenceDescription, final String sequence) throws PicardException {
        if (!ACGT_PATTERN.matcher(sequence).find()) {
            throw new PicardException(sequenceDescription + " contains non-ACGT character(s)");
        }
    }

    /**
     * Determines the chromosome and position of the record on the target build
     */
    private void liftOverToTargetBuild(final Build37ExtendedIlluminaManifestRecord build37ExtendedIlluminaManifestRecord,
                                   final IlluminaManifestRecord illuminaManifestRecord) {

        final String supportedBuildNumber = illuminaManifestRecord.getMajorGenomeBuild();
        final File chainFileToTargetBuild = chainFilesMap.get(supportedBuildNumber);
        final LiftOver liftOver = new LiftOver(chainFileToTargetBuild);
        final Interval interval = new Interval(illuminaManifestRecord.getChr(), illuminaManifestRecord.getPosition(), illuminaManifestRecord.getPosition());
        final Interval targetBuildInterval = liftOver.liftOver(interval);

        if (targetBuildInterval != null) {
            build37ExtendedIlluminaManifestRecord.b37Chr = targetBuildInterval.getContig();
            build37ExtendedIlluminaManifestRecord.b37Pos = targetBuildInterval.getStart();

            // Validate that the reference allele at the lifted over coordinates matches that of the original.
            String originalRefAllele = getSequenceAt(referenceFilesMap.get(supportedBuildNumber), illuminaManifestRecord.getChr(), illuminaManifestRecord.getPosition(), illuminaManifestRecord.getPosition());
            String newRefAllele = getSequenceAt(referenceFilesMap.get(targetBuild), build37ExtendedIlluminaManifestRecord.b37Chr, build37ExtendedIlluminaManifestRecord.b37Pos, build37ExtendedIlluminaManifestRecord.b37Pos);
            if (originalRefAllele.equals(newRefAllele)) {
                log.debug("Lifted over record " + build37ExtendedIlluminaManifestRecord);
                log.debug(" From build " + supportedBuildNumber +
                        " chr=" + illuminaManifestRecord.getChr() +
                        ", position=" + illuminaManifestRecord.getPosition() +
                        " To build " + targetBuild +
                        " chr=" + build37ExtendedIlluminaManifestRecord.b37Chr + ", position=" + build37ExtendedIlluminaManifestRecord.b37Pos);
            } else {
                build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.LIFTOVER_FAILED;
                log.error("Liftover failed for record: " + build37ExtendedIlluminaManifestRecord);
                log.error( " Sequence at lifted over position does not match that at original position");
            }
        } else {
            build37ExtendedIlluminaManifestRecord.flag = Build37ExtendedIlluminaManifestRecord.Flag.LIFTOVER_FAILED;
            log.error("Liftover failed for record: " + build37ExtendedIlluminaManifestRecord);
        }
    }

    public Set<String> getUnsupportedBuilds() {
        return unsupportedBuilds;
    }

    public boolean isRefStrandDefinedInManifest() {
        return refStrandDefinedInManifest;
    }

    /**
     * Find the sequence in a reference sequence file.
     * for the contig in the range [start,stop]
     * @param refFile ReferenceSequenceFile to use
     * @param chr Contig whose subsequence to retrieve.
     * @param startPos inclusive, 1-based start of region.
     * @param endPos inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.
     */
    private static String getSequenceAt(final ReferenceSequenceFile refFile, final String chr, final int startPos, final int endPos) {
        final int contigLength = refFile.getSequenceDictionary().getSequence(chr).getSequenceLength();
        int usedEndPos = Math.min(endPos, contigLength);
        return new String(refFile.getSubsequenceAt(chr, startPos, usedEndPos).getBases()).toUpperCase();
    }
}
