package picard.flow;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import picard.PicardException;

public class FlowBasedRead {

    public static final int MAX_CLASS = 12; //note - this is a historical value to support files with max class is not defined in the header, it is expected that mc tag exists in the CRAM

    // constants
    private final double MINIMAL_CALL_PROB = 0.1;

    final private static String FLOW_MATRIX_TAG_NAME = "tp";
    final private static String FLOW_MATRIX_T0_TAG_NAME = "t0";
    final public static String MAX_CLASS_READ_GROUP_TAG = "mc";

    /**
     * The sam record from which this flow based read originated
     */
    private SAMRecord samRecord;

    /**
     * The flow key for the read - i.e. lengths of hmers in an flow order.
     * <p>
     * For example, assuming a flow order of TGCA, and a forward sequence of GGAAT, the key will be 0,2,0,2,1
     */
    private int[] key;

    /**
     * The maximal length of an hmer that can be encoded (normally in the 10-20 range)
     */
    private int maxHmer;

    /**
     * computed minimal error probability for an hmer
     */
    private double perHmerMinErrorProbability;

    /**
     * The probability matrix for this read. [n][m] position represents that probablity that an hmer of n length will be
     * present at the m key position. Therefore, the first dimention is in the maxHmer order, where the second dimension
     * is length(key).
     */
    private double[][] flowMatrix;

    /**
     * The validity status of the key. Certain operations may produce undefined/errornous results. This is signaled by
     * the read being marked with a validKey == false
     */
    private boolean validKey;

    /**
     * The flow based argument collection under which this read was created
     */
    private final FlowBasedArgumentCollection fbargs;


    /**
     * Same as above but constructs from SAMRecord
     *
     * @param samRecord record from SAM file
     * @param flowOrder flow order (single cycle)
     * @param maxHmer   maximal hmer to keep in the flow matrix
     * @param fbargs    arguments that control resoltion of the flow matrix
     */
    public FlowBasedRead(final SAMRecord samRecord, final String flowOrder, final int maxHmer, final FlowBasedArgumentCollection fbargs) {
        this.fbargs = fbargs;
        this.maxHmer = maxHmer;
        this.samRecord = samRecord;

        // read flow matrix in.
        if (samRecord.hasAttribute(FLOW_MATRIX_TAG_NAME)) {
            readFlowMatrix(flowOrder);
        } else {
            throw new PicardException("read missing flow matrix attribute: " + FLOW_MATRIX_TAG_NAME);
        }

        //Spread boundary flow probabilities when the read is unclipped
        //in this case the value of the hmer is uncertain
        if (samRecord.getReadUnmappedFlag() || (samRecord.getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP && samRecord.getCigar().getFirstCigarElement().getLength() > 0)) {
            spreadFlowLengthProbsAcrossCountsAtFlow(findFirstNonZero(key));
        }
        if (samRecord.getReadUnmappedFlag() || (samRecord.getCigar().getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP && samRecord.getCigar().getLastCigarElement().getLength() > 0)) {
            spreadFlowLengthProbsAcrossCountsAtFlow(findLastNonZero(key));
        }


        validateSequence();
    }

    //since the last unclipped flow is uncertain (we give high probabilities to
    //also hmers higher than the called hmer)
    private void spreadFlowLengthProbsAcrossCountsAtFlow(final int flowToSpread) {
        if (flowToSpread < 0) //boundary case when all the key is zero
            return;

        final int call = key[flowToSpread];
        if (call == 0) {
            throw new IllegalStateException("Boundary key value should not be zero for the spreading");
        }

        final int numberToFill = maxHmer - call + 1;
        double total = 0;
        for (int i = call; i < maxHmer + 1; i++)
            total += flowMatrix[i][flowToSpread];
        final double fillProb = Math.max(total / numberToFill, perHmerMinErrorProbability);
        for (int i = call; i < maxHmer + 1; i++) {
            flowMatrix[i][flowToSpread] = fillProb;
        }
    }

    // This is the code for parsing the current/production BAM format (with TP tag)
    private void readFlowMatrix(final String _flowOrder) {

        // generate key (base to flow space)

        // establish min probability
        final double totalMinErrorProbability = (fbargs.fillingValue == 0) ? estimateFillingValue() : fbargs.fillingValue;
        perHmerMinErrorProbability = totalMinErrorProbability / getMaxHmer();

        // generate flow key
        key = FlowBasedKeyCodec.baseArrayToKey(samRecord.getReadBases(), _flowOrder);

        // initialize matrix
        flowMatrix = new double[maxHmer + 1][key.length];
        for (int i = 0; i < maxHmer + 1; i++) {
            for (int j = 0; j < key.length; j++) {
                flowMatrix[i][j] = perHmerMinErrorProbability;
            }
        }

        // access qual, convert to flow representation
        final byte[] quals = samRecord.getBaseQualities();
        final byte[] tp = samRecord.getSignedByteArrayAttribute(FLOW_MATRIX_TAG_NAME);
        boolean specialTreatmentForZeroCalls = false;
        final byte[] t0 = SAMUtils.fastqToPhred(samRecord.getStringAttribute(FLOW_MATRIX_T0_TAG_NAME));
        final double[] t0probs = new double[quals.length];
        if ((t0 != null) && !fbargs.ignoreT0Tag) {
            specialTreatmentForZeroCalls = true;

            if (t0.length != tp.length) {
                throw new PicardException("Illegal read len(t0)!=len(qual): " + samRecord.getReadName());
            }
        }

        final double[] probs = new double[quals.length];
        for (int i = 0; i < quals.length; i++) {
            final double q = quals[i];
            final double p = Math.pow(10, -q / 10);
            probs[i] = p;
            if (specialTreatmentForZeroCalls) {
                final double qq = t0[i];
                t0probs[i] = Math.pow(10, -qq / 10);
            }
        }

        // apply key and qual/tp to matrix
        int qualOfs = 0; //converts between base -> flow
        for (int i = 0; i < key.length; i++) {
            final int run = key[i];
            if (run > 0) {
                parseSingleHmer(probs, tp, i, run, qualOfs);
            }

            if ((run == 0) && (specialTreatmentForZeroCalls)) {
                parseZeroQuals(t0probs, i, qualOfs, totalMinErrorProbability);
            }

            double totalErrorProb = 0;

            for (int k = 0; k < maxHmer; k++) {
                totalErrorProb += flowMatrix[k][i];
            }
            final double callProb = Math.max(MINIMAL_CALL_PROB, 1 - totalErrorProb);
            // the probability in the recalibration is not divided by two for hmers of length 1
            flowMatrix[Math.min(run, maxHmer)][i] = callProb;
            qualOfs += run;
        }
        applyFilteringFlowMatrix();
    }


    //convert qualities from the single hmer to a column in a flow matrix
    private void parseSingleHmer(final double[] probs, final byte[] tp, final int flowIdx,
                                 final int flowCall, final int qualOfs) {
        for (int i = qualOfs; i < qualOfs + flowCall; i++) {
            if (tp[i] != 0) {
                final int loc = Math.max(Math.min(flowCall + tp[i], maxHmer), 0);
                if (flowMatrix[loc][flowIdx] == perHmerMinErrorProbability) {
                    flowMatrix[loc][flowIdx] = probs[i];
                } else {
                    flowMatrix[loc][flowIdx] += probs[i];
                }
            }
        }
    }

    // convert qualities from the t0 tag to the probabilities of 1->0 error.
    // This function deals with t0 tag that encodes the probability of 1->0 error
    // in this case there is no nucleotide to place the error probability on, so we
    // place it on the neighboring bases and choose the **lower** error probability between the
    // neighbors (that's how T0 encoding works). The error is placed only on the 1-mer error assuming
    // that 2->0 errors are negligibly rare.
    private void parseZeroQuals(final double[] probs, final int flowIdx, final int qualOfs, final double totalMinErrorProbability) {
        if ((qualOfs == 0) | (qualOfs == probs.length)) { // do not report zero error probability on the edge of the read
            return;
        }
        if ((Math.min(probs[qualOfs - 1], probs[qualOfs]) <= totalMinErrorProbability)) {
            flowMatrix[1][flowIdx] = Math.max(flowMatrix[1][flowIdx], perHmerMinErrorProbability);
        } else {
            flowMatrix[1][flowIdx] = Math.max(flowMatrix[1][flowIdx], Math.min(probs[qualOfs - 1], probs[qualOfs]));
        }
    }

    //Finds the quality that is being set when the probability of error is the minimal
    private double estimateFillingValue(){
        final byte[] quals = samRecord.getBaseQualities();
        double maxQual = 0;

        for (int i = 0; i < quals.length; i++) {
            if (quals[i] > maxQual){
                maxQual = quals[i];
            }
        }

        // in the very rare case when there are no qualities in the read, we set the maxQual to 40
        if (maxQual==0){
            maxQual = 40;
        }

        return Math.pow(10, -maxQual / 10);
    }

    public int getMaxHmer() {
        return maxHmer;
    }

    public int getNFlows() {
        return key.length;
    }

    private void validateSequence() {
        for (final int b : key) {
            if (b > maxHmer - 1) {
                validKey = false;
            }
        }
        validKey = true;
    }

    public boolean isValid() {
        return validKey;
    }

    /**
     * get a specific cell from the flow matrix. Each cell contains the probability
     * for an hmer of the given length to appear the given position in the flow key
     *
     * @param flow - position in the flow key (index into key[])
     * @param hmer - length of the hmer
     * @return
     */
    public double getProb(final int flow, final int hmer) {
        double prob = flowMatrix[hmer < maxHmer ? hmer : maxHmer][flow];
        return (prob <= 1) ? prob : 1;
    }

    private static int findFirstNonZero(final int[] array) {
        int result = -1;
        for (int i = 0; i < array.length; i++) {
            if (array[i] != 0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private static int findLastNonZero(final int[] array) {
        int result = -1;
        for (int i = array.length - 1; i >= 0; i--) {
            if (array[i] != 0) {
                result = i;
                break;
            }
        }
        return result;
    }

    public int[] getKey() {
        return key;
    }

    //functions that take care of simulating base format
    //they perform modifications on the flow matrix that are defined in applyFilteringFlowMatrix
    //this function was only applied when we tested what is the necessary information to be reported in the flow matrix
    private void applyFilteringFlowMatrix() {
        clipProbs();
    }

    /**
     * clip probability values to fbargs.probabilityRatioThreshold
     */
    private void clipProbs() {
        for (int i = 0; i < getMaxHmer(); i++) {
            for (int j = 0; j < getNFlows(); j++) {
                if ((flowMatrix[i][j] <= perHmerMinErrorProbability) &&
                        (key[j] != i)) {
                    flowMatrix[i][j] = perHmerMinErrorProbability;
                }
            }
        }
    }

    public String getReadName() {
        return samRecord.getReadName();
    }
}


