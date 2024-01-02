package picard.flow;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.ArrayUtils;
import picard.PicardException;

import java.util.Arrays;
import java.util.List;
import java.util.PriorityQueue;

public class FlowBasedRead {

    public static final int MAX_CLASS = 12; //note - this is a historical value to support files with max class is not defined in the header, it is expected that mc tag exists in the CRAM

    // constants
    private final double MINIMAL_CALL_PROB = 0.1;
    // constants for clippingTagContains.
    // The tag is present when the end of the read was clipped at base calling.
    // The value of the tag is a string consisting of any one or more of the following:
    // A - adaptor clipped
    // Q - quality clipped
    // Z - "three zeros" clipped

    final public static String FLOW_MATRIX_TAG_NAME = "tp";
    final public static String FLOW_MATRIX_T0_TAG_NAME = "t0";
    final public static String FLOW_MATRiX_OLD_TAG_KR = "kr";
    final public static String FLOW_MATRiX_OLD_TAG_TI = "ti";

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
     * The maximal length of an hmer that can be encoded (normally in the 10-12 range)
     */
    private int maxHmer;
    private double totalMinErrorProbability;
    private double perHmerMinErrorProbability;
    /**
     * The order in which flow key in encoded (See decription for key field). Flow order may be wrapped if a longer one
     * needed.
     */
    private byte[] flowOrder;

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
     * The direction of this read. After being red, the direction will also swing to be REFERENCE
     */
    private Direction direction = Direction.SYNTHESIS;

    /**
     * The flow based argument collection under which this read was created
     */
    private final FlowBasedArgumentCollection fbargs;

    public enum Direction {
        REFERENCE, SYNTHESIS
    }


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

        implementMatrixMods(FlowBasedReadUtils.getFlowMatrixModsInstructions(fbargs.flowMatrixMods, maxHmer));

        //Spread boundary flow probabilities when the read is unclipped
        //in this case the value of the hmer is uncertain
        if (!fbargs.keepBoundaryFlows) {
            if (samRecord.getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP && samRecord.getCigar().getFirstCigarElement().getLength() > 0) {
                spreadFlowLengthProbsAcrossCountsAtFlow(findFirstNonZero(key));
            }
            if (samRecord.getCigar().getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP && samRecord.getCigar().getLastCigarElement().getLength() > 0) {
                spreadFlowLengthProbsAcrossCountsAtFlow(findLastNonZero(key));
            }
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
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction
        if (fbargs.fillingValue == 0) {
            totalMinErrorProbability = estimateFillingValue();
        } else {
            totalMinErrorProbability = fbargs.fillingValue;
        }
        perHmerMinErrorProbability = totalMinErrorProbability/getMaxHmer();

        key = FlowBasedReadUtils.baseArrayToKey(samRecord.getReadBases(), _flowOrder);
        flowOrder = FlowBasedReadUtils.getFlowToBase(_flowOrder, key.length);
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
        if ((t0 != null) && fbargs.useT0Tag) {
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
                parseZeroQuals(t0probs, i, qualOfs);
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
    private void parseZeroQuals(final double[] probs, final int flowIdx, final int qualOfs) {
        if ((qualOfs == 0) | (qualOfs == probs.length)) { // do not report zero error probability on the edge of the read
            return;
        }
        if ((Math.min(probs[qualOfs - 1], probs[qualOfs]) <= totalMinErrorProbability)) {
            flowMatrix[1][flowIdx] = Math.max(flowMatrix[1][flowIdx], perHmerMinErrorProbability);
        } else {
            flowMatrix[1][flowIdx] = Math.max(flowMatrix[1][flowIdx], Math.min(probs[qualOfs - 1], probs[qualOfs]));
        }
    }

    //Finds the quality that is being set when the probability of error is very low
    private double estimateFillingValue(){
        final byte[] quals = samRecord.getBaseQualities();
        final byte[] tp = samRecord.getSignedByteArrayAttribute(FLOW_MATRIX_TAG_NAME);
        byte maxQual = 0;

        for (int i = 0; i < quals.length; i++) {
            if (tp[i]!=0){
                continue;
            }
            if (quals[i] > maxQual){
                maxQual = quals[i];
            }
        }
        // in the very rare case when there is no tp=0 anywhere, just return the default "filling value"
        if (maxQual==0){
            return fbargs.fillingValue;
        }
        return Math.pow(10, -maxQual / 10);
    }

    public int getMaxHmer() {
        return maxHmer;
    }

    public int getNFlows() {
        return key.length;
    }

    public Direction getDirection() {
        return direction;
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

    // execute the matrix modifications
    private void implementMatrixMods(final int[] flowMatrixModsInstructions) {

        if (flowMatrixModsInstructions != null) {
            for (int hmer = 0; hmer < flowMatrixModsInstructions.length; hmer++) {
                final int hmer2 = flowMatrixModsInstructions[hmer];
                if (hmer2 != 0) {
                    for (int pos = 0; pos < flowMatrix[0].length; pos++) {

                        if (flowMatrix[hmer][pos] > flowMatrix[hmer2][pos]) {
                            flowMatrix[hmer2][pos] = flowMatrix[hmer][pos];
                        }

                        // if we are copying bacwards, zero out source
                        if (hmer > hmer2)
                            flowMatrix[hmer][pos] = 0;
                    }
                }
            }
        }
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

    public void setDirection(final Direction dir) {
        direction = dir;
    }

    public byte[] getFlowOrderArray() {
        return flowOrder;
    }

    public int getKeyLength() {
        return key.length;
    }

    public int[] getKey() {
        return key;
    }

    //functions that take care of simulating base format
    //they perform modifications on the flow matrix that are defined in applyFilteringFlowMatrix
    //this function was only applied when we tested what is the necessary information to be reported in the flow matrix
    private void applyFilteringFlowMatrix() {

        if (fbargs.disallowLargerProbs) {
            removeLargeProbs();
        }

        if (fbargs.removeLongerThanOneIndels) {
            removeLongIndels(key);
        }

        if (fbargs.removeOneToZeroProbs) {
            removeOneToZeroProbs(key);
        }

        if ((fbargs.lumpProbs)) {
            lumpProbs();
        }
        clipProbs();

        if (fbargs.symmetricIndels) {
            smoothIndels(key);
        }
        if (fbargs.onlyInsOrDel) {
            reportInsOrDel(key);
        }

        if ((fbargs.retainMaxNProbs)) {
            reportMaxNProbsHmer(key);
        }

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

    /**
     * remove probabilities larger than 1
     */
    private void removeLargeProbs() {
        for (int i = 0; i < getNFlows(); i++) {
            for (int j = 0; j < getMaxHmer() + 1; j++) {
                if (flowMatrix[j][i] > 1) {
                    flowMatrix[j][i] = 1;
                }
            }
        }
    }

    /**
     * This is vestigial and applies only to old formats
     *
     * @param key_kh
     */
    private void removeLongIndels(final int[] key_kh) {
        for (int i = 0; i < getNFlows(); i++) {
            for (int j = 0; j < getMaxHmer() + 1; j++) {
                if (Math.abs(j - key_kh[i]) > 1) {
                    flowMatrix[j][i] = perHmerMinErrorProbability;
                }
            }
        }
    }

    /**
     * This is vestigial and applies only to old formats
     *
     * @param key_kh
     */
    private void removeOneToZeroProbs(final int[] key_kh) {
        for (int i = 0; i < getNFlows(); i++) {
            if (key_kh[i] == 0) {
                for (int j = 1; j < getMaxHmer() + 1; j++) {
                    flowMatrix[j][i] = perHmerMinErrorProbability;
                }
            }
        }
    }

    /**
     * Smooth out probabilities by averaging with neighbours
     *
     * @param kr
     */
    private void smoothIndels(final int[] kr) {
        for (int i = 0; i < kr.length; i++) {
            final int idx = kr[i];
            if ((idx > 1) && (idx < maxHmer)) {
                final double prob = (flowMatrix[idx - 1][i] + flowMatrix[idx + 1][i]) / 2;
                flowMatrix[idx - 1][i] = prob;
                flowMatrix[idx + 1][i] = prob;
            }
        }
    }

    /**
     * Only report probability of insertions or of deletions, never both
     *
     * @param kr
     */
    private void reportInsOrDel(final int[] kr) {
        for (int i = 0; i < kr.length; i++) {
            final int idx = kr[i];
            if ((idx > 1) && (idx < maxHmer)) {
                if ((flowMatrix[idx - 1][i] > perHmerMinErrorProbability) && (flowMatrix[idx + 1][i] > perHmerMinErrorProbability)) {
                    final int fixCell = flowMatrix[idx - 1][i] > flowMatrix[idx + 1][i] ? idx + 1 : idx - 1;
                    flowMatrix[fixCell][i] = perHmerMinErrorProbability;
                }
            }
        }
    }

    /**
     * Combine all probabilities of insertions together and report them as probabilities of 1-mer insertion
     * Combine all probabilities of deletions together and report them as probabilities of 1-mer deletion
     */
    private void lumpProbs() {

        for (int i = 0; i < getMaxHmer(); i++) {
            for (int j = 0; j < getNFlows(); j++) {
                final int fkey = key[j];
                if (flowMatrix[i][j] <= perHmerMinErrorProbability) {
                    continue;
                } else {
                    if ((i - fkey) < -1) {
                        flowMatrix[fkey - 1][j] += flowMatrix[i][j];
                        flowMatrix[i][j] = perHmerMinErrorProbability;
                    } else if ((i - fkey) > 1) {
                        flowMatrix[fkey + 1][j] += flowMatrix[i][j];
                        flowMatrix[i][j] = perHmerMinErrorProbability;
                    }

                }

            }
        }

    }

    /*
    Given full vector of error probabilities retain only the probabilities that can be reported in the base format
    (N+1/2 highest error probabilities)
     */
    private void reportMaxNProbsHmer(final int[] key) {
        final double[] tmpContainer = new double[maxHmer];
        for (int i = 0; i < key.length; i++) {

            for (int j = 0; j < tmpContainer.length; j++) {
                tmpContainer[j] = flowMatrix[j][i];
            }
            final int k = (key[i] + 1) / 2;
            final double kth_highest = findKthLargest(tmpContainer, k + 1);
            for (int j = 0; j < maxHmer; j++) {
                if (flowMatrix[j][i] < kth_highest){
                    flowMatrix[j][i] = perHmerMinErrorProbability;
                }
            }
        }

    }


    private static double findKthLargest(final double[] nums, final int k) {
        final PriorityQueue<Double> q = new PriorityQueue<Double>(k);
        for (final double i : nums) {
            q.offer(i);

            if (q.size() > k) {
                q.poll();
            }
        }

        return q.peek();
    }
}


