package picard.flow;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * Utility class for working with flow-based reads
 *
 * The main member static class is {@code ReadGroupInfo} that contains methods that allow
 * working with headers of flow based reads and extracting the flow order and the maximal hmer class called.
 * It also contains methods to check how the read was clipped ({@code readEndMarkedUnclipped}, {@code readEndMarkedUncertain})
 * Lastly, {@code FlowBasedReadUtils.isFlowPlatform} is **the** function to determine if the data are flow-based
 *
 */
public class FlowBasedReadUtils {
    
    private static final Map<String, ReadGroupInfo> readGroupInfo = new LinkedHashMap<>();

    /**
     * Converts base space sequence to flow space
     * @param bases base space sequence
     * @param flowOrder flow order
     * @return Array of flow values
     */

    static public int[] baseArrayToKey(final byte[] bases, final String flowOrder){

        final ArrayList<Integer> result = new ArrayList<>();
        final byte[] flowOrderBytes = flowOrder.getBytes();
        int loc = 0;
        int flowNumber = 0 ;
        final int period = flowOrderBytes.length;
        int periodGuard = 0;
        while ( loc < bases.length ) {
            final byte flowBase = flowOrderBytes[flowNumber%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!= 'N')) {
                result.add(0);
                if ( ++periodGuard > period )
                    throw new PicardException("baseArrayToKey periodGuard tripped, on " + new String(bases) + ", flowOrder: " + flowOrder
                    + " This probably indicates the presence of a base (value) in the sequence that is not included in the provided flow order");
            } else {
                int count = 0;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]== 'N')) ){
                    loc++;
                    count ++;
                }
                result.add(count);
                periodGuard = 0;
            }
            flowNumber++;
        }
        final int[] ret = new int[result.size()];
        for (int i = 0; i < result.size(); i++) {
            ret[i] = result.get(i);
        }
        return ret;
    }

    /**
     * For every flow of the key output the index of the last base that was output prior to this flow
     * @param key given key
     * @return array
     */
    static public int[] getKeyToBase(final int[] key) {
        final int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1; i < result.length; i++) {
            result[i] = result[i - 1] + key[i - 1];
        }
        return result;

    }

    /**
     * For every flow of the key output the nucleotide that is being read for this flow
     * @param flowOrder given flow order
     * @param expectedLength the length of the key (key is not provided)
     * @return array of bases
     */

    static public byte[] getFlowToBase(final String flowOrder, final int expectedLength) {
        final byte[] result = new byte[expectedLength] ;
        for ( int i = 0; i < result.length; i++ ) {
            result[i] = (byte)flowOrder.charAt(i%flowOrder.length());
        }
        return result;
    }

    /**
     * Prints the key as character-encoded string
     * @param ints (key array)
     * @return encoded string
     */
    static public String keyAsString(final int[] ints)
    {
        final StringBuilder   sb = new StringBuilder();

        for ( final int i : ints )
            sb.append((char)((i < 10) ? ('0' + i) : ('A' + i - 10)));

        return sb.toString();
    }

    /**
     * Converts a numerical array into flow space by flattening per-base elements to fill empty flows based on minimum scores.
     * This is intended to make it easier to transform per-base read statistics that are computed on the read in base-space
     * sanely into flow-space reads in a reasonable fashion that makes sense for parameters like indel-likelihoods.
     *
     * The rules are as follows:
     *  - For every run of homopolymers, the minimum score for any given base is translated to cover the whole run
     *  - For every 0 base flow, the score from the previous filled flow is copied into place.
     *  - For the beginning of the array (in the case of preceding 0-flows) the given default value is used.
     *
     * Example:
     *  A read with the following bases with a base-space calculated into flow space (ACTG) as:
     *     TTTATGC -> 0030101101
     *  With an input array like this:
     *     byte[]{1,2,3,4,5,6,7}
     *  We should expect the following result array (Which matches the key array in length)
     *     byte[]{defaultQual,defaultQual,1,1,4,4,5,6,6,7}
     *
     * @param bases base space sequence
     * @param keyLength size of the converted bases array in flow-space
     * @param baseSpacedArrayToConvert Array of base-space scores to be conformed to flow space
     * @param defaultQual default quality to use at the head of the array for non-covered flows
     * @param flowOrder flow order
     * @return Array of translated bases comparable to the keyLength
     */
    @VisibleForTesting
    static public byte[] baseArrayToKeySpace(final byte[] bases, final int keyLength, final byte[] baseSpacedArrayToConvert, final byte defaultQual, final String flowOrder){

        if (bases.length != baseSpacedArrayToConvert.length) {
            throw new IllegalArgumentException("Read and qual arrays do not match");
        }

        final byte[] result = new byte[keyLength];
        final byte[] flowOrderBytes = flowOrder.getBytes();
        int loc = 0;
        int flowNumber = 0 ;
        byte lastQual = defaultQual;
        final int period = flowOrderBytes.length;
        while ( loc < bases.length ) {
            final byte flowBase = flowOrderBytes[flowNumber%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!='N')) {
                result[flowNumber] = lastQual;
            } else {
                byte qual = Byte.MAX_VALUE;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]== 'N')) ){
                    qual = (byte)Math.min(baseSpacedArrayToConvert[loc], qual);
                    loc++;
                }
                result[flowNumber] = qual;
                lastQual = qual;

            }
            flowNumber++;
        }
        return result;
    }

    static public class ReadGroupInfo {
        final public String  flowOrder;
        final public int     maxClass;
        final public boolean isFlowPlatform;
        private String  reversedFlowOrder = null;

        public ReadGroupInfo(final SAMReadGroupRecord readGroup) {
            if (readGroup.getPlatform()==null){
                isFlowPlatform = false;
            } else if (NGSPlatform.fromReadGroupPL(readGroup.getPlatform())==NGSPlatform.UNKNOWN){
                isFlowPlatform = false;
            } else if (NGSPlatform.fromReadGroupPL(readGroup.getPlatform()) == NGSPlatform.LS454) {
                //old Ultima data can have PLATFORM==LS454 and not have FO tag
                isFlowPlatform = true;
            } else if (NGSPlatform.fromReadGroupPL(readGroup.getPlatform()) == NGSPlatform.ULTIMA){
                if (readGroup.getFlowOrder()!=null) {
                    isFlowPlatform = true;
                } else {
                    throw new RuntimeException("Malformed Ultima read group identified, aborting: " + readGroup);
                }
            } else {
                isFlowPlatform = false;
            }

            if (isFlowPlatform) {
                this.flowOrder = readGroup.getFlowOrder();
                String mc = readGroup.getAttribute(FlowBasedRead.MAX_CLASS_READ_GROUP_TAG);
                this.maxClass = (mc == null) ? FlowBasedRead.MAX_CLASS : Integer.parseInt(mc);
            } else { // not a flow platform
                this.flowOrder = null;
                this.maxClass = 0;
            }
        }

        public synchronized String getReversedFlowOrder() {
            if ( reversedFlowOrder == null ) {
                reversedFlowOrder = SequenceUtil.reverseComplement(flowOrder);
            }
            return reversedFlowOrder;
        }

    }

    public static boolean hasFlowTags(final SAMRecord rec) {
        return rec.hasAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_KR)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_TI);

    }

    public static synchronized ReadGroupInfo getReadGroupInfo(final SAMFileHeader hdr, final SAMRecord read) {

        if ( !hasFlowTags(read) ) {
            throw new IllegalArgumentException("read must be flow based: " + read);
        }

        String name = read.getReadGroup().getReadGroupId();
        ReadGroupInfo info = readGroupInfo.get(name);
        if ( info == null ) {
            readGroupInfo.put(name, info = new ReadGroupInfo(hdr.getReadGroup(name)));
        }
        return info;
    }

    /**
    /*
     * clips flows from the left to clip the input number of bases
     * Needed to trim the haplotype to the read
     * Returns number of flows to remove and the change in the left most remaining flow if necessary
     */
    static public int[] findLeftClipping(final int baseClipping, final int[] flow2base, final int[] key) {
        final int [] result = new int[2];
        if (baseClipping == 0 ){
            return result;
        }

        int stopClip = 0;
        for (int i = 0 ; i < flow2base.length; i++ ) {

            if (flow2base[i] + key[i] >= baseClipping) {
                stopClip = i;
                break;
            }
        }
        final int hmerClipped = baseClipping - flow2base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    /*
     * clips flows from the right to trim the input number of bases
     * Returns number of flows to remove and the change in the right most flow.
     */
    static public int[] findRightClipping(final int baseClipping, final int[] rFlow2Base, final int[] rKey) {
        final int [] result = new int[2];
        if (baseClipping == 0 ){
            return result;
        }

        int stopClip = 0;

        for (int i = 0; i < rFlow2Base.length; i++ ) {
            if (rFlow2Base[i] + rKey[i] >= baseClipping) {
                stopClip = i;
                break;
            }
        }

        final int hmerClipped = baseClipping - rFlow2Base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    /**
     * Retrieve flow matrix modifications matrix from its string argument format. This matrix contains
     * logic for modifying the flow matrix as it is read in. If the value of [n] is not zero,
     * then the hmer probability for hmer length n will be copied to the [n] position
     *
     * For the implementation logic, see FlowBasedRead.fillFlowMatrix
     */
    static public int[] getFlowMatrixModsInstructions(final String flowMatrixMods, final int maxHmer) {

        if ( flowMatrixMods != null ) {
            final int[] flowMatrixModsInstructions = new int[maxHmer + 1];

            final String[]    toks = flowMatrixMods.split(",");
            for ( int i = 0 ; i < toks.length - 1 ; i += 2 ) {
                final int hmer = Integer.parseInt(toks[i]);
                flowMatrixModsInstructions[hmer] = Integer.parseInt(toks[i + 1]);
            }
            return flowMatrixModsInstructions;
        } else {
            return null;
        }
    }

    /**
     * A canonical, master list of the standard NGS platforms.  These values
     * can be obtained (efficiently) from a SAMRecord object with the
     * getNGSPlatform method.
     */
    public enum NGSPlatform {
        // note the order of elements here determines the order of matching operations, and therefore the
        // efficiency of getting a NGSPlatform from a string.
        ILLUMINA(SequencerFlowClass.DISCRETE, "ILLUMINA", "SLX", "SOLEXA"),
        SOLID(SequencerFlowClass.DISCRETE, "SOLID"),
        LS454(SequencerFlowClass.FLOW, "454", "LS454"),
        COMPLETE_GENOMICS(SequencerFlowClass.DISCRETE, "COMPLETE"),
        PACBIO(SequencerFlowClass.DISCRETE, "PACBIO"),
        ION_TORRENT(SequencerFlowClass.FLOW, "IONTORRENT"),
        CAPILLARY(SequencerFlowClass.OTHER, "CAPILLARY"),
        HELICOS(SequencerFlowClass.OTHER, "HELICOS"),
        ULTIMA(SequencerFlowClass.FLOW, "ULTIMA"),
        UNKNOWN(SequencerFlowClass.OTHER, "UNKNOWN");


        /**
         * Array of the prefix names in a BAM file for each of the platforms.
         */
        protected final String[] BAM_PL_NAMES;
        protected final SequencerFlowClass sequencerType;

        NGSPlatform(final SequencerFlowClass type, final String... BAM_PL_NAMES) {
            if ( BAM_PL_NAMES.length == 0 ) throw new IllegalStateException("Platforms must have at least one name");

            for ( int i = 0; i < BAM_PL_NAMES.length; i++ )
                BAM_PL_NAMES[i] = BAM_PL_NAMES[i].toUpperCase();

            this.BAM_PL_NAMES = BAM_PL_NAMES;
            this.sequencerType = type;
        }

        /**
         * Returns the NGSPlatform corresponding to the PL tag in the read group
         * @param plFromRG -- the PL field (or equivalent) in a ReadGroup object.  Can be null => UNKNOWN
         * @return an NGSPlatform object matching the PL field of the header, or UNKNOWN if there was no match or plFromRG is null
         */
        public static NGSPlatform fromReadGroupPL(final String plFromRG) {
            if ( plFromRG == null ) return UNKNOWN;

            final String pl = plFromRG.toUpperCase();
            for ( final NGSPlatform ngsPlatform : NGSPlatform.values() ) {
                for ( final String bamPLName : ngsPlatform.BAM_PL_NAMES ) {
                    if ( pl.contains(bamPLName) )
                        return ngsPlatform;
                }
            }

            return UNKNOWN;
        }

        /**
         * In broad terms, each sequencing platform can be classified by whether it flows nucleotides in some order
         * such that homopolymers get sequenced in a single event (ie 454 or Ion) or it reads each position in the
         * sequence one at a time, regardless of base composition (Illumina or Solid).  This information is primarily
         * useful in the BQSR process
         */
        public enum SequencerFlowClass {
            DISCRETE,
            FLOW,
            OTHER //Catch-all for unknown platforms, as well as relics that GATK doesn't handle well (Capillary, Helicos)
        }
    }
}
