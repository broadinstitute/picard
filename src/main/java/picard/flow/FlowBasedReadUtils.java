package picard.flow;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;

import java.util.Arrays;
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

    public static final int FLOW_SUM_OF_BASE_QUALITY_THRESHOLD = 15;
    public static final FlowBasedArgumentCollection DEFAULT_FLOW_BASED_ARGUMENT_COLLECTION = new FlowBasedArgumentCollection();
    static final public int FLOW_BASED_INSIGNIFICANT_END = 0;

    private static final Map<String, ReadGroupInfo> readGroupInfo = new LinkedHashMap<>();

    public enum CycleSkipStatus {
        NS(0),         // no skip
        PCS(1),        // possible cycle skip
        CS(2);         // cycle skip

        private int priority;

        CycleSkipStatus(final int priority) {
            this.priority = priority;
        }

        int getPriority() {
            return this.priority;
        }

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

    public static boolean readEndMarkedUncertain(final SAMRecord rec) {
        final String        tm = rec.getStringAttribute(FlowBasedRead.CLIPPING_TAG_NAME);
        if ( tm == null ) {
            return false;
        } else {
            return tm.indexOf('Q') >= 0 || tm.indexOf('Z') >= 0;
        }
    }

    public static boolean readEndMarkedUnclipped(final SAMRecord rec, boolean FLOW_Q_IS_KNOWN_END) {
        final String        tm = rec.getStringAttribute(FlowBasedRead.CLIPPING_TAG_NAME);
        if ( tm == null ) {
            return false;
        } else {
            return (tm.indexOf('A') >= 0) || (FLOW_Q_IS_KNOWN_END && (tm.indexOf('Q') >= 0));
        }
    }

    // get flow order for a specific read
    public static byte[] getReadFlowOrder(final SAMFileHeader header, SAMRecord read) {

        // are we looking for a specific read group, as specified by the read?
        {
            final SAMReadGroupRecord rg = (read != null) ? read.getReadGroup() : null;
            if (rg != null && rg.getFlowOrder() != null) {
                return rg.getFlowOrder().getBytes();
            }
        }

        // if here, either no read was specified, or the read has no group, or the group is not found, or it has no flow
        // revert to old behavior of returning the first found
        for ( SAMReadGroupRecord rg : header.getReadGroups() ) {
            // must match read group name?
            String      flowOrder = rg.getFlowOrder();
            if ( flowOrder != null ) {
                return flowOrder.getBytes();
            }
        }
        return null;
    }

    /**
     * Computes the sum of base qualities of the given flow read.
     */
    public static int flowSumOfBaseQualities(final SAMRecord read) {
        if (read == null) {
            return 0;
        } else {
            int sum = 0;

            // access qualities and bases
            byte[]      quals = read.getBaseQualities();
            byte[]      bases = read.getReadBases();

            // loop on bases, extract qual related to homopolymer from start of homopolymer
            int         i = 0;
            byte        lastBase = 0;
            byte        effectiveQual = 0;
            for (final byte base : bases ) {
                if ( base != lastBase )
                    effectiveQual = quals[i];
                if ( effectiveQual >= FLOW_SUM_OF_BASE_QUALITY_THRESHOLD )
                    sum += effectiveQual;
                lastBase = base;
                i++;
            }

            return sum;
        }
    }

    public static boolean hasFlowTags(final SAMRecord rec) {
        return rec.hasAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_KR)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_TI);

    }

    /**
     *
     * This is the function to run if you want to ask if the data are flow-based
     *
     * @param hdr - file header
     * @param read - the read
     * @return true if the read is flow-based
     */
    public static boolean isFlowPlatform(final SAMFileHeader hdr, final SAMRecord read) {
        if (!hasFlowTags(read)){
            return false;
        }
        return getReadGroupInfo(hdr, read).isFlowPlatform;
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
     * Finds a usable FlowOrder to be used for engine calculation (when no specufic flow order already established for a specific read)
     */
    public static String findFirstUsableFlowOrder(final SAMFileHeader hdr, final FlowBasedArgumentCollection fbargs) {
        for ( final SAMReadGroupRecord rg : hdr.getReadGroups() ) {
            final String flowOrder = rg.getFlowOrder();
            if ( flowOrder != null && flowOrder.length() >= fbargs.flowOrderCycleLength ) {
                return flowOrder.substring(0, fbargs.flowOrderCycleLength);
            }
        }

        throw new PicardException("Unable to perform flow based operations without the flow order");
    }

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
     * create a FlowBasedRead from a proper SAMRecord
     */
    static public FlowBasedRead convertToFlowBasedRead(SAMRecord read, SAMFileHeader header) {
        final ReadGroupInfo readGroupInfo = getReadGroupInfo(header, read);
        return new FlowBasedRead(read, readGroupInfo.flowOrder, readGroupInfo.maxClass, DEFAULT_FLOW_BASED_ARGUMENT_COLLECTION);
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
}
