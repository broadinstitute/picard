package picard.flow;

import htsjdk.samtools.*;
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

    static public class ReadGroupInfo {
        private static final String PLATFORM_ULTIMA = "ULTIMA";
        private static final String PLATFORM_LS454 = "LS454";


        final public String  flowOrder;
        final public int     maxClass;
        final public boolean isFlowPlatform;

        public ReadGroupInfo(final SAMReadGroupRecord readGroup) {
            isFlowPlatform = PLATFORM_ULTIMA.equals(readGroup.getPlatform()) || PLATFORM_LS454.equals(readGroup.getPlatform());

            if (isFlowPlatform) {
                this.flowOrder = readGroup.getFlowOrder();
                final String mc = readGroup.getAttribute(FlowBasedRead.MAX_CLASS_READ_GROUP_TAG);
                this.maxClass = (mc == null) ? FlowBasedRead.MAX_CLASS : Integer.parseInt(mc);
            } else {
                this.flowOrder = null;
                this.maxClass = 0;
            }
        }

    }

    public static synchronized ReadGroupInfo getReadGroupInfo(final SAMFileHeader hdr, final SAMRecord read) {

        final String name = read.getReadGroup().getReadGroupId();
        ReadGroupInfo info = readGroupInfo.get(name);
        if ( info == null ) {
            readGroupInfo.put(name, info = new ReadGroupInfo(hdr.getReadGroup(name)));
        }
        return info;
    }

}
