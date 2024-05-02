package picard.flow;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.SequenceUtil;

public class FlowReadGroupInfo {
    private static final String PLATFORM_ULTIMA = "ULTIMA";
    private static final String PLATFORM_LS454 = "LS454";

    final public String  flowOrder;
    final public int     maxClass;
    final public boolean isFlowPlatform;
    private String  reversedFlowOrder = null;

    public FlowReadGroupInfo(final SAMReadGroupRecord readGroup) {
        isFlowPlatform = PLATFORM_ULTIMA.equals(readGroup.getPlatform()) || PLATFORM_LS454.equals(readGroup.getPlatform());

        if (isFlowPlatform) {
            final String mc = readGroup.getAttribute(FlowBasedRead.MAX_CLASS_READ_GROUP_TAG);
            this.maxClass = (mc == null) ? FlowBasedRead.MAX_CLASS : Integer.parseInt(mc);
        } else {
            this.maxClass = 0;
        }
        this.flowOrder = readGroup.getFlowOrder();

        if ( PLATFORM_ULTIMA.equals(readGroup.getPlatform()) && flowOrder == null ) {
            throw new RuntimeException("Malformed Ultima read group identified, aborting: " + readGroup);
        }
    }

    public synchronized String getReversedFlowOrder() {
        if ( reversedFlowOrder == null ) {
            reversedFlowOrder = SequenceUtil.reverseComplement(flowOrder);
        }
        return reversedFlowOrder;
    }
}
