package picard.flow;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.SimpleAllele;
import org.apache.commons.lang3.ArrayUtils;

/**
 * Haplotype that also keeps information on the flow space @see FlowBasedRead
 * Haplotype can't be extended, so this extends Allele
 */
public class FlowBasedHaplotype  extends SimpleAllele {
    private static final long serialVersionUID = 42L;
    private int [] key;
    private int [] rKey;
    private int[] flow2base;
    private int[] rFlow2Base;
    private Locatable genomeLoc;
    private Cigar cigar;
    private byte[] flowOrderArray;

    public FlowBasedHaplotype(final SimpleAllele sourceAllele, final String flowOrder) {
        this(sourceAllele, null, null, flowOrder);
    }

    public FlowBasedHaplotype(final SimpleAllele sourceAllele, final Locatable genomeLoc, final Cigar cigar, final String flowOrder) {
        super(sourceAllele.getBases(), sourceAllele.isReference());
        key = FlowBasedKeyCodec.baseArrayToKey(sourceAllele.getBases(), flowOrder);
        this.genomeLoc = genomeLoc;
        this.cigar = cigar;
        flow2base = FlowBasedKeyCodec.getKeyToBase(key);
        rKey = key.clone();
        ArrayUtils.reverse(rKey);
        rFlow2Base = FlowBasedKeyCodec.getKeyToBase(rKey);
        flowOrderArray = FlowBasedKeyCodec.getFlowToBase(flowOrder, key.length);
    }


    public int getKeyLength() {
        return key.length;
    }

    public int[] getKey() {
        return key;
    }

    public int getStart(){
        return genomeLoc.getStart();
    }

    public int getEnd() {
        return genomeLoc.getEnd();
    }

    public String getChr() {
        return genomeLoc.getContig();
    }

    public Cigar getCigar(){
        return cigar;
    }

    public int[] findLeftClipping(final int baseClipping) {
        return FlowBasedKeyCodec.findLeftClipping(baseClipping, flow2base, key);
    }

    public int[] findRightClipping(final int baseClipping) {
        return FlowBasedKeyCodec.findRightClipping(baseClipping, rFlow2Base, rKey);
    }

    public byte [] getFlowOrderArray() {
        return flowOrderArray;
    }
}
