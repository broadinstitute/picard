package net.sf.picard.metrics;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;

public class SAMRecordAndReference {
    private final SAMRecord samRec;
    private final ReferenceSequence refSeq;

    public SAMRecordAndReference(final SAMRecord samRec, final ReferenceSequence refSeq) {
        this.samRec = samRec;
        this.refSeq = refSeq;
    }

    public SAMRecord getSamRecord() {
        return samRec;
    }

    public ReferenceSequence getReferenceSequence() {
        return refSeq;
    }
}