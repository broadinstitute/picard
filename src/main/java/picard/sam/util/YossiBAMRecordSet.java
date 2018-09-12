package picard.sam.util;

import htsjdk.samtools.SAMRecord;

import java.util.HashSet;

public class YossiBAMRecordSet extends HashSet<SAMRecord> {
    private SAMRecord primary1 = null;

    public YossiBAMRecordSet(int initialSize) {
        super(initialSize);
    }

    @Override
    public boolean add(SAMRecord samRecord) {
        return super.add(samRecord);
    }

    private void verifyPrimary1() {
        if(primary1 != null) {
            return;
        }
        SAMRecord unmapped1 = null;

        for(SAMRecord rec: this) {
            if(rec.getFirstOfPairFlag()) {
                if(!rec.isSecondaryOrSupplementary()) {
                    primary1 = rec;
                    return;
                }
                if(rec.getReadUnmappedFlag()) {
                    unmapped1 = rec;
                }
            }
        }
        //we couldn't find a primary, first-of-pair record, so it must be unmapped. return that instead
        primary1 = unmapped1;
    }

    public String getChr() {
        verifyPrimary1();
        return primary1.getReferenceName();
    }

    public int getPos() {
        verifyPrimary1();
        return primary1.getAlignmentStart();
    }
}
