package net.sf.picard.metrics;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;

public abstract class SAMRecordAndReferenceMultiLevelCollector<BEAN extends MetricBase, HKEY extends Comparable> extends MultiLevelCollector<BEAN, HKEY, SAMRecordAndReference> {

        @Override
        protected SAMRecordAndReference makeArg(SAMRecord samRec, final ReferenceSequence refSeq) {
            return new SAMRecordAndReference(samRec, refSeq);
        }
}


