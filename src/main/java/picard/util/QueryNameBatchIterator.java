package picard.util;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekableIterator;
import picard.sam.util.YossiBAMRecordSet;

public class QueryNameBatchIterator implements CloseableIterator<YossiBAMRecordSet> {

    final PeekableIterator<SAMRecord> peekit;

    public QueryNameBatchIterator(SAMRecordIterator iter) {
        this.peekit = new PeekableIterator<SAMRecord>(iter);
    }

    @Override
    public void close() {
        peekit.close();
    }

    @Override
    public boolean hasNext() {
        return peekit.hasNext();
    }

    @Override
    public YossiBAMRecordSet next() {
        String qname = peekit.peek().getReadName();
        YossiBAMRecordSet records = new YossiBAMRecordSet(0);

        while (peekit.peek() != null && peekit.peek().getReadName().equals(qname)) {
            records.add(peekit.next());
        }

        return records;
    }
}
