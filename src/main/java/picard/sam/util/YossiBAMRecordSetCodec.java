package picard.sam.util;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SortingCollection;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class YossiBAMRecordSetCodec implements SortingCollection.Codec<YossiBAMRecordSet> {
    private final BAMRecordCodec bamRecordCodec;
    private InputStream is;
    private OutputStream os;
    private SAMFileHeader header;

    public YossiBAMRecordSetCodec(final SAMFileHeader header) {
        this.header = header;
        this.bamRecordCodec = new BAMRecordCodec(header);
    }

    @Override
    public void setOutputStream(OutputStream os) {
        this.os = os;
        bamRecordCodec.setOutputStream(os);
    }

    @Override
    public void setInputStream(InputStream is) {
        this.is = is;
        bamRecordCodec.setInputStream(is);
    }

    @Override
    public void encode(YossiBAMRecordSet val) {
        try {
            os.write(val.size());
            os.flush();

            for(SAMRecord rec: val) {
                bamRecordCodec.encode(rec);
            }
            os.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public YossiBAMRecordSet decode() {
        try {
            int numRecords = is.read();
            if(numRecords == -1) {
                return null;
            }
            YossiBAMRecordSet recs = new YossiBAMRecordSet(numRecords);
            for(int i = 0; i < numRecords; ++i) {
                recs.add(bamRecordCodec.decode());
            }
            return recs;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public SortingCollection.Codec<YossiBAMRecordSet> clone() {
        return new YossiBAMRecordSetCodec(this.header);
    }
}
