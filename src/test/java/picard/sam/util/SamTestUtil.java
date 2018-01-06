package picard.sam.util;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;

/**
 * Created by farjoun on 9/21/17.
 */
public class SamTestUtil {
    public static long countSamTotalRecord(final File samFile) {
        final SamReader reader = SamReaderFactory.make().open(samFile);
        assert reader.hasIndex();
        long total = 0;

        for (int i = 0; i < reader.getFileHeader().getSequenceDictionary().size(); i++) {
            total += reader.indexing().getIndex().getMetaData(i).getAlignedRecordCount();
            total += reader.indexing().getIndex().getMetaData(i).getUnalignedRecordCount();
        }
        return total;
    }
}
