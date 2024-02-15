package picard.sam.util;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.utils.ValidationUtils;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.nio.file.Path;

/**
 * Created by farjoun on 9/21/17.
 */
public class SamTestUtil {
    public static long countSamTotalRecord(final Path samFile) {
        final SamReader reader = SamReaderFactory.make().open(samFile);
        ValidationUtils.validateArg(reader.hasIndex(),"Reader should have an index.");
        long total = 0;

        for (int i = 0; i < reader.getFileHeader().getSequenceDictionary().size(); i++) {
            total += reader.indexing().getIndex().getMetaData(i).getAlignedRecordCount();
            total += reader.indexing().getIndex().getMetaData(i).getUnalignedRecordCount();
        }
        return total;
    }
}
