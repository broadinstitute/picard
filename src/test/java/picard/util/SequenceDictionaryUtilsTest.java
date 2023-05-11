package picard.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.List;

public class SequenceDictionaryUtilsTest {

    @Test
    public void testSequenceDictionaryPositive() {
        final SAMSequenceDictionary s1 = new SAMSequenceDictionary(
                List.of(
                        new SAMSequenceRecord("chr1", 10),
                        new SAMSequenceRecord("chr2", 15),
                        new SAMSequenceRecord("chr3", 20)
                )
        );
        SequenceDictionaryUtils.assertSequenceDictionariesEqual(
                s1,
                "s1",
                new SAMSequenceDictionary(s1.getSequences()),
                "s2");
    }

    @DataProvider
    public Object[][] sequenceDictionaryNegativeProvider() {
        return new Object[][]{
                {
                        new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 100))),
                        new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 200)))
                },
                {
                        new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 10))),
                        new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr2", 10)))
                },
        };
    }

    @Test(dataProvider = "sequenceDictionaryNegativeProvider", expectedExceptions = PicardException.class)
    public void testSequenceDictionaryNegative(final SAMSequenceDictionary s1, final SAMSequenceDictionary s2 ) {
        SequenceDictionaryUtils.assertSequenceDictionariesEqual(s1, "s1", s2, "s2");
    }

}
