package picard.analysis;

import htsjdk.samtools.util.SequenceUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class AdapterUtilityTest {
    private static final AdapterUtility adapterUtility = new AdapterUtility(AdapterUtility.DEFAULT_ADAPTER_SEQUENCE);

    @DataProvider
    Object[][] testAdapterReadsData(){
        return new Object[][]{
                {"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAAGAGGTTGTGGGTATGGGGGGGGGGGGGGGGGGGGTGGGGAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"},
                {"GATCGGAAGAGCACACGTCTGAACTCCAGTAACGCACATCTATCTCGTATTCCGTCTTCTTCTTGAACAGGGGTGGGTGGGGGGGGGGGGTGGGGGGGGGGGGGG"},
                {"CTGAACTCCAGTCACGCACATCTATCGCGTATGCCGTCTTCTGCTTGAAAAGGGGGGGGGGGTGGGGGGGGGGG"},
                {"GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAAGAGGTTGTGGGGGTGAAGGGGGGGGGGGGGGG"},
                {"AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGGGGGGGGGGGGTGGGTCCA"},
                {"AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGGGGGGGGGGGGTGGGTCCA"},
                {"CTCGTATGCCGTCTTCTGCTTGGGGGGGGGGGGTGGGTCCA"},
                {"GAGCTCGTATGCCGTCTTCTGCTTGGGGGGGGGGGGTGGGTCCA"},
                {AdapterUtility.DEFAULT_ADAPTER_SEQUENCE.get(0)},
                {AdapterUtility.DEFAULT_ADAPTER_SEQUENCE.get(1)},
                {AdapterUtility.DEFAULT_ADAPTER_SEQUENCE.get(2)},
                {AdapterUtility.DEFAULT_ADAPTER_SEQUENCE.get(3)},
                {AdapterUtility.DEFAULT_ADAPTER_SEQUENCE.get(4)},
                {AdapterUtility.DEFAULT_ADAPTER_SEQUENCE.get(5)}
        };
    }

    @Test(dataProvider = "testAdapterReadsData")
    public void testAdapterReads(final String readBases){
        Assert.assertTrue(adapterUtility.isAdapterSequence(readBases.getBytes()));
    }

    @Test(dataProvider = "testAdapterReadsData")
    public void testAdapterReadsRevComp(final String readBases){
        Assert.assertTrue(adapterUtility.isAdapterSequence(SequenceUtil.reverseComplement(readBases).getBytes(),true));
    }
}