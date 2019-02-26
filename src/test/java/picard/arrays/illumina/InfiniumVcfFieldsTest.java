package picard.arrays.illumina;

import htsjdk.samtools.util.Iso8601Date;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class InfiniumVcfFieldsTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_VCF_FILE = new File(TEST_DATA_DIR, "Test.vcf");
    private static final SimpleDateFormat autocallDateFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm");

    @Test
    public void testGetValueFromVcfOtherHeaderLine() throws ParseException {
        try (final VCFFileReader in = new VCFFileReader(TEST_VCF_FILE, false)) {
            final VCFHeader header = in.getFileHeader();
            Assert.assertEquals(InfiniumVcfFields.getIntegerFromVcfOtherHeaderLine(header, InfiniumVcfFields.ANALYSIS_VERSION_NUMBER).intValue(), 1);
            Assert.assertEquals(InfiniumVcfFields.getValueFromVcfOtherHeaderLine(header, InfiniumVcfFields.EXTENDED_ILLUMINA_MANIFEST_VERSION), "1.3");
            Assert.assertEquals(InfiniumVcfFields.getOptionalValueFromVcfOtherHeaderLine(header, InfiniumVcfFields.AUTOCALL_GENDER), "M");
            Assert.assertNull(InfiniumVcfFields.getOptionalValueFromVcfOtherHeaderLine(header, "nothing"));

            final Date truthAutocallDate = new Iso8601Date(autocallDateFormat.parse("04/18/2019 20:57"));
            final Iso8601Date autocallDate = InfiniumVcfFields.getDateFromVcfOtherHeaderLine(header, InfiniumVcfFields.AUTOCALL_DATE, autocallDateFormat);
            Assert.assertEquals(autocallDate, truthAutocallDate);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetNonExistentIntegerFromVcfOtherHeaderLine() {
        try (final VCFFileReader in = new VCFFileReader(TEST_VCF_FILE, false)) {
            final VCFHeader header = in.getFileHeader();
            InfiniumVcfFields.getIntegerFromVcfOtherHeaderLine(header, "nothing");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetNonExistentValueFromVcfOtherHeaderLine() {
        try (final VCFFileReader in = new VCFFileReader(TEST_VCF_FILE, false)) {
            final VCFHeader header = in.getFileHeader();
            InfiniumVcfFields.getValueFromVcfOtherHeaderLine(header, "nothing");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetNonExistentDateFromVcfOtherHeaderLine() {
        try (final VCFFileReader in = new VCFFileReader(TEST_VCF_FILE, false)) {
            final VCFHeader header = in.getFileHeader();
            InfiniumVcfFields.getDateFromVcfOtherHeaderLine(header, "nothing", autocallDateFormat);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetWronglyFormattedDateFromVcfOtherHeaderLine() {
        try (final VCFFileReader in = new VCFFileReader(TEST_VCF_FILE, false)) {
            final VCFHeader header = in.getFileHeader();
            final String badlyFormattedDateString = "04/18/2019";
            header.addMetaDataLine(new VCFHeaderLine("badDate", badlyFormattedDateString));
            InfiniumVcfFields.getDateFromVcfOtherHeaderLine(header, "badDate", autocallDateFormat);
        }
    }
}

