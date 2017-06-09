package picard.sam.util;

import com.beust.jcommander.internal.Nullable;
import org.testng.Assert;
import picard.sam.SamFileConverterTest;
import picard.sam.SamFormatConverter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Pattern;

/**
 * Created by farjoun on 6/8/17.
 *
 * should be moved to htsjdk. do not use
 */
public class SamTestUtils {

    /**
     * A function that checks that without the header two SAM (not bam or cram) files are the same
     * Also checks that header length is the same (in number of lines)
     *
     * @param lhs
     * @param rhs
     */

    public static void assertHeadlessSamFilesAgree(final File lhs, final File rhs) throws IOException {

        Pattern pattern = Pattern.compile("^@...*");

        BufferedReader lhsReader = new BufferedReader(new FileReader(lhs));
        BufferedReader rhsReader = new BufferedReader(new FileReader(rhs));

        //get through the headers
        String lhsLine = lhsReader.readLine();
        int lhsHeaderLength = 0;

        while (lhsLine != null && pattern.matcher(lhsLine).matches()) {
            lhsLine = lhsReader.readLine();
            lhsHeaderLength++;
        }

        String rhsLine = rhsReader.readLine();
        int rhsHeaderLength = 0;
        while (rhsLine != null && pattern.matcher(rhsLine).matches()) {
            rhsLine = rhsReader.readLine();
            rhsHeaderLength++;
        }

        Assert.assertEquals(lhsHeaderLength, rhsHeaderLength, "Header lengths of " + lhs + " and " + rhs + " differ.");

        boolean areSame = true;
        int lineNumber = 1;

        boolean lhsEnd = lhsLine == null;
        boolean rhsEnd = rhsLine == null;

        while (!lhsEnd || !rhsEnd) {
            if (rhsEnd != lhsEnd) {
                throw new AssertionError("File lengths (without headers) are unequal, " + (rhsEnd ? rhs : lhs) + " ended, but " + (rhsEnd  ? lhs : rhs) + "still had lines left");
            }
            // so now both are not null. compare contents:

            Assert.assertEquals(rhsLine, lhsLine, "Files differ at line " + lineNumber + lhsHeaderLength);
            rhsLine = rhsReader.readLine();
            lhsLine = lhsReader.readLine();

            lhsEnd = lhsLine == null;
            rhsEnd = rhsLine == null;

            lineNumber++;
        }
    }
}
