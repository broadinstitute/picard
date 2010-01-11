/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class CreateSequenceDictionaryTest {
    public static File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam");
    public static File BASIC_FASTA = new File(TEST_DATA_DIR, "basic.fasta");
    public static File DUPLICATE_FASTA = new File(TEST_DATA_DIR, "duplicate_sequence_names.fasta");

    @Test
    public void testBasic() throws Exception {
        final File outputFasta = File.createTempFile("CreateSequenceDictionaryTest.", ".fasta");
        outputFasta.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + BASIC_FASTA,
                "OUTPUT=" + outputFasta,
                "TRUNCATE_NAMES_AT_WHITESPACE=false"
        };
        Assert.assertEquals(new CreateSequenceDictionary().instanceMain(argv), 0);
    }

    /**
     * Should throw an exception because with TRUNCATE_NAMES_AT_WHITESPACE, sequence names are not unique.
     */
    @Test(expectedExceptions = PicardException.class)
    public void testNonUniqueSequenceName() throws Exception {
        final File outputFasta = File.createTempFile("CreateSequenceDictionaryTest.", ".fasta");
        outputFasta.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + DUPLICATE_FASTA,
                "OUTPUT=" + outputFasta,
                "TRUNCATE_NAMES_AT_WHITESPACE=true"
        };
        Assert.assertEquals(new CreateSequenceDictionary().instanceMain(argv), 0);
        Assert.fail("Exception should have been thrown.");
    }
}
