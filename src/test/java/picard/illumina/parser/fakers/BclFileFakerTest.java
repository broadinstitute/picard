package picard.illumina.parser.fakers;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.illumina.parser.readers.BclReader;

import java.io.File;

/**
 * The MIT License
 * <p/>
 * Copyright (c) 2014 The Broad Institute
 * <p/>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p/>
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * <p/>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
public class BclFileFakerTest {

    @Test
    public void testFileLengthMatchesHeaderLength() throws Exception {
        final File fakeFile = File.createTempFile("BclFileFakerTest", ".bcl");
        fakeFile.deleteOnExit();

        new BclFileFaker().fakeFile(fakeFile, 100000);
        // .make() has a number of checks for the file
        final BclReader bclReader = new BclReader(
            fakeFile,
            new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY), false);
        Assert.assertEquals(100000, BclReader.getNumberOfClusters(fakeFile));
        Assert.assertEquals(BclReader.getNumberOfClusters(fakeFile), fakeFile.length() - 4);
    }

    @Test
    public void testGZFileIsActuallyGZipped() throws Exception {
        final File fakeFile = File.createTempFile("BclFileFakerTest", ".bcl.gz");
        fakeFile.deleteOnExit();

        new BclFileFaker().fakeFile(fakeFile, 100000);
        new BclReader(
                fakeFile,
                new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY), false);
    }
}
