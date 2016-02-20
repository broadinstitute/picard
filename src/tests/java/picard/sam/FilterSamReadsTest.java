/*
 * The MIT License
 *
 * Pierre Lindenbaum - Institut du Thorax - 2016
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
package picard.sam;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;

public class FilterSamReadsTest extends CommandLineProgramTest {
    @Override
    public String getCommandLineProgramName() {
        return FilterSamReads.class.getSimpleName();
    }

    @DataProvider(name = "dataTestJsFilter")
    public Object[][] dataTestJsFilter() {
        return new Object[][]{
                {"testdata/picard/sam/aligned.sam","testdata/picard/sam/FilterSamReads/filterOddStarts.js",3},
                {"testdata/picard/sam/aligned.sam","testdata/picard/sam/FilterSamReads/filterReadsWithout5primeSoftClip.js", 0}
        };
    }
    
    /**
     * filters a SAM using a javascript filter
     */
    @Test(dataProvider = "dataTestJsFilter")
    public void testJavaScriptFilters(final String samFilename, final String javascriptFilename,final int expectNumber) throws Exception {
        // input as SAM file 
        final File inputSam = new File(samFilename);
        final File javascriptFile = new File(javascriptFilename);
        //loop over javascript filters
        final FilterSamReads program = new FilterSamReads();
        program.INPUT = inputSam;
        program.OUTPUT = File.createTempFile("FilterSamReads.output.", ".sam");
        program.OUTPUT.deleteOnExit();
        program.FILTER = FilterSamReads.Filter.includeJavascript;
        program.JAVASCRIPT_FILE = javascriptFile;
        Assert.assertEquals(program.doWork(),0);
        
        //count reads
        final SamReader samReader = SamReaderFactory.makeDefault().open(program.OUTPUT);
        final SAMRecordIterator iter = samReader.iterator();
        int count = 0;
        while(iter.hasNext()) { iter.next(); ++ count; }
        iter.close();
        samReader.close();
        Assert.assertEquals(count, expectNumber);
    	}
	}
