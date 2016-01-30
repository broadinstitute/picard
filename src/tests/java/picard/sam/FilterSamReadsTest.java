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
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileFilter;

public class FilterSamReadsTest extends CommandLineProgramTest {
    @Override
    public String getCommandLineProgramName() {
        return FilterSamReads.class.getSimpleName();
    }

    /**
     * filters a SAM using the available javascript filters
     */
    @Test
    public void testJavaScriptFilters() throws Exception {
        // input as SAM file 
        final File inputSam = new File("testdata/picard/sam/aligned.sam");
        //loop over javascript filters
        for(final File javascriptFile: new File("testdata/picard/sam/FilterSamReads").listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.getName().endsWith(".js");
			}}))
        	{
            final FilterSamReads program = new FilterSamReads();
            program.INPUT = inputSam;
            program.OUTPUT = File.createTempFile("FilterSamReads.output.", ".sam");
            program.OUTPUT.deleteOnExit();
            program.FILTER = FilterSamReads.Filter.includeJavascript;
            program.JAVASCRIPT_FILE = javascriptFile;
            Assert.assertEquals(program.doWork(),0);
        	}        
    	}
	}
