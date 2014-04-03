/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.vcf;

import org.broadinstitute.variant.VariantBaseTest;
import org.testng.annotations.Test;
import org.testng.Assert;


/**
 * User: ebanks
 * Date: Apr 2, 2014
 */
public class VCFCompoundHeaderLineUnitTest extends VariantBaseTest {

    @Test
    public void supportsVersionFields() {
	final String line = "<ID=FOO,Number=1,Type=Float,Description=\"foo\",Version=3>";
	final VCFCompoundHeaderLine headerline = new VCFInfoHeaderLine(line, VCFHeaderVersion.VCF4_2);
	// if we don't support version fields then we should fail before we ever get here
	Assert.assertTrue(true);
    }
}
