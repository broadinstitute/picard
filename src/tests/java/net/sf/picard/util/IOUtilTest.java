/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import net.sf.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class IOUtilTest {
    private File existingTempFile;
    private String systemTempDir;

    @BeforeClass
    public void setUp() throws IOException {
        existingTempFile = File.createTempFile("FiletypeTest.", ".tmp");
        existingTempFile.deleteOnExit();
        systemTempDir = System.getProperty("java.io.tmpdir");
        final File tmpDir = new File(systemTempDir);
        if (!tmpDir.isDirectory()) tmpDir.mkdir();
        if (!tmpDir.isDirectory())
            throw new RuntimeException("java.io.tmpdir (" + systemTempDir + ") is not a directory");
    }

    @Test(dataProvider = "fileTypeTestCases")
    public void testFileType(final String path, boolean expectedIsRegularFile) {
        final File file = new File(path);
        Assert.assertEquals(IOUtil.isRegularPath(file), expectedIsRegularFile);
    }

    @DataProvider(name = "fileTypeTestCases")
    private Object[][] fileTypeTestCases() {
        return new Object[][] {
                {"/dev/null",   Boolean.FALSE},
                {"/dev/stdout", Boolean.FALSE},
                {"/non/existent/file", Boolean.TRUE},
                {existingTempFile.getAbsolutePath(), Boolean.TRUE},
                {systemTempDir, Boolean.FALSE}

        };
    }
}
