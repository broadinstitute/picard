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
package net.sf.samtools.util;

import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.ArrayList;

public class BlockCompressedOutputStreamTest {

    @Test
    public void testBasic() throws Exception {
        final File f = File.createTempFile("BCOST.", ".gz");
        f.deleteOnExit();
        final List<String> linesWritten = new ArrayList<String>();
        System.out.println("Creating file " + f);
        final BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(f);
        String s = "Hi, Mom!\n";
        bcos.write(s.getBytes());
        linesWritten.add(s);
        s = "Hi, Dad!\n";
        bcos.write(s.getBytes());
        linesWritten.add(s);
        bcos.flush();
        final StringBuilder sb = new StringBuilder(BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE * 2);
        s = "1234567890123456789012345678901234567890123456789012345678901234567890\n";
        while (sb.length() <= BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE) {
            sb.append(s);
            linesWritten.add(s);
        }
        bcos.write(sb.toString().getBytes());
        bcos.close();
        final BlockCompressedInputStream bcis = new BlockCompressedInputStream(f);
        final BufferedReader reader = new BufferedReader(new InputStreamReader(bcis));
        String line;
        // I don't understand why I need to check ready(), but waiting for null return
        // causes an exception.
        for(int i = 0; reader.ready() && (line = reader.readLine()) != null; ++i) {
            Assert.assertEquals(line + "\n", linesWritten.get(i));
        }
    }
}
