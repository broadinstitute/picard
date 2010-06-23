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
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

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
        for(int i = 0; (line = reader.readLine()) != null; ++i) {
            Assert.assertEquals(line + "\n", linesWritten.get(i));
        }
    }

    // Writing to /dev/null doesn't work properly on some systems, because FileDescriptor.sync() fails.
    @Test(groups = "broken")
    public void testOverflow() throws Exception {
        final File f = File.createTempFile("BCOST.", ".gz");
        f.deleteOnExit();
        final List<String> linesWritten = new ArrayList<String>();
        System.out.println("Creating file " + f);
        final BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(f);
        Random r = new Random(15555);
        final int INPUT_SIZE = 64 * 1024;
        byte[] input = new byte[INPUT_SIZE];
        r.nextBytes(input);
        bcos.write(input);
        bcos.close();

        final BlockCompressedInputStream bcis = new BlockCompressedInputStream(f);
        byte[] output = new byte[INPUT_SIZE];
        int len;
        int i = 0;
        while ((len = bcis.read(output, 0, output.length)) != -1) {
            for (int j = 0; j < len; j++) {
               Assert.assertEquals(output[j], input[i++]);
            }
        }
        Assert.assertEquals(i, INPUT_SIZE);
    }

    // PIC-393 exception closing BGZF stream opened to /dev/null
    // I don't think this will work on Windows, because /dev/null doesn't work
    @Test(groups = "broken")
    public void testDevNull() throws Exception {
        final BlockCompressedOutputStream bcos = new BlockCompressedOutputStream("/dev/null");
        bcos.write("Hi, Mom!".getBytes());
        bcos.close();
    }
}
