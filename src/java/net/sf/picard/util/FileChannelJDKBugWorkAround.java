/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

import net.sf.picard.PicardException;

import java.io.File;
import java.io.FileWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.ByteBuffer;

/**
 * Work-around for the following bug
 * http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6503430
 *
 * Call the method in the main thread before launching threads that do memory-mapping.
 * @author alecw@broadinstitute.org
 */
public class FileChannelJDKBugWorkAround {
    public static byte doBugWorkAround() {
        try {
            File tmpFile = File.createTempFile("ignore-me.", ".bug-work-around");
            FileWriter writer = new FileWriter(tmpFile);
            writer.write("Hi, Mom!");
            writer.close();
            FileInputStream is = new FileInputStream(tmpFile);
            ByteBuffer buf = is.getChannel().map(FileChannel.MapMode.READ_ONLY, 0, tmpFile.length());
            is.close();
            byte ret = buf.get();
            tmpFile.delete();
            return ret;
        } catch (IOException e) {
            throw new PicardException("IOException", e);
        }
    }
}
