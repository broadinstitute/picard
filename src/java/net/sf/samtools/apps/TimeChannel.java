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
 */
package net.sf.samtools.apps;

import java.io.RandomAccessFile;
import java.io.File;
import java.io.IOException;
import java.io.FileInputStream;
import java.nio.channels.FileChannel;
import java.nio.MappedByteBuffer;

/**
 * @author alecw@broadinstitute.org
 */
public class TimeChannel {
    public static void main(String[] args) throws Exception {
        long fileSize = new File(args[0]).length();
        FileInputStream in = new FileInputStream(args[0]);
        FileChannel channel = in.getChannel();
        byte[] buf = new byte[64 * 1024];
        long totalBytesRead = 0;
        long mappedOffset = 0;
        // Map a round number of bytes that might correspond to some kind of page size boundary.
        long maxToMapAtATime = 1024 * 1024 * 1024;
        long mappedSize = Math.min(fileSize, maxToMapAtATime);
        while (totalBytesRead < fileSize) {
            System.err.println("mappedOffset: " + mappedOffset + "; mappedSize: " + mappedSize);
            System.err.println("fileSize: " + fileSize + "; totalBytesRead: " + totalBytesRead);
            MappedByteBuffer mappedBuffer = channel.map(FileChannel.MapMode.READ_ONLY, mappedOffset, mappedSize);
            while (mappedBuffer.remaining() > 0) {
                int bytesRead = Math.min(mappedBuffer.remaining(), buf.length);
                mappedBuffer.get(buf, 0, bytesRead);
                totalBytesRead += bytesRead;
            }
            mappedOffset += mappedSize;
            mappedSize = Math.min(fileSize - totalBytesRead, maxToMapAtATime);
        }
        System.out.println("Total bytes: " + totalBytesRead);
    }
}
