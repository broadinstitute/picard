/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
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
