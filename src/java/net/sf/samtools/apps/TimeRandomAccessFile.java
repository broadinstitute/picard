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

/**
 * @author alecw@broadinstitute.org
 */
public class TimeRandomAccessFile {
    public static void main(String[] args) throws Exception {
        RandomAccessFile raf = new RandomAccessFile(new File(args[0]), "r");
        byte[] buf = new byte[64 * 1024];
        long totalBytesRead = 0;
        int bytesRead;
        while ((bytesRead = readBytes(raf, buf, 0, buf.length)) > 0) {
            totalBytesRead += bytesRead;
        }
        System.out.println("Total bytes: " + totalBytesRead);
    }
    private static int readBytes(final RandomAccessFile file, final byte[] buffer, final int offset, final int length)
        throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = file.read(buffer, offset + bytesRead, length - bytesRead);
            if (count <= 0) {
                break;
            }
            bytesRead += count;
        }
        return bytesRead;
    }
}
