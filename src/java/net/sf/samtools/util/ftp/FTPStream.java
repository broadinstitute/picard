/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package net.sf.samtools.util.ftp;

import java.io.FilterInputStream;
import java.io.IOException;

/**
 * A "non-seekable" ftp stream.  This one doesn't support random access.
 *
 * It is assumed that the ftp client has been connected, put in passive mode,
 * set to binary, and otherwise prepped for reading before creating this stream.
 *
 * @author jrobinso
 * @date Oct 31, 2010
 */
public class FTPStream extends FilterInputStream {

    FTPClient ftp;

    public FTPStream(FTPClient ftp) throws IOException {
        super(ftp.getDataStream());
        this.ftp = ftp;
    }


    @Override
    public int read(byte[] bytes, int i, int i1) throws IOException {
        return super.read(bytes, i, i1);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void close() throws IOException {
        super.close();
        ftp.disconnect();
    }
}
