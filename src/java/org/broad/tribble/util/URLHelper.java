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

package org.broad.tribble.util;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * Interface defining a helper class for dealing with URL resources.  Purpose of this class is to provide the
 * flexibility to use either URLConnection or Apache HTTPClient. Also want to delegate to either HTTP or FTP
 *
 * @see HTTPHelper
 * @see FTPHelper
 * @author jrobinso
 * @date Jun 28, 2011
 */
public interface URLHelper {

    URL getUrl();

    long getContentLength() throws IOException;

    InputStream openInputStream() throws IOException;

    /**
     * May throw an OperationUnsupportedException
     * @deprecated Will be removed in a future release, as is somewhat fragile
     * and not used.
     * @param start
     * @param end
     * @return
     * @throws IOException
     */
    @Deprecated
    InputStream openInputStreamForRange(long start, long end) throws IOException;

    public boolean exists() throws IOException;

}
