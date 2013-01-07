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

import java.io.BufferedReader;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Oct 30, 2010
 */
public class FTPReply {

    String reply;
    int code;

    public FTPReply(BufferedReader inputStream) throws IOException {

        String response = null;
        do {
            response = inputStream.readLine();
        } while (response != null &&
                !(Character.isDigit(response.charAt(0)) &&
                        Character.isDigit(response.charAt(1)) &&
                        Character.isDigit(response.charAt(2)) &&
                        response.charAt(3) == ' '));
        if (response == null || response.length() < 3) {
            code = -1;
        } else {
            code = Integer.parseInt(response.substring(0, 3));
            reply = response.substring(3).trim();

        }
    }

    /**
     * Gets server reply code from the control port after an ftp command has
     * been executed.  It knows the last line of the response because it begins
     * with a 3 digit number and a space, (a dash instead of a space would be a
     * continuation).
     */

    public int getCode() throws IOException {
        return code;
    }


    /**
     * Gets server reply string from the control port after an ftp command has
     * been executed.  This consists only of the last line of the response,
     * and only the part after the response code.
     */
    public String getReplyString()
            throws IOException {

        return reply;
    }


    public boolean isSuccess() {
        return isPositiveCompletion() || isPositiveIntermediate();
    }

    /**
     * Determine if a reply code is a positive completion response.  All
     * codes beginning with a 2 are positive completion responses.
     * The FTP server will send a positive completion response on the final
     * successful completion of a command.
     * <p/>
     *
     * @return True if a reply code is a postive completion response, false
     *         if not.
     *         *
     */
    public boolean isPositiveCompletion() {
        return (code >= 200 && code < 300);
    }


    /**
     * Determine if a reply code is a positive intermediate response.  All
     * codes beginning with a 3 are positive intermediate responses.
     * The FTP server will send a positive intermediate response on the
     * successful completion of one part of a multi-part sequence of
     * commands.  For example, after a successful USER command, a positive
     * intermediate response will be sent to indicate that the server is
     * ready for the PASS command.
     * <p/>
     *
     * @return True if a reply code is a postive intermediate response, false
     *         if not.
     *         *
     */
    public boolean isPositiveIntermediate() {
        return (code >= 300 && code < 400);
    }
}
