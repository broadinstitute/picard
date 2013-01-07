/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
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

import net.sf.samtools.SAMException;

import java.io.*;
import java.net.Socket;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

/**
 * @author jrobinso
 * @date Oct 30, 2010
 */
public class FTPClient {

    private Socket commandSocket = null;

    public static int READ_TIMEOUT = 5 * 60 * 1000;

    /**
     * Stream to write commands.
     * NOTE -- a PrintStream is used no purpose (as opposed to PrintWriter).  PrintWriter will not work!
     */
    private PrintStream commandStream = null;
    private BufferedReader responseReader = null;
    private InputStream dataStream;
    private String passiveHost;
    private int passivePort;
    long restPosition = -1;
    String host;

    /**
     * Connects to the given FTP host on the default port.
     */
    public FTPReply connect(String host) throws IOException {
        this.host = host;
        commandSocket = new Socket(host, 21);
        commandSocket.setSoTimeout(READ_TIMEOUT);
        commandStream = new PrintStream(commandSocket.getOutputStream());
        responseReader = new BufferedReader(new InputStreamReader(commandSocket.getInputStream()));

        FTPReply reply = new FTPReply(responseReader);

        if (!reply.isPositiveCompletion()) {
            disconnect();
        }

        return reply;
    }


    /**
     * Executes the given FTP command on our current connection, returning the
     * three digit response code from the server.  This method only works for
     * commands that do not require an additional data port.
     */
    public FTPReply executeCommand(String command) throws IOException {
        commandStream.println(command);
        return new FTPReply(responseReader);
    }


    /**
     * Wrapper for the commands <code>user [username]</code> and <code>pass
     * [password]</code>.
     */
    public FTPReply login(String username, String password) throws IOException {
        FTPReply response = executeCommand("user " + username);
        if (!response.isPositiveIntermediate()) return response;
        response = executeCommand("pass " + password);
        return response;
    }

    public FTPReply quit() throws IOException {
        return executeCommand("QUIT");
    }

    public FTPReply binary() throws IOException {
        return executeCommand("TYPE I");
    }


    public FTPReply pasv() throws IOException {

        FTPReply reply = executeCommand("PASV");

        if (reply.getCode() == 226 || reply.getCode() == 426) {
            reply = getReply();
        }

        String response = reply.getReplyString();


        int code = reply.getCode();

        int opening = response.indexOf('(');
        int closing = response.indexOf(')', opening + 1);
        if (closing > 0) {
            String dataLink = response.substring(opening + 1, closing);
            StringTokenizer tokenizer = new StringTokenizer(dataLink, ",");
            try {
                passiveHost = tokenizer.nextToken() + "." + tokenizer.nextToken() + "."
                        + tokenizer.nextToken() + "." + tokenizer.nextToken();
                passivePort = Integer.parseInt(tokenizer.nextToken()) * 256
                        + Integer.parseInt(tokenizer.nextToken());
            } catch (NumberFormatException e) {
                throw new IOException("SimpleFTP received bad data link information: " + response);
            } catch (NoSuchElementException e){
                throw new IOException("SimpleFTP received bad data link information: " + response);
            }
        }

        if (reply.isPositiveCompletion()) {
            if (dataStream == null) {
                Socket dataSocket = new Socket(passiveHost, passivePort);
                dataSocket.setSoTimeout(READ_TIMEOUT);
                dataStream = new SocketInputStream(dataSocket, dataSocket.getInputStream());
            }
        }
        return reply;
    }

    public void setRestPosition(long position) {
        this.restPosition = position;
    }

    public FTPReply retr(String file) throws IOException {

        if (restPosition >= 0) {
            FTPReply restReply = executeCommand("REST " + restPosition);
            if (!restReply.isSuccess()) {
                return restReply;
            }
        }

        return executeCommand("RETR " + file);
    }

    public FTPReply getReply() throws IOException {
        return new FTPReply(responseReader);
    }

    /**
     * Return the size of the remote file
     *
     * @param file
     * @return
     * @throws IOException
     */
    public FTPReply size(String file) throws IOException {

        return executeCommand("SIZE " + file);

    }


    public InputStream getDataStream() throws IOException {
        return dataStream;
    }

    public void closeDataStream() throws IOException {
        //  NOTE -- some ftp servers seem to need a pause before closing the data stream
        // if (dataStream != null) {
        //    try {
        //        //
        //        Thread.sleep(3000);
        //    } catch (InterruptedException e) {
        //
        //    }
        dataStream.close();
        //}
        dataStream = null;
    }


    /**
     * Disconnects from the host to which we are currently connected.
     */
    public void disconnect() {
        try {
            //quit();
            if (commandStream != null) {
                commandStream.close();
                responseReader.close();
                commandSocket.close();

                if (dataStream != null) {
                    dataStream.close();
                }
            }
        } catch (IOException e) {
            throw new SAMException("Error disconnecting", e);
        }

        commandStream = null;
        responseReader = null;
        commandSocket = null;
    }

    class SocketInputStream extends FilterInputStream {

        Socket socket;

        SocketInputStream(Socket socket, InputStream inputStream) {
            super(inputStream);
            this.socket = socket;
        }

        @Override
        public void close() throws IOException {
            super.close();
            socket.close();
            FTPClient.this.dataStream = null;
        }
    }

}
