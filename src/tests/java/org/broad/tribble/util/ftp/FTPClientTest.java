package org.broad.tribble.util.ftp;

import net.sf.samtools.util.ftp.FTPClient;
import net.sf.samtools.util.ftp.FTPReply;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.InputStream;
import java.net.UnknownHostException;




/**
* @author Jim Robinson
* @since 10/3/11
*/
public class FTPClientTest {

    static String host = "ftp.broadinstitute.org";
    static String file = "/pub/igv/TEST/test.txt";
    static int fileSize = 27;
    static byte[] expectedBytes = "abcdefghijklmnopqrstuvwxyz\n".getBytes();
    FTPClient client;

    @BeforeMethod
    public void setUp() throws IOException {
        client = new FTPClient();
        FTPReply reply = client.connect(host);
        Assert.assertTrue(reply.isSuccess());
    }

    @AfterMethod
    public void tearDown() {
        System.out.println("Disconnecting");
        client.disconnect();
    }

    @Test
    public void testLogin() throws Exception {

    }

    @Test
    public void testPasv() throws Exception {
        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            Assert.assertTrue(reply.isSuccess());

            reply = client.pasv();
            Assert.assertTrue(reply.isSuccess());
        } finally {
            client.closeDataStream();
        }
    }

    @Test
    public void testSize() throws Exception {

        FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
        Assert.assertTrue(reply.isSuccess());

        reply = client.binary();
        Assert.assertTrue(reply.isSuccess());

        reply = client.size(file);
        String val = reply.getReplyString();
        int size = Integer.parseInt(val);
        Assert.assertEquals(fileSize, size);
    }


    @Test
    public void testDownload() throws Exception {
        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            Assert.assertTrue(reply.isSuccess());

            reply = client.binary();
            Assert.assertTrue(reply.isSuccess());

            reply = client.pasv();
            Assert.assertTrue(reply.isSuccess());

            reply = client.retr(file);
            Assert.assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();
            int idx = 0;
            int b;
            while ((b = is.read()) >= 0) {
                Assert.assertEquals(expectedBytes[idx], (byte) b);
                idx++;
            }

        } finally {
            client.closeDataStream();
            FTPReply reply = client.retr(file);
            System.out.println(reply.getCode());
            Assert.assertTrue(reply.isSuccess());

        }
    }


    @Test
    public void testRest() throws Exception {
        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            Assert.assertTrue(reply.isSuccess());

            reply = client.binary();
            Assert.assertTrue(reply.isSuccess());

            reply = client.pasv();
            Assert.assertTrue(reply.isSuccess());

            final int restPosition = 5;
            client.setRestPosition(restPosition);

            reply = client.retr(file);
            Assert.assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();
            int idx = restPosition;
            int b;
            while ((b = is.read()) >= 0) {
                Assert.assertEquals(expectedBytes[idx], (byte) b);
                idx++;
            }

        } finally {
            client.closeDataStream();
            FTPReply reply = client.retr(file);
            System.out.println(reply.getCode());
            Assert.assertTrue(reply.isSuccess());

        }
    }

    /**
     * Test accessing a non-existent file
     */
    @Test
    public void testNonExistentFile() throws Exception {

        String host = "ftp.broadinstitute.org";
        String file = "/pub/igv/TEST/fileDoesntExist.txt";
        FTPClient client = new FTPClient();

        FTPReply reply = client.connect(host);
        Assert.assertTrue(reply.isSuccess());

        reply = client.login("anonymous", "igv@broadinstitute.org");
        Assert.assertTrue(reply.isSuccess());

        reply = client.binary();
        Assert.assertTrue(reply.isSuccess());

        reply = client.executeCommand("size " + file);
        Assert.assertEquals(550, reply.getCode());

        client.disconnect();

    }

        /**
     * Test accessing a non-existent server
     */
    @Test
    public void testNonExistentServer() throws Exception {

        String host = "ftp.noSuchServer.org";
        String file = "/pub/igv/TEST/fileDoesntExist.txt";
        FTPClient client = new FTPClient();

        FTPReply reply = null;
        try {
            reply = client.connect(host);
        } catch (UnknownHostException e) {
            // This is expected
        }


        client.disconnect();

    }


    @Test
    public void testMultiplePasv() throws Exception {

        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            Assert.assertTrue(reply.isSuccess());

            reply = client.pasv();
            Assert.assertTrue(reply.isSuccess());
            client.closeDataStream();

            reply = client.pasv();
            Assert.assertTrue(reply.isSuccess());
            client.closeDataStream();


        }

        finally {

        }
    }


    @Test
    public void testMultipleRest() throws Exception {
        FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
        Assert.assertTrue(reply.isSuccess());

        reply = client.binary();
        Assert.assertTrue(reply.isSuccess());

        restRetr(5, 10);
        restRetr(2, 10);
        restRetr(15, 10);

    }

    private void restRetr(int restPosition, int length) throws IOException {

        try {

            if (client.getDataStream() == null) {
                FTPReply reply = client.pasv();
                Assert.assertTrue(reply.isSuccess());
            }

            client.setRestPosition(restPosition);

            FTPReply reply = client.retr(file);
            //assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();

            byte[] buffer = new byte[length];
            is.read(buffer);


            for (int i = 0; i < length; i++) {
                System.out.print((char) buffer[i]);
                Assert.assertEquals(expectedBytes[i + restPosition], buffer[i]);
            }
            System.out.println();
        }

        finally {
            client.closeDataStream();
            FTPReply reply = client.getReply();  // <== MUST READ THE REPLY
            System.out.println(reply.getReplyString());


        }
    }
}
