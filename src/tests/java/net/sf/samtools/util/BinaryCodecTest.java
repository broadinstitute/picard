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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;


/*
  * The Broad Institute
  * SOFTWARE COPYRIGHT NOTICE AGREEMENT
  * This software and its documentation are copyright Jan 9, 2009 by the
  * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
  *
  * This software is supplied without any warranty or guaranteed support whatsoever. Neither
  * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
  */

public class BinaryCodecTest {

    @Test
    public void testReadAndWrite() throws IOException {
        final byte[] value = new byte[2];
        value[0] = 1;
        value[1] = 2;
        //Writing to file
        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeBytes(value);
        codec.close();

        //Reading from file
        final byte[] valuesTwo = new byte[2];
        valuesTwo[0] = 1;
        valuesTwo[1] = 2;

        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        final byte[] bytesFromBinaryFile = new byte[2];
        readCodec.readBytes(bytesFromBinaryFile);
        Assert.assertEquals(valuesTwo, bytesFromBinaryFile);
        readCodec.close();
    }

    @Test
    public void testReadAndWriteString() throws IOException {
        final String value = "Test String to Write";

        //Writing to file
        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeString(value, true, false);
        codec.close();

        //Reading from file
        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        final int stringLength = readCodec.readInt();
        Assert.assertEquals(value.length(), stringLength);
        final String s = readCodec.readString(stringLength);
        Assert.assertEquals(value, s);
        readCodec.close();
    }

    @Test
    public void testReadAndWriteInt() throws IOException {
        final int value = 42;

        //Writing to file
        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeInt(value);
        codec.close();

        //Reading from file
        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        Assert.assertEquals(value, readCodec.readInt());
        readCodec.close();
    }

    @Test
    public void testReadAndWriteDouble() throws IOException {
        final double value = 54.4;

        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeDouble(value);
        codec.close();

        //Reading from file
        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        Assert.assertEquals(value, readCodec.readDouble());
        readCodec.close();
    }

    @Test
    public void testReadAndWriteLong() throws IOException {
        final long value = 42;

        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeLong(value);
        codec.close();

        //Reading from file
        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        Assert.assertEquals(value, readCodec.readLong());
        readCodec.close();

    }

    @Test
    public void testReadAndWriteFloat()  throws IOException{
        final float value = 42.5F;

        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeFloat(value);
        codec.close();

        //Reading from file
        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        Assert.assertEquals(value, readCodec.readFloat());
        readCodec.close();
    }

    @Test
    public void testReadAndWriteBoolean()  throws IOException{

        boolean values[] = {true, false};

        for (boolean value : values) {
            final File outputFile = File.createTempFile("temp", ".bin");
            final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
            final BinaryCodec codec = new BinaryCodec(stream);
            codec.writeBoolean(value);
            codec.close();

            //Reading from file
            final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
            final BinaryCodec readCodec = new BinaryCodec(instream);
            Assert.assertEquals(value, readCodec.readBoolean());
            readCodec.close();
        }
    }

    @Test
    public void testReadAndWriteMutlitpleData()  throws IOException{
        final float fValue = 42.5F;
        final String sValue = "TestString";

        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        codec.writeFloat(fValue);
        codec.writeString(sValue, true, false);
        codec.close();

        //Reading from file
        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);
        Assert.assertEquals(fValue, readCodec.readFloat());
        final int stringLength = readCodec.readInt();
        Assert.assertEquals(sValue, readCodec.readString(stringLength));
        readCodec.close();
    }

    @Test
    public void readPastEndOfFile() throws IOException{
        final long startTime = System.currentTimeMillis();
        int i = 0;

        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        while (i<100){
            codec.writeInt(i);
            i++;
        }
        codec.close();

        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);

        System.out.println((System.currentTimeMillis() - startTime) + "ms to write");
        int z = 0;
        boolean reachedStatement = false;
        while (z<1000) {
            try {
                Assert.assertEquals(z, readCodec.readInt());
            } catch (Exception e) {
               Assert.assertEquals(RuntimeEOFException.class, e.getClass());
                reachedStatement = true;
            }
            z++;
        }

        Assert.assertTrue(reachedStatement);
        readCodec.close();
    }

    @Test
    public void timeTest() throws IOException{
        final long startTime = System.currentTimeMillis();
        int i = 0;

        final File outputFile = File.createTempFile("temp", ".bin");
        final DataOutputStream stream = new DataOutputStream(new FileOutputStream(outputFile));
        final BinaryCodec codec = new BinaryCodec(stream);
        while (i<100){
            codec.writeInt(i);
            i++;
        }
        codec.close();

        final DataInputStream instream = new DataInputStream(new FileInputStream(outputFile));
        final BinaryCodec readCodec = new BinaryCodec(instream);

        System.out.println((System.currentTimeMillis() - startTime) + "ms to write");
        int z = 0;
        boolean reachedStatement = false;
        while (z<1000) {
            try {
                Assert.assertEquals(z, readCodec.readInt());
            } catch (Exception e) {
               Assert.assertEquals(RuntimeEOFException.class, e.getClass());
                reachedStatement = true;
            }
            z++;
        }

        Assert.assertTrue(reachedStatement);
        readCodec.close();
    }


}
