/*
* Adapted from example code in
* Title: Hardcore Java
* Title: Java I/O
* Second Edition: May 2006
* ISBN 10: 0-596-52750-0
* ISBN 13: 9780596527501
*
* http://www.javafaq.nu/java-example-code-1078.html
* 
*/

package org.broad.tribble.util;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;


public final class LittleEndianOutputStream extends FilterOutputStream {

    protected long written;

    public LittleEndianOutputStream(OutputStream out) {
        super(out);
    }

    public void write(int b) throws IOException {
        out.write(b);
        written++;
    }

    public void write(byte[] data, int offset, int length)
            throws IOException {
        out.write(data, offset, length);
        written += length;
    }

    public void writeBoolean(boolean b) throws IOException {
        if (b) this.write(1);
        else this.write(0);
    }

    public void writeByte(int b) throws IOException {
        out.write(b);
        written++;
    }

    public void writeShort(int s) throws IOException {
        out.write(s & 0xFF);
        out.write((s >>> 8) & 0xFF);
        written += 2;
    }

    public void writeChar(int c) throws IOException {
        out.write(c & 0xFF);
        out.write((c >>> 8) & 0xFF);
        written += 2;
    }

    public void writeInt(int i) throws IOException {

        out.write(i & 0xFF);
        out.write((i >>> 8) & 0xFF);
        out.write((i >>> 16) & 0xFF);
        out.write((i >>> 24) & 0xFF);
        written += 4;

    }

    public void writeLong(long l) throws IOException {

        out.write((int) l & 0xFF);
        out.write((int) (l >>> 8) & 0xFF);
        out.write((int) (l >>> 16) & 0xFF);
        out.write((int) (l >>> 24) & 0xFF);
        out.write((int) (l >>> 32) & 0xFF);
        out.write((int) (l >>> 40) & 0xFF);
        out.write((int) (l >>> 48) & 0xFF);
        out.write((int) (l >>> 56) & 0xFF);
        written += 8;

    }

    public final void writeFloat(float f) throws IOException {
        this.writeInt(Float.floatToIntBits(f));
    }

    public final void writeDouble(double d) throws IOException {
        this.writeLong(Double.doubleToLongBits(d));
    }

    public void writeBytes(String s) throws IOException {
        int length = s.length();
        for (int i = 0; i < length; i++) {
            out.write((byte) s.charAt(i));
        }
        written += length;
    }

    /**
     * Srite a string as a null terminated byte array.
     *
     * @param s
     * @throws IOException
     */
    public void writeString(String s) throws IOException {
        writeBytes(s);
        write((byte) 0);
    }

    public long getWrittenCount() {
        return written;
    }

    // Method provide to enable "reseting" to a previous state.
    public void setWrittenCount(long count) {
        this.written = count;
    }
}// end LittleEndianOutputStream
