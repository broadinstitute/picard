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
package net.sf.samtools;

import net.sf.samtools.SAMFileReader.ValidationStringency;
import java.io.*;

public class FixBAMFile
{
    public static void main(String[] args) {
        File inputFile = new File(args[0]);
        File outputFile = new File(args[1]);
        SAMFileReader reader = new SAMFileReader(inputFile);
        reader.setValidationStringency(ValidationStringency.SILENT);
        SAMFileHeader header = reader.getFileHeader();
        SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, outputFile);
        for (SAMRecord record : reader) {
            if (record.getIndexingBin() != null) {
                record.setIndexingBin(record.computeIndexingBin());
            }
            writer.addAlignment(record);
        }
        writer.close();
    }
}
