/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
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
