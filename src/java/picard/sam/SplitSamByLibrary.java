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
package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Command-line program to split a SAM or BAM file into separate files based on
 * library name.
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Takes a SAM or BAM file and separates all the reads " +
                "into one SAM or BAM file per library name.  Reads that do not have " +
                "a read group specified or whose read group does not have a library name " +
                "are written to a file called 'unknown.'  The format (SAM or BAM) of the  " +
                "output files matches that of the input file.  ",
        usageShort = "Splits a SAM or BAM file into individual files by library",
        programGroup = SamOrBam.class
)
public class SplitSamByLibrary extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "The SAM or BAM file to be split. ")
    public File INPUT;
    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "The directory where the library SAM or BAM files should be written " +
                    "(defaults to the current directory). ", optional = true)
    public File OUTPUT = new File(".").getAbsoluteFile();

    private static final Log log = Log.getInstance(SplitSamByLibrary.class);

    public static final int NO_LIBRARIES_SPECIFIED_IN_HEADER = 2;

    public static void main(String[] args) {
        System.exit(new SplitSamByLibrary().instanceMain(args));
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertDirectoryIsWritable(OUTPUT);

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        Map<String, SAMFileWriter> libraryToWriter = new HashMap<String, SAMFileWriter>();
        Map<String, List<SAMReadGroupRecord>> libraryToRg = new HashMap<String, List<SAMReadGroupRecord>>();
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        String extension = reader.type().equals(SamReader.Type.BAM_TYPE) ? ".bam" : ".sam";

        SAMFileHeader unknownHeader = reader.getFileHeader().clone();
        unknownHeader.setReadGroups(new ArrayList<SAMReadGroupRecord>());
        SAMFileWriter unknown = null;

        for (SAMReadGroupRecord rg : reader.getFileHeader().getReadGroups()) {
            String lib = rg.getLibrary();
            if (lib != null) {
                if (!libraryToRg.containsKey(lib)) {
                    libraryToRg.put(lib, new ArrayList<SAMReadGroupRecord>());
                }
                libraryToRg.get(lib).add(rg);
            } else {
                unknownHeader.addReadGroup(rg);
            }
        }

        if (libraryToRg.size() == 0) {
            log.error("No individual libraries are " +
                    "specified in the header of " + INPUT.getAbsolutePath());
            return NO_LIBRARIES_SPECIFIED_IN_HEADER;
        }

        for (Map.Entry<String, List<SAMReadGroupRecord>> entry : libraryToRg.entrySet()) {
            String lib = entry.getKey();
            SAMFileHeader header = reader.getFileHeader().clone();
            header.setReadGroups(entry.getValue());
            libraryToWriter.put(lib, factory.makeSAMOrBAMWriter(header, true,
                    new File(OUTPUT, IOUtil.makeFileNameSafe(lib) + extension)));
        }

        for (Iterator<SAMRecord> it = reader.iterator(); it.hasNext(); ) {
            SAMRecord sam = it.next();
            SAMReadGroupRecord rg = sam.getReadGroup();
            if (rg != null && rg.getLibrary() != null) {
                libraryToWriter.get(rg.getLibrary()).addAlignment(sam);
            } else {
                if (unknown == null) {
                    unknown = factory.makeSAMOrBAMWriter(unknownHeader, true,
                            new File(OUTPUT, "unknown" + extension));
                }
                unknown.addAlignment(sam);
            }
        }

        // Close the reader and writers
        CloserUtil.close(reader);

        if (unknown != null) {
            unknown.close();
        }

        for (SAMFileWriter writer : libraryToWriter.values()) {
            writer.close();
        }

        return 0;
    }
}
