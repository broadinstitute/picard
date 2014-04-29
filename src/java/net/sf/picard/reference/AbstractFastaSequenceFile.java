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

package net.sf.picard.reference;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.util.BufferedLineReader;

import java.io.File;
import java.io.FileInputStream;

/**
 * Provide core sequence dictionary functionality required by all fasta file readers.
 * @author Matt Hanna
 */
abstract class AbstractFastaSequenceFile implements ReferenceSequenceFile {
    protected final File file;
    protected SAMSequenceDictionary sequenceDictionary;

    /**
     * Finds and loads the sequence file dictionary.
     * @param file Fasta file to read.  Also acts as a prefix for supporting files.
     */
    AbstractFastaSequenceFile(final File file) {
        this.file = file;
        final File dictionary = findSequenceDictionary(file);

        if (dictionary != null) {
            IoUtil.assertFileIsReadable(dictionary);

            try {
                final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                final SAMFileHeader header = codec.decode(new BufferedLineReader(new FileInputStream(dictionary)),
                        dictionary.toString());
                if (header.getSequenceDictionary() != null && header.getSequenceDictionary().size() > 0) {
                    this.sequenceDictionary = header.getSequenceDictionary();
                }
            }
            catch (Exception e) {
                throw new PicardException("Could not open sequence dictionary file: " + dictionary, e);
            }
        }
    }

    protected static File findSequenceDictionary(final File file) {
        // Try and locate the dictionary
        String dictionaryName = file.getAbsolutePath();
        String dictionaryNameExt = file.getAbsolutePath();
        boolean fileTypeSupported = false;
        for (final String extension : ReferenceSequenceFileFactory.FASTA_EXTENSIONS) {
            if (dictionaryName.endsWith(extension)) {
                  dictionaryNameExt = new String(dictionaryName);
                  dictionaryNameExt += ".dict";
                  dictionaryName = dictionaryName.substring(0, dictionaryName.lastIndexOf(extension));
                  dictionaryName += ".dict";
                  fileTypeSupported = true;
                  break;
            }
        }
        if (!fileTypeSupported)
            throw new IllegalArgumentException("File is not a supported reference file type: " + file.getAbsolutePath());

        final File dictionary = new File(dictionaryName);
        if (dictionary.exists())
            return dictionary;
        // try without removing the file extension
        final File dictionaryExt = new File(dictionaryNameExt);
        if (dictionaryExt.exists())
            return dictionaryExt;
        else return null;
    }

    /**
     * Returns the list of sequence records associated with the reference sequence if found
     * otherwise null.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return this.sequenceDictionary;
    }

    /** Returns the full path to the reference file. */
    public String toString() {
        return this.file.getAbsolutePath();
    }

    /** default implementation -- override if index is supported */
    public boolean isIndexed() {return false;}

    /** default implementation -- override if index is supported */
    public ReferenceSequence getSequence( String contig ) {
        throw new UnsupportedOperationException();
    }

    /** default implementation -- override if index is supported */
    public ReferenceSequence getSubsequenceAt( String contig, long start, long stop ) {
        throw new UnsupportedOperationException();
    }

}
