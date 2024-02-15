/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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
package picard.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;

import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.math.BigInteger;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Class with helper methods for generating and writing SequenceDictionary objects.
 */
public class SequenceDictionaryUtils {

    static public class SamSequenceRecordsIterator implements Iterator<SAMSequenceRecord> {
        final private boolean truncateNamesAtWhitespace;
        final private ReferenceSequenceFile refSeqFile;
        private String genomeAssembly;
        private String uri;

        public void setGenomeAssembly(final String genomeAssembly) {
            this.genomeAssembly = genomeAssembly;
        }

        public void setUri(final String uri) {
            this.uri = uri;
        }

        public void setSpecies(final String species) {
            this.species = species;
        }

        private String species;
        private ReferenceSequence nextRefSeq;
        private final MessageDigest md5;

        public SamSequenceRecordsIterator(final Path referenceSequence, final boolean truncateNamesAtWhitespace) {
            this.truncateNamesAtWhitespace = truncateNamesAtWhitespace;
            this.refSeqFile = ReferenceSequenceFileFactory.
                    getReferenceSequenceFile(referenceSequence, truncateNamesAtWhitespace);

            this.nextRefSeq = refSeqFile.nextSequence();
            try {
                md5 = MessageDigest.getInstance("MD5");
            } catch (NoSuchAlgorithmException e) {
                throw new PicardException("MD5 algorithm not found", e);
            }
        }

        private String md5Hash(final byte[] bytes) {
            md5.reset();
            md5.update(bytes);
            String s = new BigInteger(1, md5.digest()).toString(16);
            if (s.length() != 32) {
                final String zeros = "00000000000000000000000000000000";
                s = zeros.substring(0, 32 - s.length()) + s;
            }
            return s;
        }

        /**
         * Create one SAMSequenceRecord from a single fasta sequence
         */
        private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) {
            final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName(), refSeq.length());

            // Compute MD5 of upcased bases
            final byte[] bases = refSeq.getBases();
            for (int i = 0; i < bases.length; ++i) {
                bases[i] = StringUtil.toUpperCase(bases[i]);
            }

            ret.setAttribute(SAMSequenceRecord.MD5_TAG, md5Hash(bases));
            if (genomeAssembly != null) {
                ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, genomeAssembly);
            }
            ret.setAttribute(SAMSequenceRecord.URI_TAG, uri);
            if (species != null) {
                ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, species);
            }
            return ret;
        }

        @Override
        public boolean hasNext() {
            return nextRefSeq != null;
        }

        @Override
        public SAMSequenceRecord next() {
            if (!hasNext()) {
                throw new NoSuchElementException("next() was called when hasNext() was false.");
            }
            final SAMSequenceRecord samSequenceRecord = makeSequenceRecord(nextRefSeq);
            nextRefSeq = refSeqFile.nextSequence();
            return samSequenceRecord;
        }
    }

    /**
     * Encodes a sequence dictionary
     *
     * @param writer                    a Buffered writer into which the dictionary will be written
     * @param samSequenceRecordIterator an iterator that produces SAMSequenceRecords
     * @throws IllegalArgumentException if the iterator produces two SAMSequenceRecord with the same name
     */
    public static void encodeDictionary(final BufferedWriter writer, Iterator<SAMSequenceRecord> samSequenceRecordIterator) {
        final SAMSequenceDictionaryCodec samDictCodec = new SAMSequenceDictionaryCodec(writer);

        samDictCodec.encodeHeaderLine(false);
        // SortingCollection is used to check uniqueness of sequence names
        final SortingCollection<String> sequenceNames = makeSortingCollection();

        // read reference sequence one by one and write its metadata
        while (samSequenceRecordIterator.hasNext()) {
            final SAMSequenceRecord samSequenceRecord = samSequenceRecordIterator.next();
            samDictCodec.encodeSequenceRecord(samSequenceRecord);
            sequenceNames.add(samSequenceRecord.getSequenceName());
        }

        // check uniqueness of sequences names
        final CloseableIterator<String> iterator = sequenceNames.iterator();

        if (!iterator.hasNext()) {
            return;
        }

        String current = iterator.next();
        while (iterator.hasNext()) {
            final String next = iterator.next();
            if (current.equals(next)) {
                throw new PicardException("Sequence name " + current +
                        " appears more than once in reference file");
            }
            current = next;
        }
    }

    public static SortingCollection<String> makeSortingCollection() {
        final File tmpDir = IOUtil.createTempDir("SamDictionaryNames").toFile();
        tmpDir.deleteOnExit();
        // 256 byte for one name, and 1/10 part of all memory for this, rough estimate
        long maxNamesInRam = Runtime.getRuntime().maxMemory() / 256 / 10;
        return SortingCollection.newInstance(
                String.class,
                new StringCodec(),
                String::compareTo,
                (int) Math.min(maxNamesInRam, Integer.MAX_VALUE),
                tmpDir.toPath()
        );
    }

    /**
     * Throw an exception if the two provided sequence dictionaries are not equal.
     *
     * @param firstDict first dictionary to compare
     * @param firstDictSource a user-recognizable message identifying the source of the first dictionary, preferably a file path
     * @param secondDict second dictionary to compare
     * @param secondDictSource a user-recognizable message identifying the source of the second dictionary,  preferably a file path
     */
    public static void assertSequenceDictionariesEqual(
            final SAMSequenceDictionary firstDict,
            final String firstDictSource,
            final SAMSequenceDictionary secondDict,
            final String secondDictSource) {
        try {
            SequenceUtil.assertSequenceDictionariesEqual(firstDict, secondDict);
        } catch (final SequenceUtil.SequenceListsDifferException e) {
            throw new PicardException(
                    String.format("Sequence dictionary for (%s) does not match sequence dictionary for (%s)",
                            firstDictSource,
                            secondDictSource),
                    e);
        }
    }

    private static class StringCodec implements SortingCollection.Codec<String> {
        private DataInputStream dis;
        private DataOutputStream dos;

        public StringCodec clone() {
            return new StringCodec();
        }

        public void setOutputStream(final OutputStream os) {
            dos = new DataOutputStream(os);
        }

        public void setInputStream(final InputStream is) {
            dis = new DataInputStream(is);
        }

        public void encode(final String str) {
            try {
                dos.writeUTF(str);
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        public String decode() {
            try {
                return dis.readUTF();
            } catch (EOFException e) {
                return null;
            } catch (IOException e) {
                throw new PicardException("Exception reading sequence name from temporary file.", e);
            }
        }
    }
}
