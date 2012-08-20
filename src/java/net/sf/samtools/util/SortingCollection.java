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

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Collection to which many records can be added.  After all records are added, the collection can be
 * iterated, and the records will be returned in order defined by the comparator.  Records may be spilled
 * to a temporary directory if there are more records added than will fit in memory.  As a result of this,
 * the objects returned may not be identical to the objects added to the collection, but they should be
 * equal as determined by the codec used to write them to disk and read them back.
 *
 * When iterating over the collection, the number of file handles required is numRecordsInCollection/maxRecordsInRam.
 * If this becomes a limiting factor, a file handle cache could be added.
 *
 * If Snappy DLL is available and snappy.disable system property is not set to true, then Snappy is used
 * to compress temporary files.
 */
public class SortingCollection<T> implements Iterable<T> {

    /**
     * Client must implement this class, which defines the way in which records are written to and
     * read from file.
     */
    public interface Codec<T> extends Cloneable {
        /**
         * Where to write encoded output
         * @param os
         */
        void setOutputStream(OutputStream os);

        /**
         * Where to read encoded input from
         * @param is
         */
        void setInputStream(InputStream is);
        /**
         * Write object to output stream
         * @param val what to write
         */
        void encode(T val);

        /**
         * Read the next record from the input stream and convert into a java object.
         * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
         * a record.
         */
        T decode();

        /**
         * Must return a cloned copy of the codec that can be used independently of
         * the original instance.  This is required so that multiple codecs can exist simultaneously
         * that each is reading a separate file.
         */
        Codec<T> clone();
    }

    /** Directories where files of sorted records go. */
    private final File[] tmpDirs;

    /** The minimum amount of space free on a temp filesystem to write a file there. */
    private final long TMP_SPACE_FREE = IOUtil.FIVE_GBS;

    /**
     * Used to write records to file, and used as a prototype to create codecs for reading.
     */
    private final SortingCollection.Codec<T> codec;

    /**
     * For sorting, both when spilling records to file, and merge sorting.
     */
    private final Comparator<T> comparator;
    private final int maxRecordsInRam;
    private int numRecordsInRam = 0;
    private T[] ramRecords;
    private boolean iterationStarted = false;
    private boolean doneAdding = false;

    /**
     * Set to true when all temp files have been cleaned up
     */
    private boolean cleanedUp = false;

    /**
     * List of files in tmpDir containing sorted records
     */
    private final List<File> files = new ArrayList<File>();

    private boolean destructiveIteration = true;

    private TempStreamFactory tempStreamFactory = new TempStreamFactory();

    /**
     * Prepare to accumulate records to be sorted
     * @param componentType Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec For writing records to file and reading them back into RAM
     * @param comparator Defines output sort order
     * @param maxRecordsInRam how many records to accumulate before spilling to disk
     * @param tmpDir Where to write files of records that will not fit in RAM
     */
    private SortingCollection(final Class<T> componentType, final SortingCollection.Codec<T> codec,
                             final Comparator<T> comparator, final int maxRecordsInRam, final File... tmpDir) {
        if (maxRecordsInRam <= 0) {
            throw new IllegalArgumentException("maxRecordsInRam must be > 0");
        }

        if (tmpDir == null || tmpDir.length == 0) {
            throw new IllegalArgumentException("At least one temp directory must be provided.");
        }

        this.tmpDirs = tmpDir;
        this.codec = codec;
        this.comparator = comparator;
        this.maxRecordsInRam = maxRecordsInRam;
        this.ramRecords = (T[])Array.newInstance(componentType, maxRecordsInRam);
    }

    public void add(final T rec) {
        if (doneAdding) {
            throw new IllegalStateException("Cannot add after calling doneAdding()");
        }
        if (iterationStarted) {
            throw new IllegalStateException("Cannot add after calling iterator()");
        }
        if (numRecordsInRam == maxRecordsInRam) {
            spillToDisk();
        }
        ramRecords[numRecordsInRam++] = rec;
    }

    /**
     * This method can be called after caller is done adding to collection, in order to possibly free
     * up memory.  If iterator() is called immediately after caller is done adding, this is not necessary,
     * because iterator() triggers the same freeing.
     */
    public void doneAdding() {
        if (this.cleanedUp) {
            throw new IllegalStateException("Cannot call doneAdding() after cleanup() was called.");
        }
        if (doneAdding) {
            return;
        }

        doneAdding = true;

        if (this.files.isEmpty()) {
            return;
        }

        if (this.numRecordsInRam > 0) {
            spillToDisk();
        }

        // Facilitate GC
        this.ramRecords = null;
    }

    /**
     * @return True if this collection is allowed to discard data during iteration in order to reduce memory
     * footprint, precluding a second iteration over the collection.
     */
    public boolean isDestructiveIteration() {
        return destructiveIteration;
    }

    /**
     * Tell this collection that it is allowed to discard data during iteration in order to reduce memory footprint,
     * precluding a second iteration.  This is true by default.
     */
    public void setDestructiveIteration(boolean destructiveIteration) {
        this.destructiveIteration = destructiveIteration;
    }

    /**
     * Sort the records in memory, write them to a file, and clear the buffer of records in memory.
     */
    private void spillToDisk() {
        try {
            Arrays.sort(this.ramRecords, 0, this.numRecordsInRam, this.comparator);
            final File f = newTempFile();
            OutputStream os = null;
            try {
                os = tempStreamFactory.wrapTempOutputStream(new FileOutputStream(f), IOUtil.STANDARD_BUFFER_SIZE);
                this.codec.setOutputStream(os);
                for (int i = 0; i < this.numRecordsInRam; ++i) {
                    this.codec.encode(ramRecords[i]);
                    // Facilitate GC
                    this.ramRecords[i] = null;
                }

                os.flush();
            }
            finally {
                if (os != null) {
                    os.close();
                }
            }

            this.numRecordsInRam = 0;
            this.files.add(f);

        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Creates a new tmp file on one of the available temp filesystems, registers it for deletion
     * on JVM exit and then returns it.
     */
    private File newTempFile() throws IOException {
        return IOUtil.newTempFile("sortingcollection.", ".tmp", this.tmpDirs, TMP_SPACE_FREE);
    }

    /**
     * Prepare to iterate through the records in order.  This method may be called more than once,
     * but add() may not be called after this method has been called.
     */
    public CloseableIterator<T> iterator() {
        if (this.cleanedUp) {
            throw new IllegalStateException("Cannot call iterator() after cleanup() was called.");
        }
        doneAdding();

        this.iterationStarted = true;
        if (this.files.isEmpty()) {
            return new InMemoryIterator();
        } else {
            return new MergingIterator();
        }
    }

    /**
     * Delete any temporary files.  After this method is called, iterator() may not be called.
     */
    public void cleanup() {
        this.iterationStarted = true;
        this.cleanedUp = true;

        IOUtil.deleteFiles(this.files);
    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec For writing records to file and reading them back into RAM
     * @param comparator Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     * @param tmpDir Where to write files of records that will not fit in RAM
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final File... tmpDir) {
        return new SortingCollection<T>(componentType, codec, comparator, maxRecordsInRAM, tmpDir);

    }

    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters
     *
     * @param componentType Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec For writing records to file and reading them back into RAM
     * @param comparator Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     * @param tmpDirs Where to write files of records that will not fit in RAM
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM,
                                                       final Collection<File> tmpDirs) {
        return new SortingCollection<T>(componentType,
                                        codec,
                                        comparator,
                                        maxRecordsInRAM,
                                        tmpDirs.toArray(new File[tmpDirs.size()]));

    }


    /**
     * Syntactic sugar around the ctor, to save some typing of type parameters.  Writes files to java.io.tmpdir
     *
     * @param componentType Class of the record to be sorted.  Necessary because of Java generic lameness.
     * @param codec For writing records to file and reading them back into RAM
     * @param comparator Defines output sort order
     * @param maxRecordsInRAM how many records to accumulate in memory before spilling to disk
     */
    public static <T> SortingCollection<T> newInstance(final Class<T> componentType,
                                                       final SortingCollection.Codec<T> codec,
                                                       final Comparator<T> comparator,
                                                       final int maxRecordsInRAM) {

        final File tmpDir = new File(System.getProperty("java.io.tmpdir"));
        return new SortingCollection<T>(componentType, codec, comparator, maxRecordsInRAM, tmpDir);
    }

    /**
     * For iteration when number of records added is less than the threshold for spilling to disk.
     */
    class InMemoryIterator implements CloseableIterator<T> {
        private int iterationIndex = 0;

        InMemoryIterator() {
            Arrays.sort(SortingCollection.this.ramRecords,
                        0,
                        SortingCollection.this.numRecordsInRam,
                        SortingCollection.this.comparator);
        }

        public void close() {
            // nothing to do
        }

        public boolean hasNext() {
            return this.iterationIndex < SortingCollection.this.numRecordsInRam;
        }

        public T next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            T ret = SortingCollection.this.ramRecords[iterationIndex];
            if (destructiveIteration) SortingCollection.this.ramRecords[iterationIndex] = null;
            ++iterationIndex;
            return ret;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    /**
     * For iteration when spilling to disk has occurred.
     * Each file is has records in sort order within the file.
     * This iterator automatically closes when it iterates to the end, but if not iterating
     * to the end it is a good idea to call close().
     *
     * Algorithm: MergingIterator maintains a PriorityQueue of PeekFileRecordIterators.
     * Each PeekFileRecordIterator iterates through a file in which the records are sorted.
     * The comparator for PeekFileRecordIterator used by the PriorityQueue peeks at the next record from
     * the file, so the first element in the PriorityQueue is the file that has the next record to be emitted.
     * In order to get the next record, the first PeekFileRecordIterator in the PriorityQueue is popped,
     * the record is obtained from that iterator, and then if that iterator is not empty, it is pushed back into
     * the PriorityQueue.  Because it now has a different record as its next element, it may go into another
     * location in the PriorityQueue
     */
    class MergingIterator implements CloseableIterator<T> {
        private final PollableTreeSet<PeekFileRecordIterator> queue;

        MergingIterator() {
            this.queue = new PollableTreeSet<PeekFileRecordIterator>(new PeekFileRecordIteratorComparator());
            int n = 0;
            for (final File f : SortingCollection.this.files) {
                final FileRecordIterator it = new FileRecordIterator(f);
                if (it.hasNext()) {
                    this.queue.add(new PeekFileRecordIterator(it, n++));
                }
                else {
                    it.close();
                }
            }
        }

        public boolean hasNext() {
            return !this.queue.isEmpty();
        }

        public T next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            final PeekFileRecordIterator fileIterator = queue.poll();
            final T ret = fileIterator.next();
            if (fileIterator.hasNext()) {
                this.queue.add(fileIterator);
            }
            else {
                ((CloseableIterator<T>)fileIterator.getUnderlyingIterator()).close();
            }

            return ret;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }

        public void close() {
            while (!this.queue.isEmpty()) {
                final PeekFileRecordIterator it = this.queue.poll();
                ((CloseableIterator<T>)it.getUnderlyingIterator()).close();
            }
        }
    }

    /**
     * Read a file of records in format defined by the codec
     */
    class FileRecordIterator implements CloseableIterator<T> {
        private final File file;
        private final FileInputStream is;
        private final Codec<T> codec;
        private T currentRecord = null;

        FileRecordIterator(final File file) {
            this.file = file;
            try {
                this.is = new FileInputStream(file);
                this.codec = SortingCollection.this.codec.clone();
                this.codec.setInputStream(tempStreamFactory.wrapTempInputStream(this.is, IOUtil.STANDARD_BUFFER_SIZE));
                advance();
            }
            catch (FileNotFoundException e) {
                throw new RuntimeIOException(e);
            }
        }

        public boolean hasNext() {
            return this.currentRecord != null;
        }

        public T next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            final T ret = this.currentRecord;
            advance();
            return ret;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }

        private void advance() {
            this.currentRecord = this.codec.decode();
        }

        public void close() {
            CloserUtil.close(this.is);
        }
    }


    /**
     * Just a typedef
     */
    class PeekFileRecordIterator extends PeekIterator<T> {
        final int n; // A serial number used for tie-breaking in the sort
        PeekFileRecordIterator(final Iterator<T> underlyingIterator, final int n) {
            super(underlyingIterator);
            this.n = n;
        }
    }

    class PeekFileRecordIteratorComparator implements Comparator<PeekFileRecordIterator> {

        public int compare(final PeekFileRecordIterator lhs, final PeekFileRecordIterator rhs) {
            final int result = comparator.compare(lhs.peek(), rhs.peek());
            if (result == 0) return lhs.n - rhs.n;
            else return result;
        }
    }

    /** Little class that provides the Java 1.5 TreeSet with a poll() method */
    static class PollableTreeSet<T> extends TreeSet<T> {
        PollableTreeSet(final Comparator<? super T> comparator) {
            super(comparator);
        }

        public T poll() {
            if (isEmpty()) {
                return null;
            }
            else {
                final T t = first();
                remove(t);
                return t;
            }
        }
    }
}
