/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.tribble.index;

import org.broad.tribble.TribbleException;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * User: jrobinso
 * <p/>
 * An abstract implementation of the index class.  This class takes care of the basics that are common
 * to all of the current indexing classes; including the version information, common header properties,
 * and reading and writing the header to disk.
 */
public abstract class AbstractIndex implements Index {

    // todo -- up to version 4 and use ETag to detect out of date
    // todo -- inode number + size in bytes + modification time
    // todo -- remove MD5

    // the current version of the index
    public static int VERSION = 3;

    private final static String NO_MD5 = "";
    private final static long NO_FILE_SIZE = -1L;
    private final static long NO_TS = -1L;

    protected int version;                    // Our version value
    protected File indexedFile = null;         // The file we've created this index for
    protected long indexedFileSize = NO_FILE_SIZE; // The size of the indexed file
    protected long indexedFileTS = NO_TS;      // The timestamp
    protected String indexedFileMD5 = NO_MD5;        // The MD5 value, generally not filled in (expensive to calc)
    protected int flags;

    public boolean hasFileSize() {
        return indexedFileSize != NO_FILE_SIZE;
    }

    public boolean hasTimestamp() {
        return indexedFileTS != NO_TS;
    }

    public boolean hasMD5() {
        return indexedFileMD5 != NO_MD5;
    }

    // out listing of properties, order preserved
    private LinkedHashMap<String, String> properties;

    // the hashmap of our chromosome bins
    protected LinkedHashMap<String, ChrIndex> chrIndices;

    // any flags we're using:
    private static final int SEQUENCE_DICTIONARY_FLAG = 0x8000; // if we have a sequence dictionary in our header

    /**
     * Returns true if this and obj are 'effectively' equivalent data structures.
     *
     * @param obj
     * @return
     */
    public boolean equalsIgnoreProperties(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof AbstractIndex)) {
            System.err.printf("equals: %s not instance of AbstractIndex", obj);
            return false;
        }

        AbstractIndex other = (AbstractIndex) obj;

        if (version != other.version) {
            System.err.printf("equals version: this %d != other %d%n", version, other.version);
            return false;
        }

        if (indexedFile != other.indexedFile && (indexedFile == null || !indexedFile.equals(other.indexedFile))) {
            System.err.printf("equals indexedFile: this %s != other %s%n", indexedFile, other.indexedFile);
            return false;
        }

        if (indexedFileSize != other.indexedFileSize) {
            System.err.printf("equals indexedFileSize: this %d != other %d%n", indexedFileSize, other.indexedFileSize);
            return false;
        }

        if (!indexedFileMD5.equals(other.indexedFileMD5)) {
            System.err.printf("equals indexedFileMD5: this %s != other %s%n", indexedFileMD5, other.indexedFileMD5);
            return false;
        }

        if (flags != other.flags) {
            System.err.printf("equals flags: this %d != other %d%n", flags, other.flags);
            return false;
        }

        if (!chrIndices.equals(other.chrIndices)) {
            System.err.printf("equals chrIndeces: this %s != other %s%n", chrIndices, other.chrIndices);
            return false;
        }

        return true;
    }

    /**
     * create an abstract index, with defaults for the version value, and empty properties and chromosome lists
     */
    public AbstractIndex() {
        this.version = VERSION; // <= is overriden when file is read
        this.properties = new LinkedHashMap<String, String>();
        chrIndices = new LinkedHashMap();
    }

    /**
     * create an index file from the target feature file
     *
     * @param featureFile the feature file to create an index from
     */
    public AbstractIndex(String featureFile) {
        this(new File(featureFile));
    }

    public AbstractIndex(File featureFile) {
        this();
        this.indexedFile = featureFile;
    }

    public AbstractIndex(AbstractIndex parent) {
        this();
        this.version = parent.version;
        this.indexedFile = parent.indexedFile;
        this.indexedFileSize = parent.indexedFileSize;
        this.indexedFileTS = parent.indexedFileTS;
        this.indexedFileMD5 = parent.indexedFileMD5;
        this.flags = parent.flags;
        this.properties = (LinkedHashMap<String, String>) parent.properties.clone();
    }


    /**
     * check the current version against the version we read in
     *
     * @return true if we're up to date, false otherwise
     */
    public boolean isCurrentVersion() {
        return version == VERSION;
    }

    public File getIndexedFile() {
        return indexedFile;
    }

    public long getIndexedFileSize() {
        return indexedFileSize;
    }

    public long getIndexedFileTS() {
        return indexedFileTS;
    }

    public String getIndexedFileMD5() {
        return indexedFileMD5;
    }

    public int getFlags() {
        return flags;
    }

    public int getVersion() {
        return version;
    }

    /**
     * set the MD5 value
     *
     * @param md5 the MD5 sum
     */
    public void setMD5(String md5) {
        this.indexedFileMD5 = md5;
    }

    /**
     * do we have an entry for the target chromosome?
     *
     * @param chr the chromosome (or contig) name
     * @return true if we have an entry; false otherwise
     */
    public boolean containsChromosome(String chr) {
        return chrIndices.containsKey(chr);
    }

    /**
     * Should be called after the index is created to finalize all of the header information
     */
    public void finalizeIndex() {
        // these two functions must be called now because the file may be being written during on the fly indexing
        if (indexedFile != null) {
            this.indexedFileSize = indexedFile.length();
            this.indexedFileTS = indexedFile.lastModified();
        }
    }

    /**
     * write the header to the target output stream
     *
     * @param dos the little endian output stream
     * @throws IOException an exception when we can't write to the file
     */
    private void writeHeader(LittleEndianOutputStream dos) throws IOException {
        int magicNumber = 1480870228;   //  byte[]{'T', 'I', 'D', 'X'};

        dos.writeInt(magicNumber);
        dos.writeInt(getType());
        dos.writeInt(version);
        dos.writeString(indexedFile.getAbsolutePath());
        dos.writeLong(indexedFileSize);
        dos.writeLong(indexedFileTS);
        dos.writeString(indexedFileMD5);
        dos.writeInt(flags);

        // Properties (Version 3 and later)
        dos.writeInt(properties.size());
        for (Map.Entry<String, String> prop : properties.entrySet()) {
            dos.writeString(prop.getKey());
            dos.writeString(prop.getValue());
        }
    }

    /**
     * read the header from the input stream
     *
     * @param dis the little endian input stream
     * @throws IOException if we fail to read from the file at any point
     */
    private void readHeader(LittleEndianInputStream dis) throws IOException {

        version = dis.readInt();
        indexedFile = new File(dis.readString());
        indexedFileSize = dis.readLong();
        indexedFileTS = dis.readLong();
        indexedFileMD5 = dis.readString();
        flags = dis.readInt();
        if (version < 3 && (flags & SEQUENCE_DICTIONARY_FLAG) == SEQUENCE_DICTIONARY_FLAG) {
            readSequenceDictionary(dis);
        }

        if (version >= 3) {
            int nProperties = dis.readInt();
            while (nProperties-- > 0) {
                String key = dis.readString();
                String value = dis.readString();
                properties.put(key, value);
            }
        }
    }

    /**
     * Kept to maintain backward compatibility with pre version 3 indexes.  The sequence dictionary is no longer
     * used,  use getSequenceNames() instead.
     *
     * @param dis
     * @throws IOException
     */
    private void readSequenceDictionary(LittleEndianInputStream dis) throws IOException {
        int size = dis.readInt();
        if (size < 0) throw new IllegalStateException("Size of the sequence dictionary entries is negative");
        for (int x = 0; x < size; x++) {
            dis.readString();
            dis.readInt();
        }
    }


    /**
     * get the sequence names
     *
     * @return a linked hash set of sequence names (Linked to ensure ordering)
     */
    public LinkedHashSet<String> getSequenceNames() {
        return new LinkedHashSet(chrIndices.keySet());
    }

    /**
     * @param chr   the chromosome
     * @param start the start position, one based
     * @param end   the end position, one based
     * @return a list of blocks that include the region defined from start to stop.  Can never return null
     */
    public List<Block> getBlocks(String chr, int start, int end) {
        return getChrIndex(chr).getBlocks(start, end);
    }

    public List<Block> getBlocks(String chr) {
        return getChrIndex(chr).getBlocks();
    }

    /**
     * @param chr
     * @return return the ChrIndex associated with chr, or throw an IllegalArgumentException if not found
     */
    private final ChrIndex getChrIndex(final String chr) {
        ChrIndex chrIdx = chrIndices.get(chr);
        if (chrIdx == null) {
            throw new IllegalArgumentException("getBlocks() called with of unknown contig " + chr);
        } else {
            return chrIdx;
        }
    }

    /**
     * Store self to a file
     *
     * @param stream the input stream
     * @throws java.io.IOException
     */
    public void write(LittleEndianOutputStream stream) throws IOException {
        writeHeader(stream);

        //# of chromosomes
        stream.writeInt(chrIndices.size());
        for (ChrIndex chrIdx : chrIndices.values()) {
            chrIdx.write(stream);
        }
    }

    /**
     * read the index from a little endian input stream
     *
     * @param dis the little endian input stream
     * @throws IOException if we fail to successfully read from disk
     */
    public void read(LittleEndianInputStream dis) throws IOException {
        try {
            readHeader(dis);

            int nChromosomes = dis.readInt();
            chrIndices = new LinkedHashMap<String, ChrIndex>(nChromosomes);

            while (nChromosomes-- > 0) {
                ChrIndex chrIdx = (ChrIndex) getChrIndexClass().newInstance();
                chrIdx.read(dis);
                chrIndices.put(chrIdx.getName(), chrIdx);
            }

        } catch (InstantiationException e) {
            throw new TribbleException.UnableToCreateCorrectIndexType("Unable to create class " + getChrIndexClass(), e);
        } catch (IllegalAccessException e) {
            throw new TribbleException.UnableToCreateCorrectIndexType("Unable to create class " + getChrIndexClass(), e);
        } finally {
            dis.close();
        }

        //printIndexInfo();
    }

    protected void printIndexInfo() {
        System.out.println(String.format("Index for %s with %d indices", indexedFile, chrIndices.size()));
        BlockStats stats = getBlockStats(true);
        System.out.println(String.format("  total blocks %d", stats.total));
        System.out.println(String.format("  total empty blocks %d", stats.empty));
    }

    protected static class BlockStats {
        long total = 0, empty = 0, objects = 0, size = 0;
    }

    protected BlockStats getBlockStats(boolean logDetails) {
        BlockStats stats = new BlockStats();
        for (Map.Entry<String, ChrIndex> elt : chrIndices.entrySet()) {
            List<Block> blocks = elt.getValue().getBlocks();

            if (blocks != null) {
                int nBlocks = blocks.size();

                int nEmptyBlocks = 0;
                for (Block b : elt.getValue().getBlocks()) {
                    if (b.getSize() == 0) nEmptyBlocks++;
                }
                stats.empty += nEmptyBlocks;
                stats.total += nBlocks;

                if (logDetails)
                    System.out.println(String.format("  %s => %d blocks, %d empty, %.2f", elt.getKey(), nBlocks, nEmptyBlocks, (100.0 * nEmptyBlocks) / nBlocks));
            }
        }

        return stats;
    }

    protected String statsSummary() {
        BlockStats stats = getBlockStats(false);
        return String.format("%12d blocks (%12d empty (%.2f%%))", stats.total, stats.empty, (100.0 * stats.empty) / stats.total);
    }

    /**
     * add properties to the the existing property listing;
     *
     * @param key   the key
     * @param value the value, stored as a string, though it may represent an different underlying type
     */
    public void addProperty(String key, String value) {
        properties.put(key, value);
    }

    /**
     * return a mapping of name to property value
     *
     * @return the mapping of values as an unmodifiable map
     */
    public Map<String, String> getProperties() {
        return Collections.unmodifiableMap(properties);
    }

    /**
     * get the index type
     *
     * @return
     */
    protected abstract int getType();

    /**
     * returns the class for the index type
     *
     * @return a Class, from which a new instance can be created
     */
    public abstract Class getChrIndexClass();
}
