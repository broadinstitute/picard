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
 * THE SOFTWARE.
 */
package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.util.StringLineReader;

import java.io.*;
import java.util.*;

/**
 * Represents a list of intervals against a reference sequence that can be written to
 * and read from a file.  The file format is relatively simple and reflects the SAM
 * alignment format to a degree.
 *
 * A SAM style header must be present in the file which lists the sequence records
 * against which the intervals are described.  After the header the file then contains
 * records one per line in text format with the following values tab-separated:
 *   - Sequence name
 *   - Start position (1-based)
 *   - End position (1-based, end inclusive)
 *   - Strand (either + or -)
 *   - Interval name (an, ideally unique, name for the interval)
 *
 * @author Tim Fennell
 */
public class IntervalList implements Iterable<Interval> {
    private final SAMFileHeader header;
    private final List<Interval> intervals = new ArrayList<Interval>();

    private static final Log log = Log.getInstance(IntervalList.class);

    /** Constructs a new interval list using the supplied header information. */
    public IntervalList(final SAMFileHeader header) {
        if (header == null) {
            throw new IllegalArgumentException("SAMFileHeader must be supplied.");
        }
        this.header = header;
    }

    /** Gets the header (if there is one) for the interval list. */
    public SAMFileHeader getHeader() { return header; }

    /** Returns an iterator over the intervals. */
    public Iterator<Interval> iterator() { return this.intervals.iterator(); }

    /** Adds an interval to the list of intervals. */
    public void add(final Interval interval) { this.intervals.add(interval); }

    /** Sorts the internal collection of intervals by coordinate. */
    public void sort() {
        Collections.sort(this.intervals, new IntervalCoordinateComparator(this.header));
        this.header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    }

    /** Sorts and uniques the list of intervals held within this interval list. */
    public void unique() {
        sort();
        final List<Interval> tmp = getUniqueIntervals();
        this.intervals.clear();
        this.intervals.addAll(tmp);
    }

    /** Gets the set of intervals as held internally. */
    public List<Interval> getIntervals() {
        return Collections.unmodifiableList(this.intervals);
    }

    /**
     * Merges the list of intervals and then reduces them down where regions overlap
     * or are directly adjacent to one another.  During this process the "merged" interval
     * will retain the strand and name of the 5' most interval merged.
     *
     * Note: has the side-effect of sorting the stored intervals in coordinate order if not already sorted.
     *
     * @return the set of unique intervals condensed from the contained intervals
     */
    public List<Interval> getUniqueIntervals() {
        if (getHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            sort();
        }
        final List<Interval> unique = new ArrayList<Interval>();
        final ListIterator<Interval> iterator = this.intervals.listIterator();
        Interval previous = iterator.next();

        while (iterator.hasNext()) {
            final Interval next = iterator.next();
            if (previous.intersects(next) || previous.abuts(next)) {
                previous = new Interval(previous.getSequence(),
                                        previous.getStart(),
                                        Math.max(previous.getEnd(), next.getEnd()),
                                        previous.isNegativeStrand(),
                                        previous.getName());
            }
            else {
                unique.add(previous);
                previous = next;
            }
        }

        if (previous != null) unique.add(previous);

        return unique;
    }

    /** Gets the (potentially redundant) sum of the length of the intervals in the list. */
    public long getBaseCount() {
        return Interval.countBases(this.intervals);
    }

    /** Gets the count of unique bases represented by the intervals in the list. */
    public long getUniqueBaseCount() {
        return Interval.countBases(getUniqueIntervals());
    }

    /** Returns the count of intervals in the list. */
    public int size() {
        return this.intervals.size();
    }

    /**
     * Parses an interval list from a file.
     * @param file the file containing the intervals
     * @return an IntervalList object that contains the headers and intervals from the file
     */
    public static IntervalList fromFile(final File file) {
        return fromReader(new BufferedReader(new InputStreamReader(IoUtil.openFileForReading(file))));
    }

    /**
     * Parses an interval list from a reader in a stream based fashion.
     * @param in a BufferedReader that can be read from
     * @return an IntervalList object that contains the headers and intervals from the file
     */
    public static IntervalList fromReader(final BufferedReader in) {
        try {
            // Setup a reader and parse the header
            final StringBuilder builder = new StringBuilder(4096);
            String line = null;

            while ((line = in.readLine()) != null) {
                if (line.startsWith("@")) {
                    builder.append(line).append('\n');
                }
                else {
                    break;
                }
            }

            if (builder.length() == 0) {
                throw new IllegalStateException("Interval list file must contain header. ");
            }

            final StringLineReader headerReader = new StringLineReader(builder.toString());
            final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            final IntervalList list = new IntervalList(codec.decode(headerReader, "BufferedReader"));
            final SAMSequenceDictionary dict = list.getHeader().getSequenceDictionary();

            // Then read in the intervals
            final FormatUtil format = new FormatUtil();
            do {
                if (line.trim().length() == 0) continue; // skip over blank lines

                // Make sure we have the right number of fields
                final String[] fields = line.split("\t");
                if (fields.length != 5) {
                    throw new PicardException("Invalid interval record contains " +
                                              fields.length + " fields: " + line);
                }

                // Then parse them out
                final String seq = fields[0];
                final int start = format.parseInt(fields[1]);
                final int end   = format.parseInt(fields[2]);

                final boolean negative;
                if (fields[3].equals("-")) negative = true;
                else if (fields[3].equals("+")) negative = false;
                else throw new IllegalArgumentException("Invalid strand field: " + fields[3]);

                final String name = fields[4];

                final Interval interval = new Interval(seq, start, end, negative, name);
                if (dict.getSequence(seq) == null) {
                    log.warn("Ignoring interval for unknown reference: " + interval);
                }
                else {
                    list.intervals.add(interval);
                }
            }
            while ((line = in.readLine()) != null);

            return list;
        }
        catch (IOException ioe) {
            throw new PicardException("Error parsing interval list.", ioe);
        }
        finally {
            try { in.close(); } catch (Exception e) { /* do nothing */ }
        }
    }

    /**
     * Writes out the list of intervals to the supplied file.
     * @param file a file to write to.  If exists it will be overwritten.
     */
    public void write(final File file) {
        try {
            final BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IoUtil.openFileForWriting(file)));
            final FormatUtil format = new FormatUtil();

            // Write out the header
            if (this.header != null) {
                final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                codec.encode(out, this.header);
            }

            // Write out the intervals
            for (final Interval interval : this) {
                out.write(interval.getSequence());
                out.write('\t');
                out.write(format.format(interval.getStart()));
                out.write('\t');
                out.write(format.format(interval.getEnd()));
                out.write('\t');
                out.write(interval.isPositiveStrand() ? '+' : '-');
                out.write('\t');
                out.write(interval.getName());
                out.newLine();
            }

            out.flush();
            out.close();
        }
        catch (IOException ioe) {
            throw new PicardException("Error writing out interval list to file: " + file.getAbsolutePath(), ioe);
        }
    }
}

/**
 * Comparator that orders intervals based on their sequence index, by coordinate
 * then by strand and finally by name.
 */
class IntervalCoordinateComparator implements Comparator<Interval> {
    private final SAMFileHeader header;

    /** Constructs a comparator using the supplied sequence header. */
    IntervalCoordinateComparator(final SAMFileHeader header) {
        this.header = header;
    }

    public int compare(final Interval lhs, final Interval rhs) {
        final int lhsIndex = this.header.getSequenceIndex(lhs.getSequence());
        final int rhsIndex = this.header.getSequenceIndex(rhs.getSequence());
        int retval = lhsIndex - rhsIndex;

        if (retval == 0) retval = lhs.getStart() - rhs.getStart();
        if (retval == 0) retval = lhs.getEnd()   - rhs.getEnd();
        if (retval == 0) {
            if (lhs.isPositiveStrand() && rhs.isNegativeStrand()) retval = -1;
            else if (lhs.isNegativeStrand() && rhs.isPositiveStrand()) retval = 1;
        }
        if (retval == 0) {
            retval = lhs.getName().compareTo(rhs.getName());
        }

        return retval;
    }
}