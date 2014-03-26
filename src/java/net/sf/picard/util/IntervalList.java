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
import net.sf.samtools.*;
import net.sf.samtools.util.CollectionUtil;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringLineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFFileReader;
import net.sf.samtools.util.StringUtil;

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
 *    Sequence name,
 *    Start position (1-based),
 *    End position (1-based, end inclusive),
 *    Strand (either + or -),
 *    Interval name (an, ideally unique, name for the interval),
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
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

    /** Adds a Collection of intervals to the list of intervals. */
    public void addall(final Collection<Interval> intervals) {
        this.intervals.addAll(intervals);
    }

    /** Sorts the internal collection of intervals by coordinate. */
    @Deprecated // Use sorted() instead of sort(). The sort() function modifies the object in-place and
    // is therefore difficult to work with. sorted() returns a new IntervalList that is sorted
    public void sort() {
        Collections.sort(this.intervals, new IntervalCoordinateComparator(this.header));
        this.header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    }

    /** returns an independent sorted IntervalList*/
    public IntervalList sorted() {
        final IntervalList sorted = IntervalList.copyOf(this);
        Collections.sort(sorted.intervals, new IntervalCoordinateComparator(sorted.header));
        sorted.header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return sorted;
    }

    /** Returned an independent IntervalList that is sorted and uniquified. */
    public IntervalList uniqued() {
        return uniqued(true);
    }

    /**
     * Returned an independent IntervalList that is sorted and uniquified.
     * @param concatenateNames If false, interval names are not concatenated when merging intervals to save space.
     */
    public IntervalList uniqued(final boolean concatenateNames) {
        final List<Interval> tmp = getUniqueIntervals(sorted(), concatenateNames);
        final IntervalList value = new IntervalList(this.header.clone());
        value.intervals.addAll(tmp);
        return value;
    }

    /** Sorts and uniques the list of intervals held within this interval list. */
    @Deprecated//use uniqued() instead. This function modifies the object in-place and
    // is therefore difficult to work with.
    public void unique() {
        unique(true);
    }

    /**
     * Sorts and uniques the list of intervals held within this interval list.
     * @param concatenateNames If false, interval names are not concatenated when merging intervals to save space.
     */
    @Deprecated//use uniqued() instead. This function modifies the object in-place and
    // is therefore difficult to work with
    public void unique(final boolean concatenateNames) {
        sort();
        final List<Interval> tmp = getUniqueIntervals(concatenateNames);
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
    @Deprecated//use uniqued().getIntervals() instead. This function modifies the object in-place and
    // is therefore difficult to work with
    public List<Interval> getUniqueIntervals() {
        return getUniqueIntervals(true);
    }

    //NO SIDE EFFECTS HERE!
    /**
     * Merges list of intervals and reduces them like net.sf.picard.util.IntervalList#getUniqueIntervals()
     * @param concatenateNames If false, the merged interval has the name of the earlier interval.  This keeps name shorter.
     */
    public static List<Interval> getUniqueIntervals(final IntervalList list, final boolean concatenateNames) {

        final List<Interval> intervals;
        if (list.getHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            intervals = list.sorted().intervals;
        }
        else {
            intervals = list.intervals;
        }

        final List<Interval> unique = new ArrayList<Interval>();
        final TreeSet<Interval> toBeMerged = new TreeSet<Interval>();
        Interval current = null;

        for (final Interval next : intervals) {
            if (current == null) {
                toBeMerged.add(next);
                current = next;
            }
            else if (current.intersects(next) || current.abuts(next)) {
                toBeMerged.add(next);
                current = new Interval(current.getSequence(), current.getStart(), Math.max(current.getEnd(), next.getEnd()), false , "");
            }
            else {
                // Emit merged/unique interval
                unique.add(merge(toBeMerged, concatenateNames));

                // Set current == next for next iteration
                toBeMerged.clear();
                current = next;
                toBeMerged.add(current);
            }
        }

        if (toBeMerged.size() > 0) unique.add(merge(toBeMerged, concatenateNames));
        return unique;
    }

    /**
     * Merges list of intervals and reduces them like net.sf.picard.util.IntervalList#getUniqueIntervals()
     * @param concatenateNames If false, the merged interval has the name of the earlier interval.  This keeps name shorter.
     */
    @Deprecated //use uniqued(concatenateNames).getIntervals() or the static version instead to avoid changing the underlying object.
    /**
     * Merges list of intervals and reduces them like net.sf.picard.util.IntervalList#getUniqueIntervals()
     * @param concatenateNames If false, the merged interval has the name of the earlier interval.  This keeps name shorter.
     */
    public List<Interval> getUniqueIntervals(final boolean concatenateNames) {
        if (getHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            sort();
        }

        return getUniqueIntervals(this, concatenateNames);
    }

    /** Merges a sorted collection of intervals and optionally concatenates unique names or takes the first name. */
    static Interval merge(final SortedSet<Interval> intervals, final boolean concatenateNames) {
        final String chrom = intervals.first().getSequence();
        int start = intervals.first().getStart();
        int end   = intervals.last().getEnd();
        final boolean neg  = intervals.first().isNegativeStrand();
        final LinkedHashSet<String> names = new LinkedHashSet<String>();
        final String name;

        for (final Interval i : intervals) {
            if (i.getName() != null) names.add(i.getName());
            start = Math.min(start, i.getStart());
            end   = Math.max(end, i.getEnd());
        }

        if (concatenateNames) { name = StringUtil.join("|", names); }
        else { name = names.iterator().next(); }

        return new Interval(chrom, start, end, neg, name);
    }

    /** Gets the (potentially redundant) sum of the length of the intervals in the list. */
    public long getBaseCount() {
        return Interval.countBases(this.intervals);
    }

    /** Gets the count of unique bases represented by the intervals in the list. */
    public long getUniqueBaseCount() {
        return uniqued().getBaseCount();
    }

    /** Returns the count of intervals in the list. */
    public int size() {
        return this.intervals.size();
    }

    /**
     * Parse a VCF file and convert to an IntervalList The name field of the IntervalList is taken from the ID field of the variant, if it exists. if not,
     * creates a name of the format interval-n where n is a running number that increments only on un-named intervals
     * @param file
     * @return
     */
    public static IntervalList fromVcf(final File file){
        final VCFFileReader vcfFileReader = new VCFFileReader(file, false);
        final IntervalList intervalList = IntervalList.fromVcf(vcfFileReader);
        vcfFileReader.close();
        return intervalList;
    }

    /**
     * Converts a vcf to an IntervalList. The name field of the IntervalList is taken from the ID field of the variant, if it exists. if not,
     * creates a name of the format interval-n where n is a running number that increments only on un-named intervals
     * @param vcf the vcfReader to be used for the conversion
     * @return an IntervalList constructed from input vcf
     */
    public static IntervalList fromVcf(final VCFFileReader vcf){

        //grab the dictionary from the VCF and use it in the IntervalList
        final SAMSequenceDictionary dict = vcf.getFileHeader().getSequenceDictionary();
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        samFileHeader.setSequenceDictionary(dict);
        final IntervalList list = new IntervalList( samFileHeader);

        int intervals=0;
        for(final VariantContext vc : vcf){
            if(!vc.isFiltered()){
                String name = vc.getID();
                if(".".equals(name) || name == null)
                    name = "interval-" + (++intervals);
                list.add(new Interval(vc.getChr(), vc.getStart(), vc.getEnd(), false, name));
            }
        }

        return list;
    }

    /** creates a independent copy of the given IntervalList
     *
     * @param list
     * @return
     */
    public static IntervalList copyOf(final IntervalList list){
        final IntervalList clone = new IntervalList(list.header.clone());
        clone.intervals.addAll(list.intervals);
        return clone;
    }

    /**
     * Parses an interval list from a file.
     * @param file the file containing the intervals
     * @return an IntervalList object that contains the headers and intervals from the file
     */
    public static IntervalList fromFile(final File file) {
        final BufferedReader reader=IoUtil.openFileForBufferedReading(file);
        final IntervalList list = fromReader(reader);
        try {
            reader.close();
        } catch (final IOException e) {
            throw new PicardException(String.format("Failed to close file %s after reading",file));
        }

        return list;
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

            //there might not be any lines after the header, in which case we should return an empty list
            if(line == null) return list;

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
        catch (final IOException ioe) {
            throw new PicardException("Error parsing interval list.", ioe);
        }
        finally {
            try { in.close(); } catch (final Exception e) { /* do nothing */ }
        }
    }

    /**
     * Writes out the list of intervals to the supplied file.
     * @param file a file to write to.  If exists it will be overwritten.
     */
    public void write(final File file) {
        try {
            final BufferedWriter out = IoUtil.openFileForBufferedWriting(file);
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
                if(interval.getName() != null){
                    out.write(interval.getName());
                }
                else{
                    out.write(".");
                }
                out.newLine();
            }

            out.flush();
            out.close();
        }
        catch (final IOException ioe) {
            throw new PicardException("Error writing out interval list to file: " + file.getAbsolutePath(), ioe);
        }
    }

    /**
     * A utility function for generating the intersection of two IntervalLists, checks for equal dictionaries.
     *
     * @param list1 the first IntervalList
     * @param list2 the second IntervalList
     * @return the intersection of list1 and list2.
     */

    public static IntervalList intersection(final IntervalList list1, final IntervalList list2) {

        final IntervalList result;
        // Ensure that all the sequence dictionaries agree and merge the lists
        SequenceUtil.assertSequenceDictionariesEqual(list1.getHeader().getSequenceDictionary(),
                list2.getHeader().getSequenceDictionary());

        result = new IntervalList(list1.getHeader().clone());

        final OverlapDetector<Interval> detector = new OverlapDetector<Interval>(0, 0);

        detector.addAll(list1.getIntervals(), list1.getIntervals());

        for (final Interval i : list2.getIntervals()) {
            final Collection<Interval> as = detector.getOverlaps(i);
            for (final Interval j : as) {
                final Interval tmp = i.intersect(j);

                result.add(tmp);
            }
        }
        return result.uniqued();

    }

    /**
     * A utility function for intersecting a list of IntervalLists, checks for equal dictionaries.
     *
     * @param lists the list of IntervalList
     * @return the intersection of all the IntervalLists in lists.
     */


    public static IntervalList intersection(final Collection<IntervalList> lists) {

        IntervalList intersection = null;
        for (final IntervalList list : lists) {
            if(intersection == null){
                intersection = list;
            }
            else{
                intersection = intersection(intersection, list);
            }
        }
        return intersection;
    }




    /**
     * A utility function for merging a list of IntervalLists, checks for equal dictionaries.
     * Merging does not look for overlapping intervals nor uniquify
     *
     * @param lists a list of IntervalList
     * @return the union of all the IntervalLists in lists.
     */
    public static IntervalList concatenate(final Collection<IntervalList> lists) {
        if(lists.isEmpty()){
            throw new PicardException("Cannot concatenate an empty list of IntervalLists.");
        }

        // Ensure that all the sequence dictionaries agree and merge the lists
        final SAMFileHeader header = lists.iterator().next().getHeader().clone();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        final IntervalList merged = new IntervalList(header);

        for (final IntervalList in : lists) {
            SequenceUtil.assertSequenceDictionariesEqual(merged.getHeader().getSequenceDictionary(),
                    in.getHeader().getSequenceDictionary());

            merged.addall(in.intervals);
            }
        return merged;
    }



    /**
     * A utility function for finding the union of a list of IntervalLists, checks for equal dictionaries.
     * also looks for overlapping intervals, uniquifies, and sorts (by coordinate)
     *
     * @param lists the list of IntervalList
     * @return the union of all the IntervalLists in lists.
     */
    public static IntervalList union(final Collection<IntervalList> lists) {
        final IntervalList merged = concatenate(lists);
        return merged.uniqued();
    }


    public static IntervalList union(final IntervalList list1, final IntervalList list2) {
        final Collection<IntervalList> duo = CollectionUtil.makeList(list1, list2);
        return IntervalList.union(duo);
    }


    /** inverts an IntervalList and returns one that has exactly all the bases in the dictionary that the original one does not.
     *
     * @param list an IntervalList
     * @return an IntervalList that is complementary to list
     */
    public static IntervalList invert(final IntervalList list) {
        final IntervalList inverse = new IntervalList(list.header.clone());

        final ListMap<Integer,Interval> map = new ListMap<Integer,Interval>();

        //add all the intervals (uniqued and therefore also sorted) to a ListMap from sequenceIndex to a list of Intervals
        for(final Interval i : list.uniqued().getIntervals()){
            map.add(list.getHeader().getSequenceIndex(i.getSequence()),i);
        }

        int intervals = 0; // a counter to supply newly-created intervals with a name


        //iterate over the contigs in the dictionary
        for (final SAMSequenceRecord samSequenceRecord : list.getHeader().getSequenceDictionary().getSequences()) {
            final Integer sequenceIndex = samSequenceRecord.getSequenceIndex();
            final String sequenceName   = samSequenceRecord.getSequenceName();
            final int sequenceLength    = samSequenceRecord.getSequenceLength();

            Integer lastCoveredPosition = 0; //start at beginning of sequence
            //iterate over list of intervals that are in sequence
            if (map.containsKey(sequenceIndex)) // if there are intervals in the ListMap on this contig, iterate over them (in order)
                for (final Interval i : map.get(sequenceIndex)) {
                    if (i.getStart() > lastCoveredPosition + 1) //if there's space between the last interval and the current one, add an interval between them
                        inverse.add(new Interval(sequenceName, lastCoveredPosition + 1, i.getStart() - 1, false, "interval-" + (++intervals)));
                    lastCoveredPosition = i.getEnd(); //update the last covered position
                }
            //finally, if there's room between the last covered position and the end of the sequence, add an interval
            if (sequenceLength > lastCoveredPosition) //if there's space between the last interval and the next
                // one, add an interval. This also covers the case that there are no intervals in the ListMap for a contig.
                inverse.add(new Interval(sequenceName, lastCoveredPosition + 1, sequenceLength, false, "interval-" + (++intervals)));
        }

        return inverse;
    }


    /**
     * A utility function for subtracting a collection of IntervalLists from another. Resulting loci are those that are in the first collection
     * but not the second.
     *
     * @param listsToSubtractFrom the collection of IntervalList from which to subtract intervals
     * @param listsToSubtract the collection of intervals to subtract
     * @return an IntervalLists comprising all loci that are in first collection but not second.
     */
    public static IntervalList subtract(final Collection<IntervalList> listsToSubtractFrom, final Collection<IntervalList> listsToSubtract) {
        return intersection(
                union(listsToSubtractFrom),
                invert(union(listsToSubtract)));
    }


    /**
     * A utility function for finding the difference between two IntervalLists.
     *
     * @param lists1 the first collection of IntervalLists
     * @param lists2 the second collection of IntervalLists
     * @return the difference between the two intervals, i.e. the loci that are only in one IntervalList but not both
     */
    public static IntervalList difference(final Collection<IntervalList> lists1, final Collection<IntervalList> lists2) {
        return union(
                subtract(lists1, lists2),
                subtract(lists2, lists1));
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
            if (lhs.getName() == null) {
                if (rhs.getName() == null) return 0;
                else return -1;
            } else if (rhs.getName() == null) {
                return 1;
            }
            else {
                return lhs.getName().compareTo(rhs.getName());
            }
        }

        return retval;
    }
}