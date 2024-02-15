package picard.sam.util;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.sam.SamComparisonMetric;

import java.io.File;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

/**
 * Compare two SAM/BAM files.  Compares headers, and if headers are compatible enough, compares SAMRecords,
 * looking at alignment and duplicate marking info.  Can perform either a naive comparison for which each alignment must be identical, or a more sophisticated check of "equivalence", where mapping quality
 * 0 reads are allowed to have different alignments, and duplicate marks are allowed to differ to account for ambiguities in selecting the representative read of a duplicate set.  Results of comparison are
 * summarised in an output metrics file.
 */
public final class SamComparison {
    private final SamReader leftReader;
    private final SamReader rightReader;

    private boolean sequenceDictionariesDiffer;

    private final SamComparisonMetric comparisonMetric = new SamComparisonMetric();

    private final SAMComparisonArgumentCollection samComparisonArgumentCollection;

    // Ideally, this should be a Histogram<Pair<Integer, Integer>>, however, the MetricsFile class cannot read
    // arbitrary types, therefore, it must be converted to a String, which may be slower
    private final Histogram<String> mappingQualityHistogram = new Histogram<>();

    private SortingCollection<SAMRecord> markDuplicatesCheckLeft;
    private SortingCollection<SAMRecord> markDuplicatesCheckRight;

    public enum AlignmentComparison {
        UNMAPPED_BOTH, UNMAPPED_LEFT, UNMAPPED_RIGHT, MAPPINGS_DIFFER, MAPPINGS_MATCH
    }

    /**
     * Note: the caller must make sure the SamReaders are closed properly.
     */
    public SamComparison(final SamReader leftReader, final SamReader rightReader) {
        this(leftReader, rightReader, null, null, new SAMComparisonArgumentCollection());
    }

    public SamComparison(final SamReader leftReader, final SamReader rightReader, final String leftName, final String rightName, final SAMComparisonArgumentCollection samComparisonArgumentCollection) {
        this.leftReader = leftReader;
        this.rightReader = rightReader;
        this.samComparisonArgumentCollection = samComparisonArgumentCollection;
        comparisonMetric.LEFT_FILE = leftName;
        comparisonMetric.RIGHT_FILE = rightName;
        if (samComparisonArgumentCollection.LENIENT_DUP) {
            setupLenientDuplicateChecking();
        }
        comparisonMetric.ARE_EQUAL = compareHeaders();
        comparisonMetric.ARE_EQUAL &= compareAlignmentsAndCatalogDuplicateMarkingDifferences();
        if (samComparisonArgumentCollection.LENIENT_DUP) {
            countLenientDuplicateMarkingDifferences();
        }
        comparisonMetric.ARE_EQUAL &= comparisonMetric.DUPLICATE_MARKINGS_DIFFER == 0;
    }

    public void writeReport(final File output) {
        writeReport(output, Collections.EMPTY_LIST);
    }

    public void writeReport(final File output, final List<Header> headers) {
        final MetricsFile<SamComparisonMetric, String> comparisonMetricFile = new MetricsFile<>();

        headers.forEach(comparisonMetricFile::addHeader);
        comparisonMetricFile.addAllMetrics(Collections.singletonList(comparisonMetric));
        if (samComparisonArgumentCollection.COMPARE_MQ) {
            comparisonMetricFile.addHistogram(mappingQualityHistogram);
        }
        comparisonMetricFile.write(output);
    }

    private void setupLenientDuplicateChecking() {
        /* Setup for duplicate marking checks.  Recs which disagree on duplicate marking will be added to markDuplicatesCheckLeft/Right and then fed to a DuplicateSetIterator to check if
         * differences are due only to choice of representative read within duplicate set.
         */
        final SAMFileHeader leftHeader = leftReader.getFileHeader();
        final SAMFileHeader rightHeader = rightReader.getFileHeader();
        final SAMRecordDuplicateComparator leftDupComparator = new SAMRecordDuplicateComparator(Collections.singletonList(leftHeader));
        final SAMRecordDuplicateComparator rightDupComparator = new SAMRecordDuplicateComparator(Collections.singletonList(rightHeader));
        markDuplicatesCheckLeft = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(leftHeader),
                leftDupComparator, SAMFileWriterImpl.getDefaultMaxRecordsInRam());
        markDuplicatesCheckRight = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(rightHeader),
                rightDupComparator, SAMFileWriterImpl.getDefaultMaxRecordsInRam());
    }

    private void countLenientDuplicateMarkingDifferences() {
        /* Count duplicate marking differences allowing swaps within duplicate sets.  Should only be used when lenientDup == true.
           Reads which have different duplicate markings will have been saved in markDuplicatesCheckLeft and markDuplicatesCheckRight
           in catalogDuplicateDifferences.
         */

        if (!samComparisonArgumentCollection.LENIENT_DUP) {
            throw new PicardException("Should only use countLenientDuplicateMarkingDifferences when in lenient duplicate marking mode.");
        }

        final DuplicateSetIterator duplicateSetLeftIterator = new DuplicateSetIterator(markDuplicatesCheckLeft.iterator(), leftReader.getFileHeader(), true);
        final DuplicateSetIterator duplicateSetRightIterator = new DuplicateSetIterator(markDuplicatesCheckRight.iterator(), rightReader.getFileHeader(), true);
        /* Iterate through duplicateSetIterators.  This will iterate through all reads with differing duplicate markings in left and right files.  Passing through the first
         * iterator (i.e. left file), for each duplicate set we create a set of duplicate reads that the representative read could "swap" with.  On the pass through the second iterator, for each
         * duplicate read, we look for the reads it is allowed to swap with in the duplicate set from the first iterator.  If the allowed swap can also be performed in this duplicate set,
         * we do not count the read or the read it could be swapped with as having mismatched duplicates.  Note however, that once a swap is allowed, the reads swapped cannot be used for any
         * further swaps.
         */
        // count reads which are not saved
        final HashMap<String, Set<String>> swapTargetsForLeftReps = new HashMap<>();
        for (final DuplicateSet duplicateSet : new IterableAdapter<>(duplicateSetLeftIterator)) {
           // do not want to redo duplicate marking, want to keep marks assigned in input files
           final List<SAMRecord> allRecords = duplicateSet.getRecords(false);
           if (allRecords.size() > 1) {
               // for this set, get a list of non-rep reads that can potentially swap with the rep read
               final Set<String> nonRepReads = allRecords.stream()
                       .filter(SAMRecord::getDuplicateReadFlag)
                       .map(SAMRecord::getReadName)
                       .collect(Collectors.toSet());
               allRecords.stream().filter(r -> !r.getDuplicateReadFlag())
                       .forEach(leftRep -> swapTargetsForLeftReps.put(leftRep.getReadName(), nonRepReads));
           }
        }
        // set of fragments which have been "matched", i.e. will not be counted as having mismatched duplicate marks.
        final HashSet<String> matchedSwaps = new HashSet<>();
        for (final DuplicateSet duplicateSet : new IterableAdapter<>(duplicateSetRightIterator)) {
           // do not want to redo duplicate marking, want to keep marking assigned in files
           final List<SAMRecord> records = duplicateSet.getRecords(false);
           if (records.size() > 1) {
               // create a list of non-rep reads for this set that were also rep-reads in the left file
               final List<String> nonRepSwapCandidates = records.stream()
                       .filter(SAMRecord::getDuplicateReadFlag)
                       .map(SAMRecord::getReadName)
                       .filter(swapTargetsForLeftReps::containsKey)
                       .collect(Collectors.toList());
               // now try to match up the rep reads in this set with a non-rep swap candidate
               records.stream().filter(r -> !r.getDuplicateReadFlag()).map(SAMRecord::getReadName)
                       .filter(repRead -> !matchedSwaps.contains(repRead))
                       .forEach(unMatchedRep ->
                           nonRepSwapCandidates.stream()
                                   .filter(nonRep -> !matchedSwaps.contains(nonRep) && swapTargetsForLeftReps.get(nonRep).contains(unMatchedRep))
                                   .findFirst()
                                   .ifPresent(matchedRep -> matchedSwaps.addAll(Arrays.asList(matchedRep, unMatchedRep)))
                   );
           }
           // count reads which differ, but for which no swap match was found
           comparisonMetric.DUPLICATE_MARKINGS_DIFFER += records.stream().filter(n -> !matchedSwaps.contains(n.getReadName())).count();
        }

    }

    private boolean compareAlignmentsAndCatalogDuplicateMarkingDifferences() {
        /* Compares alignments of primary sam records.  Eventual comparison of sam records found in both files is performed in tallyAlignments.  tallyAlignments
         * includes duplicate marking cataloging, which counts all duplicate marking differences when in strict mode, or adds reads which differ in duplicate marking
         * to sam record iterators to be counted while allowing for swaps later in countLenientDuplicateMarkingDifferences.
         */
        if (!compareValues(leftReader.getFileHeader().getSortOrder(), rightReader.getFileHeader().getSortOrder(),
                "Sort Order")) {
            System.out.println("Cannot compare alignments if sort orders differ.");
            return false;
        }
        switch (leftReader.getFileHeader().getSortOrder()) {
            case coordinate:
                if (sequenceDictionariesDiffer) {
                    System.out.println("Cannot compare coordinate-sorted SAM files because sequence dictionaries differ.");
                    return false;
                }
                return compareCoordinateSortedAlignments();
            case queryname:
                return compareQueryNameSortedAlignments();
            case duplicate:
            case unsorted:
                return compareUnsortedAlignments();
            default:
                // unreachable
                throw new PicardException(String.format("Unrecognized sort order (%s) found.", leftReader.getFileHeader().getSortOrder()));
        }
    }

    private boolean compareCoordinateSortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator itLeft =
                new SecondaryOrSupplementarySkippingIterator(leftReader.iterator());
        final SecondaryOrSupplementarySkippingIterator itRight =
                new SecondaryOrSupplementarySkippingIterator(rightReader.iterator());

        // Save any reads which haven't been matched during in-order scan.
        final Map<PrimaryAlignmentKey, SAMRecord> leftUnmatched = new LinkedHashMap<>();
        final Map<PrimaryAlignmentKey, SAMRecord> rightUnmatched = new LinkedHashMap<>();

        while (itLeft.hasCurrent()) {
            if (!itRight.hasCurrent()) {
                // Exhausted right side.  See if any of the remaining left reads match
                // any of the saved right reads.
                for (; itLeft.hasCurrent(); itLeft.advance()) {
                    final SAMRecord left = itLeft.getCurrent();
                    final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(left);
                    final SAMRecord right = rightUnmatched.remove(leftKey);
                    if (right == null) {
                        ++comparisonMetric.MISSING_RIGHT;
                    } else {
                        tallyAlignmentRecords(left, right);
                    }
                }
                break;
            }
            // Don't assume stability of order beyond the coordinate.  Therefore grab all the
            // reads from the left that has the same coordinate.
            final SAMRecord left = itLeft.getCurrent();
            final Map<PrimaryAlignmentKey, SAMRecord> leftCurrentCoordinate = new LinkedHashMap<>();
            final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(left);
            leftCurrentCoordinate.put(leftKey, left);
            while (itLeft.advance()) {
                final SAMRecord nextLeft = itLeft.getCurrent();
                if (compareAlignmentCoordinates(left, nextLeft) == 0) {
                    final PrimaryAlignmentKey nextLeftKey = new PrimaryAlignmentKey(nextLeft);
                    leftCurrentCoordinate.put(nextLeftKey, nextLeft);
                } else {
                    break;
                }
            }
            // Advance the right iterator until it is >= the left reads that have just been grabbed
            while (itRight.hasCurrent() && compareAlignmentCoordinates(left, itRight.getCurrent()) > 0) {
                final SAMRecord right = itRight.getCurrent();
                final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(right);
                rightUnmatched.put(rightKey, right);
                itRight.advance();
            }
            // For each right read that has the same coordinate as the current left reads,
            // see if there is a matching left read.  If so, process and discard.  If not,
            // save the right read for later.
            for (; itRight.hasCurrent() && compareAlignmentCoordinates(left, itRight.getCurrent()) == 0; itRight.advance()) {
                final SAMRecord right = itRight.getCurrent();
                final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(right);
                final SAMRecord matchingLeft = leftCurrentCoordinate.remove(rightKey);
                if (matchingLeft != null) {
                    tallyAlignmentRecords(matchingLeft, right);
                } else {
                    rightUnmatched.put(rightKey, right);
                }
            }

            // Anything left in leftCurrentCoordinate has not been matched
            for (final SAMRecord samRecord : leftCurrentCoordinate.values()) {
                final PrimaryAlignmentKey recordKey = new PrimaryAlignmentKey(samRecord);
                leftUnmatched.put(recordKey, samRecord);
            }
        }
        // The left iterator has been exhausted.  See if any of the remaining right reads
        // match any of the saved left reads.
        consumeUnmatchedRights(itRight, leftUnmatched);

        // Look up reads that were unmatched from left, and see if they are in rightUnmatched.
        // If found, remove from rightUnmatched and tally.
        for (final Map.Entry<PrimaryAlignmentKey, SAMRecord> leftEntry : leftUnmatched.entrySet()) {
            final PrimaryAlignmentKey leftKey = leftEntry.getKey();
            final SAMRecord left = leftEntry.getValue();
            final SAMRecord right = rightUnmatched.remove(leftKey);
            if (right == null) {
                ++comparisonMetric.MISSING_RIGHT;
                continue;
            }
            tallyAlignmentRecords(left, right);
        }

        // Any elements remaining in rightUnmatched are guaranteed not to be in leftUnmatched.
        comparisonMetric.MISSING_LEFT += rightUnmatched.size();

        return comparisonMetric.allVisitedAlignmentsEqual();
    }

    private int compareAlignmentCoordinates(final SAMRecord left, final SAMRecord right) {
        final String leftReferenceName = left.getReferenceName();
        final String rightReferenceName = right.getReferenceName();
        if (leftReferenceName == null && rightReferenceName == null) {
            return 0;
        } else if (leftReferenceName == null) {
            return 1;
        } else if (rightReferenceName == null) {
            return -1;
        }
        final int leftReferenceIndex = leftReader.getFileHeader().getSequenceIndex(leftReferenceName);
        final int rightReferenceIndex = rightReader.getFileHeader().getSequenceIndex(rightReferenceName);

        if (leftReferenceIndex != rightReferenceIndex) {
            return leftReferenceIndex - rightReferenceIndex;
        }
        return left.getAlignmentStart() - right.getAlignmentStart();
    }

    private boolean compareQueryNameSortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(leftReader.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(rightReader.iterator());

        while (it1.hasCurrent()) {
            if (!it2.hasCurrent()) {
                comparisonMetric.MISSING_RIGHT += countRemaining(it1);
            }
            final PrimaryAlignmentKey leftKey = new PrimaryAlignmentKey(it1.getCurrent());
            final PrimaryAlignmentKey rightKey = new PrimaryAlignmentKey(it2.getCurrent());
            final int cmp = leftKey.compareTo(rightKey);
            if (cmp < 0) {
                ++comparisonMetric.MISSING_RIGHT;
                it1.advance();
            } else if (cmp > 0) {
                ++comparisonMetric.MISSING_LEFT;
                it2.advance();
            } else {
                tallyAlignmentRecords(it1.getCurrent(), it2.getCurrent());
                it1.advance();
                it2.advance();
            }
        }
        if (it2.hasCurrent()) {
            comparisonMetric.MISSING_LEFT += countRemaining(it2);
        }
        return comparisonMetric.allVisitedAlignmentsEqual();
    }

    /**
     * For unsorted alignments, assume nothing about the order. Determine which records to compare solely on the
     * basis of their PrimaryAlignmentKey.
     */
    private boolean compareUnsortedAlignments() {
        final SecondaryOrSupplementarySkippingIterator it1 = new SecondaryOrSupplementarySkippingIterator(leftReader.iterator());
        final SecondaryOrSupplementarySkippingIterator it2 = new SecondaryOrSupplementarySkippingIterator(rightReader.iterator());

        final Map<PrimaryAlignmentKey, SAMRecord> leftUnmatched = new LinkedHashMap<>();

        // Consume all of the lefts by adding them to the leftUnmatched map, then consume all of the
        // rights by finding and removing their left mate from the "leftUnmatched" map. Anything remaining
        // after that is unmatched.
        consumeAll(it1, (alignmentRecord, primaryKey) -> leftUnmatched.put(primaryKey, alignmentRecord));
        consumeUnmatchedRights(it2, leftUnmatched);

        comparisonMetric.MISSING_RIGHT += leftUnmatched.size();

        return comparisonMetric.allVisitedAlignmentsEqual();
    }

    /**
     * Consume every record in the "right" iterator by either finding its matched "left", or
     * failing that, incrementing the number of missing "lefts".
     */
    private void consumeUnmatchedRights(
            final SecondaryOrSupplementarySkippingIterator rightIt,
            final Map<PrimaryAlignmentKey, SAMRecord> leftUnmatched) {
        consumeAll(rightIt,
                (alignmentRecord, primaryKey) -> {
                    final SAMRecord left = leftUnmatched.remove(primaryKey);
                    if (left != null) {
                        tallyAlignmentRecords(left, alignmentRecord);
                    } else {
                        ++comparisonMetric.MISSING_LEFT;
                    }
                });
    }

    /**
     * Consume every record in the iterator, passing it and the corresponding PrimaryAlignmentKey.
     * to the provided handler.
     */
    private void consumeAll(
            final SecondaryOrSupplementarySkippingIterator it,
            final BiConsumer<SAMRecord, PrimaryAlignmentKey> recordKeyHandler) {
        for (; it.hasCurrent(); it.advance()) {
            final SAMRecord record = it.getCurrent();
            recordKeyHandler.accept(record, new PrimaryAlignmentKey(record));
        }
    }

    private int countRemaining(final SecondaryOrSupplementarySkippingIterator it) {
        int i;
        for (i = 0; it.hasCurrent(); ++i) {
            it.advance();
        }
        return i;
    }

    private AlignmentComparison compareAlignmentRecords(final SAMRecord s1, final SAMRecord s2) {
        if (s1.getReadUnmappedFlag() && s2.getReadUnmappedFlag()) {
            return AlignmentComparison.UNMAPPED_BOTH;
        } else if (s1.getReadUnmappedFlag()) {
            return AlignmentComparison.UNMAPPED_LEFT;
        } else if (s2.getReadUnmappedFlag()) {
            return AlignmentComparison.UNMAPPED_RIGHT;
        } else if (alignmentsMatch(s1, s2)) {
            return AlignmentComparison.MAPPINGS_MATCH;
        } else {
            return AlignmentComparison.MAPPINGS_DIFFER;
        }
    }

    private boolean alignmentsMatch(final SAMRecord s1, final SAMRecord s2) {
        return ((s1.getReferenceName().equals(s2.getReferenceName()) &&
                s1.getAlignmentStart() == s2.getAlignmentStart() &&
                s1.getReadNegativeStrandFlag() == s1.getReadNegativeStrandFlag()) ||
                (samComparisonArgumentCollection.LENIENT_LOW_MQ_ALIGNMENT && s1.getMappingQuality() <= samComparisonArgumentCollection.LOW_MQ_THRESHOLD &&
                        s2.getMappingQuality() <= samComparisonArgumentCollection.LOW_MQ_THRESHOLD) ||
                (samComparisonArgumentCollection.LENIENT_UNKNOWN_MQ_ALIGNMENT && s1.getMappingQuality() == SAMRecord.UNKNOWN_MAPPING_QUALITY && s2.getMappingQuality() == SAMRecord.UNKNOWN_MAPPING_QUALITY));
    }

    private void compareAndUpdateMappingQualityConcordance(final SAMRecord s1, final SAMRecord s2) {
        mappingQualityHistogram.increment(String.format("%d,%d", s1.getMappingQuality(), s2.getMappingQuality()));
    }

    /**
     * Compare the mapping information for two SAMRecords.  Makes comparison of alignments, and also catalogs duplicate marking differences.
     */
    private void tallyAlignmentRecords(final SAMRecord s1, final SAMRecord s2) {
        if (!s1.getReadName().equals(s2.getReadName())) {
            throw new PicardException("Read names do not match: " + s1.getReadName() + " : " + s2.getReadName());
        }
        catalogDuplicateDifferences(s1, s2);
        final AlignmentComparison comp = compareAlignmentRecords(s1, s2);
        comparisonMetric.updateMetric(comp);

        if (samComparisonArgumentCollection.COMPARE_MQ) {
            compareAndUpdateMappingQualityConcordance(s1, s2);
        }
    }

    private void catalogDuplicateDifferences(final SAMRecord s1, final SAMRecord s2) {
        // if strict, reads with differing duplicate marking are counted by DUPLICATE_MARKINGS_DIFFER.
        // if lenient, reads with differing duplicate marking are added to markDuplicatesCheckLeft/Right
        // to later be counted while allowing for swaps withing duplicate sets by updateLenientDuplicateMarkingDifferences
        if (s1.getDuplicateReadFlag() != s2.getDuplicateReadFlag()) {
            if (samComparisonArgumentCollection.LENIENT_DUP) {
                markDuplicatesCheckLeft.add(s1);
                markDuplicatesCheckRight.add(s2);
            } else {
                comparisonMetric.DUPLICATE_MARKINGS_DIFFER++;
            }
        }
    }

    private boolean compareHeaders() {
        final SAMFileHeader h1 = leftReader.getFileHeader();
        final SAMFileHeader h2 = rightReader.getFileHeader();
        boolean ret = compareValues(h1.getVersion(), h2.getVersion(), "File format version");
        ret = compareValues(h1.getCreator(), h2.getCreator(), "File creator") && ret;
        ret = compareValues(h1.getAttribute("SO"), h2.getAttribute("SO"), "Sort order") && ret;
        if (!compareSequenceDictionaries(h1, h2)) {
            ret = false;
            sequenceDictionariesDiffer = true;
        }
        ret = compareReadGroups(h1, h2) && ret;
        if (!samComparisonArgumentCollection.LENIENT_HEADER) {
            ret = compareProgramRecords(h1, h2) && ret;
        }
        return ret;
    }

    private boolean compareProgramRecords(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMProgramRecord> l1 = h1.getProgramRecords();
        final List<SAMProgramRecord> l2 = h2.getProgramRecords();
        if (!compareValues(l1.size(), l2.size(), "Number of program records")) {
            return false;
        }
        boolean ret = true;
        for (int i = 0; i < l1.size(); ++i) {
            ret = compareProgramRecord(l1.get(i), l2.get(i)) && ret;
        }
        return ret;
    }

    private static boolean compareProgramRecord(final SAMProgramRecord programRecord1, final SAMProgramRecord programRecord2) {
        if (programRecord1 == null && programRecord2 == null) {
            return true;
        }
        if (programRecord1 == null) {
            reportDifference("null", programRecord2.getProgramGroupId(), "Program Record");
            return false;
        }
        if (programRecord2 == null) {
            reportDifference(programRecord1.getProgramGroupId(), "null", "Program Record");
            return false;
        }
        boolean ret = compareValues(programRecord1.getProgramGroupId(), programRecord2.getProgramGroupId(),
                "Program Name");
        final String[] attributes = {"VN", "CL"};
        for (final String attribute : attributes) {
            ret = compareValues(programRecord1.getAttribute(attribute), programRecord2.getAttribute(attribute),
                    attribute + " Program Record attribute") && ret;
        }
        return ret;
    }

    private static boolean compareReadGroups(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMReadGroupRecord> l1 = h1.getReadGroups();
        final List<SAMReadGroupRecord> l2 = h2.getReadGroups();
        if (!compareValues(l1.size(), l2.size(), "Number of read groups")) {
            return false;
        }
        boolean ret = true;
        for (int i = 0; i < l1.size(); ++i) {
            ret = compareReadGroup(l1.get(i), l2.get(i)) && ret;
        }
        return ret;
    }

    private static boolean compareReadGroup(final SAMReadGroupRecord samReadGroupRecord1, final SAMReadGroupRecord samReadGroupRecord2) {
        boolean ret = compareValues(samReadGroupRecord1.getReadGroupId(), samReadGroupRecord2.getReadGroupId(),
                "Read Group ID");
        ret = compareValues(samReadGroupRecord1.getSample(), samReadGroupRecord2.getSample(),
                "Sample for read group " + samReadGroupRecord1.getReadGroupId()) && ret;
        ret = compareValues(samReadGroupRecord1.getLibrary(), samReadGroupRecord2.getLibrary(),
                "Library for read group " + samReadGroupRecord1.getReadGroupId()) && ret;
        final String[] attributes = {"DS", "PU", "PI", "CN", "DT", "PL"};
        for (final String attribute : attributes) {
            ret = compareValues(samReadGroupRecord1.getAttribute(attribute), samReadGroupRecord2.getAttribute(attribute),
                    attribute + " for read group " + samReadGroupRecord1.getReadGroupId()) && ret;
        }
        return ret;
    }

    private boolean compareSequenceDictionaries(final SAMFileHeader h1, final SAMFileHeader h2) {
        final List<SAMSequenceRecord> s1 = h1.getSequenceDictionary().getSequences();
        final List<SAMSequenceRecord> s2 = h2.getSequenceDictionary().getSequences();
        if (s1.size() != s2.size()) {
            reportDifference(s1.size(), s2.size(), "Length of sequence dictionaries");
            return false;
        }
        boolean ret = true;
        for (int i = 0; i < s1.size(); ++i) {
            ret = compareSequenceRecord(s1.get(i), s2.get(i), i + 1) && ret;
        }
        return ret;
    }

    private boolean compareSequenceRecord(final SAMSequenceRecord sequenceRecord1, final SAMSequenceRecord sequenceRecord2, final int which) {
        if (!sequenceRecord1.getSequenceName().equals(sequenceRecord2.getSequenceName())) {
            reportDifference(sequenceRecord1.getSequenceName(), sequenceRecord2.getSequenceName(),
                    "Name of sequence record " + which);
            return false;
        }
        boolean ret = compareValues(sequenceRecord1.getSequenceLength(), sequenceRecord2.getSequenceLength(), "Length of sequence " +
                sequenceRecord1.getSequenceName());
        if (!samComparisonArgumentCollection.LENIENT_HEADER) {
            ret = compareValues(sequenceRecord1.getSpecies(), sequenceRecord2.getSpecies(), "Species of sequence " +
                    sequenceRecord1.getSequenceName()) && ret;
            ret = compareValues(sequenceRecord1.getAssembly(), sequenceRecord2.getAssembly(), "Assembly of sequence " +
                    sequenceRecord1.getSequenceName()) && ret;
            ret = compareValues(sequenceRecord1.getAttribute("M5"), sequenceRecord2.getAttribute("M5"), "MD5 of sequence " +
                    sequenceRecord1.getSequenceName()) && ret;
            ret = compareValues(sequenceRecord1.getAttribute("UR"), sequenceRecord2.getAttribute("UR"), "URI of sequence " +
                    sequenceRecord1.getSequenceName()) && ret;
        }
        return ret;
    }

    private static <T> boolean compareValues(final T v1, final T v2, final String label) {
        boolean eq = Objects.equals(v1, v2);
        if (eq) {
            return true;
        } else {
            reportDifference(v1, v2, label);
            return false;
        }
    }

    private static void reportDifference(final String s1, final String s2, final String label) {
        System.out.println(label + " differs.");
        System.out.println("File 1: " + s1);
        System.out.println("File 2: " + s2);
    }

    private static void reportDifference(Object o1, Object o2, final String label) {
        if (o1 == null) {
            o1 = "null";
        }
        if (o2 == null) {
            o2 = "null";
        }
        reportDifference(o1.toString(), o2.toString(), label);
    }

    public int getMappingsMatch() {
        return comparisonMetric.MAPPINGS_MATCH;
    }

    public int getUnmappedBoth() {
        return comparisonMetric.UNMAPPED_BOTH;
    }

    public int getUnmappedLeft() {
        return comparisonMetric.UNMAPPED_LEFT;
    }

    public int getUnmappedRight() {
        return comparisonMetric.UNMAPPED_RIGHT;
    }

    public int getMappingsDiffer() {
        return comparisonMetric.MAPPINGS_DIFFER;
    }

    public int getMissingLeft() {
        return comparisonMetric.MISSING_LEFT;
    }

    public int getMissingRight() {
        return comparisonMetric.MISSING_RIGHT;
    }

    public int getDuplicateMarkingsDiffer() {
        return comparisonMetric.DUPLICATE_MARKINGS_DIFFER;
    }

    public boolean areEqual() {
        return comparisonMetric.ARE_EQUAL;
    }
}
