/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollectionUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.SamOrBam;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.ArgsPreparer;
import picard.sam.markduplicates.util.DiskBasedReadEndsForMarkDuplicatesMap;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.markduplicates.util.ReadEnds;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesWithBarcodes;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.RepresentativeReadIndexer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Parallel implementation of MarkDuplicates.
 */
@CommandLineProgramProperties(
        summary = ParallelMarkDuplicates.USAGE_SUMMARY + ParallelMarkDuplicates.USAGE_DETAILS,
        oneLineSummary = ParallelMarkDuplicates.USAGE_SUMMARY,
        programGroup = SamOrBam.class)
@DocumentedFeature
public class ParallelMarkDuplicates extends MarkDuplicates {
    static final String USAGE_SUMMARY = "Parallel implementation of MarkDuplicates.  ";
    static final String USAGE_DETAILS = "<p>This tool fully inherits original MarkDuplicates interface except the fact " +
            "that ParallelMarkDuplicates could work only with a single indexed BAM file. So, only coordinate sorted BAM " +
            "is valid.</p>" +
            "" +
            "<p>Parallel version will split input file by chromosomes to process them in parallel. ParallelMarkDuplicates " +
            "will try to split original file by PARALLELISM * SPLIT_FACTOR to scale up to PARALLELISM threads more smoothly.</p>" +
            "" +
            "<p>After step of parallel reading and reads sorting there is a sequential part of resolving cross-chromosome " +
            "paired reads between parts processed. On the next step for each split ParallelMarkDuplicates will produce " +
            "temp output file with marked duplicates. All those outputs will be merged to OUTPUT at the end and be deleted. " +
            "So, parallel version will use twice more disc space while working then original.</p>" +
            "<p>To get more details see original MarkDuplicates documentation." +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar ParallelMarkDuplicates \\<br />" +
            "      TH=8 \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=marked_duplicates.bam \\<br />" +
            "      M=marked_dup_metrics.txt" +
            "</pre>" +
            "" +
            "Please see " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics'>MarkDuplicates</a> " +
            "for detailed explanations of the output metrics." +
            "<hr />";

    private final Log log = Log.getInstance(ParallelMarkDuplicates.class);

    @Argument(shortName = "TH", doc = "Number of threads, which the program will try to utilize")
    public int PARALLELISM = 2;

    @Argument(doc = "ParallelMarkDuplicates will try to split original file by PARALLELISM * SPLIT_FACTOR to scale " +
            "up to PARALLELISM threads more smoothly. In general this value should not be changed.")
    public float SPLIT_FACTOR = 6;

    @Argument(doc = "Coefficient of SortingCollections size in each MarkDuplicatesInstance working in parallel." +
            "In general this value should not be changed.")
    private float SORTING_COLLECTION_SIZE_FACTOR = 0.5f;

    private ExecutorService threadPool;

    private List<MarkDuplicatesInstance> mds = new ArrayList<>();
    private NavigableMap<Integer, MarkDuplicatesInstance> mdsMap = new TreeMap<>();

    public static void main(String[] args) {
        new ParallelMarkDuplicates().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        if (!assertInputs()) return 1;

        threadPool = new ForkJoinPool(PARALLELISM);
        useBarcodes = (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG);

        log.info("Waiting building read ends...");
        buildSortedReadEndsListsParallel().join();

        log.info("Merging pairs...");
        mergeReadEnds();

        log.info("Waiting for writing...");
        writeBamFilesParallel().join();

        log.info("Merging output files...");
        mergeOutputs();

        log.info("Merging and writing metrics...");
        mergeMetrics();

        finalizeAndWriteMetrics(libraryIdGenerator);

        //thread pool actually should be empty
        try {
            threadPool.shutdown();
            threadPool.awaitTermination(1, TimeUnit.MINUTES);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }

        return 0;
    }

    private boolean assertInputs() {
        if (INPUT.size() > 1) {
            log.error("For parallel implementation only single coordinate sorted INPUT is valid.");
            return false;
        }

        if (PARALLELISM < 1) {
            log.error("Number of thread must be a positive number.");
            return false;
        }

        if (PARALLELISM == 1) {
            log.warn("You are running ParallelMarkDuplicates with only one thread allowed. In that case better run " +
                    "sequential version (MarkDuplicates) to save disk space usage during processing.");
        }

        IOUtil.assertInputsAreValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        return true;
    }

    private CompletableFuture<Void> buildSortedReadEndsListsParallel() {
        try (SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(INPUT.get(0)))) {
            SAMFileHeader header = samReader.getFileHeader();

            this.libraryIdGenerator = new LibraryIdGenerator(header);

            final long recordsInFile = header.getSequenceDictionary().getSequences().stream()
                    .mapToLong(seq -> {
                        long seqRecordCount = samReader.indexing().getIndex().getMetaData(seq.getSequenceIndex()).getAlignedRecordCount();
                        seqRecordCount += samReader.indexing().getIndex().getMetaData(seq.getSequenceIndex()).getUnalignedRecordCount();

                        return seqRecordCount;
                    })
                    .sum();
            final long maxRecordsInPart = (long) (recordsInFile / (PARALLELISM * SPLIT_FACTOR));

            List<CompletableFuture<Void>> buildSortedReadEndListsFutures = new ArrayList<>();

            int prevPartsRecordsCount = 0;
            int currentRecordsCount = 0;
            List<SAMSequenceRecord> seqBatch = new ArrayList<>();

            for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
                int seqRecordCount = samReader.indexing().getIndex().getMetaData(seq.getSequenceIndex()).getAlignedRecordCount();
                seqRecordCount += samReader.indexing().getIndex().getMetaData(seq.getSequenceIndex()).getUnalignedRecordCount();

                if (seqRecordCount == 0) {
                    log.debug("Skip seq ", seq.getSequenceName(), " because there is no reads");
                    continue;
                }

                seqBatch.add(seq);
                currentRecordsCount += seqRecordCount;

                if (currentRecordsCount - prevPartsRecordsCount >= maxRecordsInPart) {
                    buildSortedReadEndListsFutures.add(
                            launchMDInstance(new MarkDuplicatesInstance(seqBatch, prevPartsRecordsCount)));

                    seqBatch = new ArrayList<>();
                    prevPartsRecordsCount = currentRecordsCount;
                }
            }

            //last one if presented
            if (!seqBatch.isEmpty()) {
                buildSortedReadEndListsFutures.add(
                        launchMDInstance(new MarkDuplicatesInstance(seqBatch, prevPartsRecordsCount)));
            }

            // unmapped part
            buildSortedReadEndListsFutures.add(
                    launchMDInstance(new MarkDuplicatesInstance()));

            return CompletableFuture.allOf(
                    buildSortedReadEndListsFutures.toArray(new CompletableFuture[buildSortedReadEndListsFutures.size()])
            );
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private CompletableFuture<Void> launchMDInstance(MarkDuplicatesInstance md) {
        md.instanceMain(prepareArgsForMDInstances(getOriginalArgs(), INPUT.get(0)));
        mds.add(md);
        mdsMap.put(md.leftmostSeqIndex, md);
        log.debug("Queueing MarkDuplicates instance on region ", md.fileNamePrefix(), "...");

        return CompletableFuture.runAsync(md::buildSortedReadEndLists, threadPool);
    }

    private String[] prepareArgsForMDInstances(String[] argv, String input) {
        return new ArgsPreparer(argv)
                /* those only for parallel version */
                .exclude("PARALLELISM").exclude("TH")
                .exclude("SORTING_COLLECTION_SIZE_FACTOR")
                .exclude("SPLIT_FACTOR")
                /* to be sure QUIET=true */
                .exclude("QUIET")
                .add("QUIET", "true")
                /* input should be changed */
                .exclude("INPUT").exclude("I")
                .add("INPUT", input)
                .toArray();
    }

    /**
     * Algorithm of pair ends combining is optimized to minimize count of loading in RAM and Unloading to disk of
     * sequence maps in CoordinateSortedPairInfoMap (DiskBasedReadEndsForMarkDuplicatesMap).
     * Although there are 3 nested loops, merging has linear complexity.
     */
    private void mergeReadEnds() {
        int endsCount = 0;
        long t1 = System.nanoTime();
        for (int mdIndex = 0; mdIndex < mds.size() - 2 /*skip unmapped*/; mdIndex++) {
            MarkDuplicatesInstance md1 = mds.get(mdIndex);
            endsCount += md1.pairedWithNoMate.size();
            for (int seqId = md1.leftmostSeqIndex; seqId <= md1.rightmostSeqIndex; seqId++) {
                mergeAllReadEndsFromMDInstanceOfSequence(mdIndex, seqId);
            }
        }
        log.debug("Merge time ", (System.nanoTime() - t1) / 1e9, " sec");
        log.debug(endsCount, " ends were merged.");
    }

    private void mergeAllReadEndsFromMDInstanceOfSequence(int mdIndex, int sequenceId) {
        MarkDuplicatesInstance md1 = mds.get(mdIndex);
        MarkDuplicatesInstance md2 = mds.get(mdIndex + 1);

        try (CloseableIterator<Map.Entry<String, ReadEndsForMarkDuplicates>> map1Iterator =
                     md1.pairedWithNoMate.iterator()) {

            while (map1Iterator.hasNext()) {
                Map.Entry<String, ReadEndsForMarkDuplicates> entry1 = map1Iterator.next();
                ReadEndsForMarkDuplicates end1 = entry1.getValue();
                String key1 = entry1.getKey();

                if (end1.read1ReferenceIndex != sequenceId) {
                    // skip all mates which are not from the current sequence
                    // note: iteration is in arbitrary order
                    continue;
                }

                if (!md2.containsSequenceInProcess(end1.read2ReferenceIndex)) {
                    md2 = pickMDForReferenceIndex(end1.read2ReferenceIndex);
                }
                ReadEndsForMarkDuplicates end2 =
                        md2.pairedWithNoMate.remove(end1.read1ReferenceIndex, key1);

                if (end2 == null) {
                    throw new IllegalStateException("Missing of a mate");
                }

                mergePairToBothMDInstances(md1, md2, end1, end2);
            }
        }
    }

    private MarkDuplicatesInstance pickMDForReferenceIndex(int seqIndex) {
        return mdsMap.floorEntry(seqIndex).getValue();
    }

    private void mergePairToBothMDInstances(MarkDuplicatesInstance md1, MarkDuplicatesInstance md2,
                                            ReadEndsForMarkDuplicates end1, ReadEndsForMarkDuplicates end2) {
        assert md1.leftmostSeqIndex < md2.leftmostSeqIndex;
        assert end1.read1ReferenceIndex < end2.read1ReferenceIndex;

        ReadEndsForMarkDuplicates mergedPair = combinePairEnds(end1, end2);

        ReadEndsForMarkDuplicates readEndsFor1;
        ReadEndsForMarkDuplicates readEndsFor2;

        if (useBarcodes) {
            readEndsFor1 = new ReadEndsForMarkDuplicatesWithBarcodes(
                    (ReadEndsForMarkDuplicatesWithBarcodes) mergedPair
            );
            readEndsFor2 = new ReadEndsForMarkDuplicatesWithBarcodes(
                    (ReadEndsForMarkDuplicatesWithBarcodes) mergedPair
            );
        } else {
            readEndsFor1 = mergedPair; // avoid creating excess objects
            readEndsFor2 = mergedPair.clone();
        }

        // read end #2 will not be processed by md1
        readEndsFor1.read2IndexInFile += (md2.getRecordsInFileBefore() - md1.getRecordsInFileBefore());
        md1.pairSort.add(readEndsFor1);

        // read end #1 with negative file indexes will be skipped by md2
        readEndsFor2.read1IndexInFile -= (md2.getRecordsInFileBefore() - md1.getRecordsInFileBefore());
        md2.pairSort.add(readEndsFor2);
    }

    private CompletableFuture<Void> writeBamFilesParallel() {
        CompletableFuture<Void> writeBamsFeature = CompletableFuture.completedFuture(null);
        for (MarkDuplicatesInstance md : mds) {
            writeBamsFeature = CompletableFuture.allOf(
                    writeBamsFeature,
                    CompletableFuture.runAsync(
                            () -> {
                                md.generateDuplicateIndexes(this.REMOVE_SEQUENCING_DUPLICATES
                                        || this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag);
                                md.writeFileAndCalculateMetric();
                            }, threadPool
                    )
            );
        }

        return writeBamsFeature;
    }

    private void mergeMetrics() {
        for (MarkDuplicatesInstance md : mds) {
            LibraryIdGenerator partLibIdGenerator = md.libraryIdGenerator;
            for (Map.Entry<String, DuplicationMetrics> toMergeMetricsByLibrary : partLibIdGenerator.getMetricsByLibraryMap().entrySet()) {
                String lib = toMergeMetricsByLibrary.getKey();
                DuplicationMetrics toMergeMetrics = toMergeMetricsByLibrary.getValue();

                libraryIdGenerator.getLibraryIdsMap().putIfAbsent(lib, partLibIdGenerator.getLibraryIdsMap().get(lib));

                DuplicationMetrics metrics = libraryIdGenerator.getMetricsByLibrary(lib);
                if (metrics == null) {
                    libraryIdGenerator.addMetricsByLibrary(lib, toMergeMetrics);
                } else {
                    metrics.merge(toMergeMetrics);
                }

            }

            libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap().addHistogram(
                    partLibIdGenerator.getOpticalDuplicatesByLibraryIdMap()
            );
        }
    }

    private void mergeOutputs() {
        List<File> parts = mds.stream()
                .map(md -> md.OUTPUT)
                .collect(Collectors.toList());

        if (BamFileIoUtils.isBamFile(OUTPUT)) {
            BamFileIoUtils.gatherWithBlockCopying(parts, OUTPUT, CREATE_INDEX, CREATE_MD5_FILE);
        } else {
            mergeSamFiles(parts);
        }

        mds.forEach(md -> md.OUTPUT.delete());
    }

    private class MarkDuplicatesInstance extends MarkDuplicates {

        // to safely skip all duplicate indexes from previous part
        public static final int DELIMITER_DUPLICATE_INDEX = -1;

        private List<SAMSequenceRecord> sequences;
        private int recordsInFileBefore;

        private final int leftmostSeqIndex;
        private final int rightmostSeqIndex;

        // to safely skip all ReadIndexers from previous part
        private final RepresentativeReadIndexer delimiterReadIndexer = new RepresentativeReadIndexer();
        {
            delimiterReadIndexer.setSize = DELIMITER_DUPLICATE_INDEX;
            delimiterReadIndexer.readIndexInFile = DELIMITER_DUPLICATE_INDEX;
            delimiterReadIndexer.representativeReadIndexInFile = DELIMITER_DUPLICATE_INDEX;
        }

        private MarkDuplicatesInstance() {
            //for unmapped part processing
            leftmostSeqIndex = Integer.MAX_VALUE;
            rightmostSeqIndex = Integer.MAX_VALUE;
        }

        MarkDuplicatesInstance(List<SAMSequenceRecord> sequences, int recordsInFileBefore) {
            this.sequences = sequences;
            this.recordsInFileBefore = recordsInFileBefore;

            this.leftmostSeqIndex = sequences.get(0).getSequenceIndex();
            this.rightmostSeqIndex = sequences.get(sequences.size() - 1).getSequenceIndex();
        }

        public List<SAMSequenceRecord> getSequences() {
            return sequences;
        }

        public int getRecordsInFileBefore() {
            return recordsInFileBefore;
        }

        public String fileNamePrefix() {
            String prefix;

            if (sequences == null) {
                prefix = "unmapped";
            } else if (sequences.size() == 1) {
                prefix = sequences.get(0).getSequenceName();
            } else {
                prefix = sequences.get(0).getSequenceName() + "-" + sequences.get(sequences.size() - 1).getSequenceName();
            }

            return prefix + "__";
        }

        public boolean containsSequenceInProcess(int index) {
            return index >= leftmostSeqIndex && index <= rightmostSeqIndex;
        }

        @Override
        protected int doWork() {
            try {
                OUTPUT = IOUtil.newTempFile(
                        fileNamePrefix(),
                        ParallelMarkDuplicates.this.OUTPUT.getName(),
                        TMP_DIR.toArray(new File[TMP_DIR.size()])
                );
                OUTPUT.deleteOnExit();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }

            // is used on fragSort and pairSort SortingCollections creation
            SORTING_COLLECTION_SIZE_RATIO *= SORTING_COLLECTION_SIZE_FACTOR / ParallelMarkDuplicates.this.PARALLELISM ;

            useBarcodes = ParallelMarkDuplicates.this.useBarcodes;

            return 0;
        }

        @Override
        int maxInMemoryForDuplicateIndexes(int entryOverhead) {
            return (int) Math.min(
                    (Runtime.getRuntime().maxMemory() * 0.25 * SORTING_COLLECTION_SIZE_FACTOR / PARALLELISM) / (double) entryOverhead,
                    (double) (Integer.MAX_VALUE - 5)
            );
        }

        @Override
        void freeUpSortingCollectionsMemory() {
            // Free up memory
            SortingCollectionUtil.freeUpMemory(pairSort);
            SortingCollectionUtil.freeUpMemory(fragSort);
            freeUpMemory(pairedWithNoMate);

            /* will call doneAdding later before generating of duplicate indices */
        }

        @Override
        void generateDuplicateIndexes(boolean indexOpticalDuplicates) {
            // Tell these collections to free up memory if possible.
            this.pairSort.doneAdding();
            this.fragSort.doneAdding();

            super.generateDuplicateIndexes(indexOpticalDuplicates);
        }

        @Override
        void writeFileAndCalculateMetric() {
            skipDuplicateIndexesFromPreviousPart();
            super.writeFileAndCalculateMetric();

            closeResources();
        }

        @Override
        CloseableIterator<RepresentativeReadIndexer> representativeReadIndicesIterator() {
            CloseableIterator<RepresentativeReadIndexer> iterator = skipRepresentativeEndsFromPrevPart();
            return new RepresentativeReadIndicesIteratorWrapper(iterator);
        }

        private void skipDuplicateIndexesFromPreviousPart() {
            while (duplicateIndexes.hasNext()) {
                if (duplicateIndexes.next() == DELIMITER_DUPLICATE_INDEX) {
                    break;
                }
            }
        }

        private CloseableIterator<RepresentativeReadIndexer> skipRepresentativeEndsFromPrevPart() {
            CloseableIterator<RepresentativeReadIndexer> iterator = representativeReadIndicesForDuplicates.iterator();
            //skip mates from previous part
            while (iterator.hasNext()) {
                RepresentativeReadIndexer next = iterator.next();
                if (next.representativeReadIndexInFile == DELIMITER_DUPLICATE_INDEX
                        && next.readIndexInFile == DELIMITER_DUPLICATE_INDEX
                        && next.setSize == DELIMITER_DUPLICATE_INDEX) {
                    break;
                }
            }
            return iterator;
        }

        @Override
        protected void doneAddingDuplicateIndexes() {
            duplicateIndexes.add(DELIMITER_DUPLICATE_INDEX);

            if (TAG_DUPLICATE_SET_MEMBERS) {
                representativeReadIndicesForDuplicates.add(delimiterReadIndexer);
            }

            super.doneAddingDuplicateIndexes();
        }

        @Override
        protected SamHeaderAndIterator openInputs(boolean eagerlyDecode) {
            fetchReaders(eagerlyDecode);

            final SamReader reader = readersByEagerlyDecode.get(eagerlyDecode).get(0);
            SAMFileHeader header = reader.getFileHeader();

            if (ASSUME_SORT_ORDER != null || ASSUME_SORTED) {
                if (ASSUME_SORT_ORDER == null) {
                    ASSUME_SORT_ORDER = SAMFileHeader.SortOrder.coordinate;
                    ASSUME_SORTED = false; // to maintain the "mutex" regarding these two arguments.
                }

                //if we assume a particular order, then the output will have that order in the header
                header.setSortOrder(ASSUME_SORT_ORDER);
            }

            SAMRecordIterator chrIterator = sequences != null ?
                    reader.query(
                            sequences.stream()
                                    .map(seq -> new QueryInterval(seq.getSequenceIndex(), 1, 0))
                                    .toArray(QueryInterval[]::new)
                            , true) :
                    reader.queryUnmapped();

            return new SamHeaderAndIterator(header, chrIterator);

        }

        @Override
        public void setupOpticalDuplicateFinder() {
            this.opticalDuplicateFinder = new OpticalDuplicateFinder(READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, MAX_OPTICAL_DUPLICATE_SET_SIZE, LOG) {
                //Not to count optical duplicates twice
                @Override
                public boolean[] findOpticalDuplicates(List<? extends PhysicalLocation> list, PhysicalLocation keeper) {
                    boolean[] flags = super.findOpticalDuplicates(list, keeper);

                    int flagIndex = 0;
                    for (PhysicalLocation physicalLocation : list) {
                        ReadEndsForMarkDuplicates end = (ReadEndsForMarkDuplicates) physicalLocation;

                        // we will not count ends from other parts as optical duplicates in this MD instance.
                        // But to tag them as optical (if needed) after cross-part paired ends merging,
                        // we set isOpticalDuplicate = true
                        if (!containsSequenceInProcess(end.read1ReferenceIndex)) {
                            if (flags[flagIndex]) {
                                flags[flagIndex] = false;
                                ((ReadEnds) physicalLocation).isOpticalDuplicate = true;
                            }
                        }
                        flagIndex++;
                    }

                    return flags;
                }
            };
        }

        /**
         * Corrects representativeReadIndexInFile by incrementing it's value by recordsInFileBefore
         * of the current MarkDuplicatesInstance.
         * Original value means index in the current part file, corrected value is an index in the result file.
         */
        private class RepresentativeReadIndicesIteratorWrapper implements CloseableIterator<RepresentativeReadIndexer> {

            private final CloseableIterator<RepresentativeReadIndexer> wrappedIterator;

            RepresentativeReadIndicesIteratorWrapper(CloseableIterator<RepresentativeReadIndexer> wrappedIterator) {
                this.wrappedIterator = wrappedIterator;
            }

            @Override
            public boolean hasNext() {
                return wrappedIterator.hasNext();
            }

            @Override
            public RepresentativeReadIndexer next() {
                RepresentativeReadIndexer next = wrappedIterator.next();
                next.representativeReadIndexInFile += recordsInFileBefore;

                return next;
            }

            @Override
            public void remove() {
                wrappedIterator.remove();
            }

            @Override
            public void close() {
                wrappedIterator.close();
            }

            @Override
            public List<RepresentativeReadIndexer> toList() {
                return wrappedIterator.toList();
            }

            @Override
            public Stream<RepresentativeReadIndexer> stream() {
                return wrappedIterator.stream();
            }
        }
    }

    private void mergeSamFiles(final List<File> parts) {
        final SAMFileHeader header = SamReaderFactory.makeDefault().getFileHeader(parts.get(0));

        try (SAMFileWriter out = new SAMFileWriterFactory()
                .setCreateIndex(CREATE_INDEX)
                .setCreateMd5File(CREATE_MD5_FILE)
                .makeSAMOrBAMWriter(header, true, OUTPUT)) {

            for (final File f : parts) {
                try (SamReader in = SamReaderFactory.makeDefault().open(f)) {
                    for (final SAMRecord rec : in) {
                        out.addAlignment(rec);
                    }
                }
            }

        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private static void freeUpMemory(DiskBasedReadEndsForMarkDuplicatesMap MDReadEndsMap) {
        MDReadEndsMap.remove(-1, "spillToDisk"); // will spill current sequence map to disk
    }
}
