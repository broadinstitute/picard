/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.EdgingRecordAndOffset;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SequenceUtil;

import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.HashMap;

/**
 * Class represents fast algorithm for collecting data from <code>AbstractLocusInfo</code>
 * with a list of aligned {@link htsjdk.samtools.util.EdgingRecordAndOffset} objects. According to the algorithm
 * we receive only two {@link htsjdk.samtools.util.EdgingRecordAndOffset} objects for each alignment block of a read:
 * one for the start of block and one for the end. When meeting a {@link htsjdk.samtools.util.EdgingRecordAndOffset}
 * with type {@link htsjdk.samtools.util.EdgingRecordAndOffset.Type#BEGIN}, all information from the alignment block is accumulated in the collector,
 * read name is added to a map holding the names of processed reads for detecting overlapping positions.
 * When meeting a {@link htsjdk.samtools.util.EdgingRecordAndOffset} with type {@link htsjdk.samtools.util.EdgingRecordAndOffset.Type#END},
 * the read name is removed from the map with names of processed reads.
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */

public class FastWgsMetricsCollector extends AbstractWgsMetricsCollector<EdgingRecordAndOffset> {

    /**
     * Index of the sequence of the previous processed <code>AbstractLocusInfo</code>
     */
    private int previousSequenceIndex;

    /**
     * Manager for {@link this#pileupSize} Counter.
     */
    private final CounterManager counterManager;
    /**
     * Counter, accumulating the data on coverage for the calculation of base quality histogram.
     */
    private final CounterManager.Counter pileupSize;
    /**
     * Counter, accumulating the data on unfiltered coverage (includes all but quality 2 bases)
     * for the calculation of theoretical het sensitivity.
     */
    private final CounterManager.Counter unfilteredDepthSize;

    /**
     * Map, holding information on currently processed reads, that possibly have overlapping regions.
     */
    private Map<String, Set<? extends AbstractRecordAndOffset>> readsNames;

    /**
     * Determines the size of created {@link picard.analysis.CounterManager.Counter} objects. The bigger {@link picard.analysis.CounterManager.Counter} objects
     * are created, the less rebasing in {@link picard.analysis.CounterManager.Counter} will occur.
     */
    private final int ARRAY_SIZE_PER_READ_LENGTH = 2000;

    /**
     * Creates a collector and initializes the inner data structures
     *
     * @param collectWgsMetrics CollectWgsMetrics, that creates this collector
     * @param coverageCap       coverage cap
     */
    public FastWgsMetricsCollector(CollectWgsMetrics collectWgsMetrics, int coverageCap, final IntervalList intervals) {
        super(collectWgsMetrics, coverageCap, intervals);
        this.previousSequenceIndex = -1;
        this.counterManager = new CounterManager(collectWgsMetrics.READ_LENGTH * ARRAY_SIZE_PER_READ_LENGTH, collectWgsMetrics.READ_LENGTH);
        this.pileupSize = counterManager.newCounter();
        this.unfilteredDepthSize = counterManager.newCounter();
    }

    @Override
    public void addInfo(final AbstractLocusInfo<EdgingRecordAndOffset> info, final ReferenceSequence ref, boolean referenceBaseN) {
        prepareCollector(info);
        for (final EdgingRecordAndOffset record : info.getRecordAndOffsets()) {
            final String readName = record.getReadName();
            Optional<Set<EdgingRecordAndOffset>> recordsAndOffsetsForName = Optional.ofNullable((Set<EdgingRecordAndOffset>)readsNames.get(readName));
            if (record.getType() == EdgingRecordAndOffset.Type.BEGIN) {
                processRecord(info.getPosition(), ref, record, recordsAndOffsetsForName.orElse(new HashSet<>()));
            } else {
                recordsAndOffsetsForName.ifPresent(
                        edgingRecordAndOffsets -> removeRecordFromMap(record,
                                edgingRecordAndOffsets));
            }
        }
        if (!referenceBaseN) {
            final int readNamesSize = pileupSize.get(info.getPosition());
            final int highQualityDepth = Math.min(readNamesSize, coverageCap);
            if (highQualityDepth < readNamesSize) {
                basesExcludedByCapping += readNamesSize - coverageCap;
            }
            highQualityDepthHistogramArray[highQualityDepth]++;
            unfilteredDepthHistogramArray[unfilteredDepthSize.get(info.getPosition())]++;
        }
    }

    private void processRecord(int position, ReferenceSequence ref, EdgingRecordAndOffset record, Set<EdgingRecordAndOffset> recordsAndOffsetsForName) {
        long processedLoci = counter;
        readsNames.put(record.getReadName(), recordsAndOffsetsForName);
        final byte[] qualities = record.getBaseQualities();
        final byte[] bases = record.getRecord().getReadBases();
        for (int i = 0; i < record.getLength(); i++) {
            final int index = i + record.getRefPos();
            if (isReferenceBaseN(index, ref)) {
                continue;
            }
            final byte quality = qualities[i + record.getOffset()];
            if (quality <= 2) {
                basesExcludedByBaseq++;
            } else {
                if (unfilteredDepthSize.get(index) < coverageCap) {
                    unfilteredBaseQHistogramArray[quality]++;
                    unfilteredDepthSize.increment(index);
                }
                if (quality < collectWgsMetrics.MINIMUM_BASE_QUALITY || SequenceUtil.isNoCall(bases[i + record.getOffset()])){
                    basesExcludedByBaseq++;
                } else {
                    final int bsq = excludeByQuality(recordsAndOffsetsForName, index);
                    if (recordsAndOffsetsForName.size() - bsq > 0) {
                        basesExcludedByOverlap++;
                    } else {
                        pileupSize.increment(index);
                    }
                }
            }
            if (isTimeToStop(++processedLoci)) {
                break;
            }
        }
        recordsAndOffsetsForName.add(record);
    }

    private void removeRecordFromMap(EdgingRecordAndOffset record, Set<EdgingRecordAndOffset> recordsAndOffsetsForName) {
        if (recordsAndOffsetsForName.size() == 1) {
            readsNames.remove(record.getReadName());
        } else {
            recordsAndOffsetsForName.remove(record.getStart());
        }
    }

    /**
     * Prepares the accumulator objects to process a new {@link htsjdk.samtools.util.AbstractLocusInfo}.
     * If we switch to a new sequence, all accumulators are cleared.
     *
     * @param info the next {@link htsjdk.samtools.util.AbstractLocusInfo} to process
     */
    private void prepareCollector(AbstractLocusInfo<EdgingRecordAndOffset> info) {
        if (readsNames == null) {
            readsNames = new HashMap<>();
        }
        if (previousSequenceIndex != info.getSequenceIndex()) {
            readsNames.clear();
            counterManager.clear();
            previousSequenceIndex = info.getSequenceIndex();
        }
        counterManager.checkOutOfBounds(info.getPosition());
    }

    private int excludeByQuality(final Set<EdgingRecordAndOffset> setForName, int position) {
        int bsq = 0;
        for (EdgingRecordAndOffset recordAndOffset : setForName) {
            if (position - recordAndOffset.getRefPos() >= recordAndOffset.getLength()
                    || recordAndOffset.getBaseQuality(position) < collectWgsMetrics.MINIMUM_BASE_QUALITY) {
                bsq++;
            }
        }
        return bsq;
    }
}
