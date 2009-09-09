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

package net.sf.picard.sam;

import java.io.*;
import java.util.*;

import net.sf.samtools.*;
import net.sf.samtools.util.*;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMValidationError.Type;
import net.sf.picard.PicardException;
import net.sf.picard.util.Histogram;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileWalker;

/**
 * Validates SAM files as follows:
 * <ul>
 * <li>checks sam file header for sequence dictionary</li>
 * <li>checks sam file header for read groups</li>
 * <li>for each sam record
 * <ul>
 *     <li>reports error detected by SAMRecord.isValid()</li>
 *     <li>validates NM (nucleotide differences) exists and matches reality</li>
 *     <li>validates mate fields agree with data in the mate record</li>
 * </ul>
 * </li>
 * </ul>
 *
 * @see SAMRecord#isValid()
 * @author Doug Voet
 */
public class SamFileValidator {
    private Histogram<Type> errorsByType = new Histogram<Type>();
    private final PrintWriter out;
    private Map<String, PairEndInfo> pairEndInfoByName;
    private ReferenceSequenceFileWalker refFileWalker = null;
    private boolean verbose = false;
    private int maxVerboseOutput = 100;
    private SAMRecordComparator recordComparator;
    private SortOrder sortOrder;
    private Set<Type> errorsToIgnore = EnumSet.noneOf(Type.class);
    private boolean ignoreWarnings = false;

    public SamFileValidator(final PrintWriter out) {
        this.out = out;
    }

    /** Sets one or more error types that should not be reported on. */
    public void setErrorsToIgnore(final Collection<Type> types) {
        if (!types.isEmpty()) {
            this.errorsToIgnore = EnumSet.copyOf(types);
        }
    }

    public void setIgnoreWarnings(final boolean ignoreWarnings) {
        this.ignoreWarnings = ignoreWarnings;
    }

    /**
     * Outputs validation summary report to out.
     * 
     * @param samReader records to validate
     * @param reference if null, NM tag validation is skipped
     */
    public void validateSamFileSummary(final SAMFileReader samReader, final ReferenceSequenceFile reference) {
        init(reference);

        validateSamFile(samReader, out);
        
        if (errorsByType.getCount() > 0) {
            // Convert to a histogram with String IDs so that WARNING: or ERROR: can be prepended to the error type.
            final Histogram<String> errorsAndWarningsByType = new Histogram<String>("Error Type", "Count");
            for (final Histogram<SAMValidationError.Type>.Bin bin : errorsByType.values()) {
                errorsAndWarningsByType.increment(bin.getId().getHistogramString(), bin.getValue());
            }
            final MetricsFile<ValidationMetrics, String> metricsFile = new MetricsFile<ValidationMetrics, String>();
            errorsByType.setBinLabel("Error Type");
            errorsByType.setValueLabel("Count");
            metricsFile.setHistogram(errorsAndWarningsByType);
            metricsFile.write(out);
        }
        cleanup();
    }

    /**
     * Outputs validation error details to out.
     * 
     * @param samReader records to validate
     * @param reference if null, NM tag validation is skipped
     * processing will stop after this threshold has been reached
     */
    public void validateSamFileVerbose(final SAMFileReader samReader, final ReferenceSequenceFile reference) {
        init(reference);

        try {
            validateSamFile(samReader, out);
        } catch (MaxOutputExceededException e) {
            out.println("Maximum output of [" + maxVerboseOutput + "] errors reached.");
        }
        
        cleanup();
    }

    public void validateBamFileTermination(final File inputFile) {
        BufferedInputStream inputStream = null;
        try {
            inputStream = IOUtil.toBufferedStream(new FileInputStream(inputFile));
            if (!BlockCompressedInputStream.isValidFile(inputStream)) {
                return;
            }
            final BlockCompressedInputStream.FileTermination terminationState =
                    BlockCompressedInputStream.checkTermination(inputFile);
            if (terminationState.equals(BlockCompressedInputStream.FileTermination.DEFECTIVE)) {
                addError(new SAMValidationError(Type.TRUNCATED_FILE, "BAM file has defective last gzip block",
                        inputFile.getPath()));
            } else if (terminationState.equals(BlockCompressedInputStream.FileTermination.HAS_HEALTHY_LAST_BLOCK)) {
                addError(new SAMValidationError(Type.BAM_FILE_MISSING_TERMINATOR_BLOCK,
                        "Older BAM file -- does not have terminator block",
                        inputFile.getPath()));

            }
        } catch (IOException e) {
            throw new PicardException("IOException", e);
        } finally {
            if (inputStream != null) {
                CloserUtil.close(inputStream);
            }
        }
    }

    private void validateSamFile(final SAMFileReader samReader, final PrintWriter out) {
        try {
            samReader.setValidationStringency(ValidationStringency.SILENT);
            validateHeader(samReader.getFileHeader());
            initComparator(samReader.getFileHeader());
            validateSamRecords(samReader);
            
            if (errorsByType.isEmpty()) {
                out.println("No errors found");
            }
        } finally {
            out.flush();
        }
    }

    private void initComparator(final SAMFileHeader fileHeader) {
        this.sortOrder = fileHeader.getSortOrder();
        switch (this.sortOrder) {
        case coordinate:
            this.recordComparator = new SAMRecordCoordinateComparator();
            break;
        case queryname:
            this.recordComparator = new SAMRecordQueryNameComparator();
            break;
        default:
            // dummy SAMRecordComparator that always says records are equal
            this.recordComparator = new SAMRecordComparator() {
                public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) { return 0; }
                public int compare(final SAMRecord o1, final SAMRecord o2) { return 0; }    
            };
            break;
        }
    }

    private void validateSamRecords(final Iterable<SAMRecord> samRecords) {
        long recordNumber = 1;
        SAMRecord lastRecord = null;
        try {
            for (final SAMRecord record : samRecords) {
                final Collection<SAMValidationError> errors = record.isValid();
                if (errors != null) {
                    for (final SAMValidationError error : errors) {
                        error.setRecordNumber(recordNumber);
                        addError(error);
                    }
                }
                
                validateMateFields(record, recordNumber);
                validateNmTag(record, recordNumber);
                validateSortOrder(lastRecord, record, recordNumber);
                validateCigar(record, recordNumber);
                validateSecondaryBaseCalls(record, recordNumber);

                recordNumber++;
                lastRecord = record;
            }
        } catch (SAMFormatException e) {
            // increment record number because the iterator behind the SAMFileReader
            // reads one record ahead so we will get this failure one record ahead
            out.println("SAMFormatException on record " + ++recordNumber);
            throw new PicardException("SAMFormatException on record " + recordNumber, e);
        } catch (FileTruncatedException e) {
            addError(new SAMValidationError(Type.TRUNCATED_FILE, "File is truncated", null));
        }
    }

    private void validateSecondaryBaseCalls(final SAMRecord record, final long recordNumber) {
        final String e2 = (String)record.getAttribute(SAMTag.E2.name());
        if (e2 != null) {
            if (e2.length() != record.getReadLength()) {
                addError(new SAMValidationError(Type.MISMATCH_READ_LENGTH_AND_E2_LENGTH,
                        String.format("E2 tag length (%d) != read length (%d)", e2.length(), record.getReadLength()),
                        record.getReadName(), recordNumber));
            }
            final byte[] bases = record.getReadBases();
            final byte[] secondaryBases = StringUtil.stringToBytes(e2);
            for (int i = 0; i < Math.min(bases.length, secondaryBases.length); ++i) {
                if (SequenceUtil.isNoCall(bases[i]) || SequenceUtil.isNoCall(secondaryBases[i])) {
                    continue;
                }
                if (SequenceUtil.basesEqual(bases[i], secondaryBases[i])) {
                    addError(new SAMValidationError(Type.E2_BASE_EQUALS_PRIMARY_BASE,
                            String.format("Secondary base call  (%c) == primary base call (%c)",
                                    (char)secondaryBases[i], (char)bases[i]),
                            record.getReadName(), recordNumber));
                    break;
                }
            }
        }
        final String u2 = (String)record.getAttribute(SAMTag.U2.name());
        if (u2 != null && u2.length() != record.getReadLength()) {
            addError(new SAMValidationError(Type.MISMATCH_READ_LENGTH_AND_U2_LENGTH,
                    String.format("U2 tag length (%d) != read length (%d)", u2.length(), record.getReadLength()),
                    record.getReadName(), recordNumber));
        }
    }

    private void validateCigar(final SAMRecord record, final long recordNumber) {
        if (record.getReadUnmappedFlag()) {
            return;
        }
        final ValidationStringency savedStringency = record.getValidationStringency();
        record.setValidationStringency(ValidationStringency.LENIENT);
        final List<SAMValidationError> errors = record.validateCigar(recordNumber);
        record.setValidationStringency(savedStringency);
        if (errors == null) {
            return;
        }
        for (final SAMValidationError error : errors) {
            addError(error);
        }
    }

    private void validateSortOrder(final SAMRecord lastRecord, final SAMRecord record, final long recordNumber) {
        if (lastRecord != null && recordComparator.fileOrderCompare(lastRecord, record) > 0) {
            addError(new SAMValidationError(
                    Type.RECORD_OUT_OF_ORDER, 
                    String.format(
                            "The record is out of [%s] order, prior read name [%s], prior coodinates [%d:%d]",
                            this.sortOrder.name(),
                            lastRecord.getReadName(),
                            lastRecord.getReferenceIndex(),
                            lastRecord.getAlignmentStart()),
                    record.getReadName(), 
                    recordNumber));
        }
    }
    
    private void init(final ReferenceSequenceFile reference) {
        this.pairEndInfoByName = new HashMap<String, PairEndInfo>();
        if (reference != null) {
            this.refFileWalker = new ReferenceSequenceFileWalker(reference);
        }
    }

    private void cleanup() {
        this.errorsByType = null;
        this.pairEndInfoByName = null;
        this.refFileWalker = null;
    }

    private void validateNmTag(final SAMRecord record, final long recordNumber) {
        if (!record.getReadUnmappedFlag()) {
            final Integer tagNucleotideDiffs = record.getIntegerAttribute(ReservedTagConstants.NM);
            if (tagNucleotideDiffs == null) {
                addError(new SAMValidationError(
                        Type.MISSING_TAG_NM, 
                        "NM tag (nucleotide differences) is missing", 
                        record.getReadName(),
                        recordNumber));
            } else if (refFileWalker != null) {
                final ReferenceSequence refSequence = refFileWalker.get(record.getReferenceIndex());
                final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSequence.getBases());

                if (!tagNucleotideDiffs.equals(actualNucleotideDiffs)) {
                    addError(new SAMValidationError(
                            Type.INVALID_TAG_NM, 
                            "NM tag (nucleotide differences) in file [" + tagNucleotideDiffs +
                            "] does not match reality [" + actualNucleotideDiffs + "]", 
                            record.getReadName(),
                            recordNumber));
                }
            }
        }
    }
    
    private void validateMateFields(final SAMRecord record, final long recordNumber) {
        if (!record.getReadPairedFlag()) {
            return;
        }
        
        final PairEndInfo pairEndInfo = pairEndInfoByName.remove(record.getReadName());
        if (pairEndInfo == null) {
            pairEndInfoByName.put(record.getReadName(), new PairEndInfo(record, recordNumber));
        } else {
            final List<SAMValidationError> errors =
                pairEndInfo.validateMates(new PairEndInfo(record, recordNumber), record.getReadName());
            for (final SAMValidationError error : errors) {
                addError(error);
            }
        }
    }

    private void validateHeader(final SAMFileHeader fileHeader) {
        if (fileHeader.getVersion() == null) {
            addError(new SAMValidationError(Type.MISSING_VERSION_NUMBER, "Header has no version number", null));
        } else if (!fileHeader.getVersion().equals(SAMFileHeader.CURRENT_VERSION)) {
            addError(new SAMValidationError(Type.INVALID_VERSION_NUMBER, "Header version: " +
                    fileHeader.getVersion() + " does not match expected version: " + SAMFileHeader.CURRENT_VERSION,
                    null));
        }
        if (fileHeader.getSequenceDictionary().isEmpty()) {
            addError(new SAMValidationError(Type.MISSING_SEQUENCE_DICTIONARY, "Sequence dictionary is empty", null));
        }
        if (fileHeader.getReadGroups().isEmpty()) {
            addError(new SAMValidationError(Type.MISSING_READ_GROUP, "Read groups is empty", null));
        }
    }

    private void addError(final SAMValidationError error) {
        // Just ignore an error if it's of a type we're not interested in
        if (this.errorsToIgnore.contains(error.getType())) return;

        if (this.ignoreWarnings && error.getType().severity == SAMValidationError.Severity.WARNING) return;

        this.errorsByType.increment(error.getType());
        if (verbose) {
            out.println(error);
            if (this.errorsByType.getCount() >= maxVerboseOutput) {
                throw new MaxOutputExceededException();
            }
        }
    }

    /**
     * Control verbosity
     * @param verbose True in order to emit a message per error or warning.
     * @param maxVerboseOutput If verbose, emit no more than this many messages.  Ignored if !verbose.
     */
    public void setVerbose(final boolean verbose, final int maxVerboseOutput) {
        this.verbose = verbose;
        this.maxVerboseOutput = maxVerboseOutput;
    }

    public static class ValidationMetrics extends MetricBase {
    }
    
    /** 
     * This class is used so we don't have to store the entire SAMRecord in memory while we wait
     * to find a record's mate and also to store the record number.
     */
    private static class PairEndInfo {
        private final int readAlignmentStart;
        private final int readReferenceIndex;
        private final boolean readNegStrandFlag;
        private final boolean readUnmappedFlag;

        private final int mateAlignmentStart;
        private final int mateReferenceIndex;
        private final boolean mateNegStrandFlag;
        private final boolean mateUnmappedFlag;
        
        private final long recordNumber;
        
        public PairEndInfo(final SAMRecord record, final long recordNumber) {
            this.recordNumber = recordNumber;
            
            this.readAlignmentStart = record.getAlignmentStart();
            this.readNegStrandFlag = record.getReadNegativeStrandFlag();
            this.readReferenceIndex = record.getReferenceIndex();
            this.readUnmappedFlag = record.getReadUnmappedFlag();
            
            this.mateAlignmentStart = record.getMateAlignmentStart();
            this.mateNegStrandFlag = record.getMateNegativeStrandFlag();
            this.mateReferenceIndex = record.getMateReferenceIndex();
            this.mateUnmappedFlag = record.getMateUnmappedFlag();
        }
        
        public List<SAMValidationError> validateMates(final PairEndInfo mate, final String readName) {
            final List<SAMValidationError> errors = new ArrayList<SAMValidationError>();
            validateMateFields(this, mate, readName, errors);
            validateMateFields(mate, this, readName, errors);
            return errors;
        }
        
        private void validateMateFields(final PairEndInfo end1, final PairEndInfo end2, final String readName, final List<SAMValidationError> errors) {
            if (end1.mateAlignmentStart != end2.readAlignmentStart) {
                errors.add(new SAMValidationError(
                        Type.MISMATCH_MATE_ALIGNMENT_START, 
                        "Mate alignment does not match alignment start of mate", 
                        readName,
                        end1.recordNumber));
            }
            if (end1.mateNegStrandFlag != end2.readNegStrandFlag) {
                errors.add(new SAMValidationError(
                        Type.MISMATCH_FLAG_MATE_NEG_STRAND, 
                        "Mate negative strand flag does not match read negative strand flag of mate", 
                        readName,
                        end1.recordNumber));
            }
            if (end1.mateReferenceIndex != end2.readReferenceIndex) {
                errors.add(new SAMValidationError(
                        Type.MISMATCH_MATE_REF_INDEX, 
                        "Mate reference index (MRNM) does not match reference index of mate", 
                        readName,
                        end1.recordNumber));
            }
            if (end1.mateUnmappedFlag != end2.readUnmappedFlag) {
                errors.add(new SAMValidationError(
                        Type.MISMATCH_FLAG_MATE_UNMAPPED, 
                        "Mate unmapped flag does not match read unmapped flag of mate", 
                        readName,
                        end1.recordNumber));
            }
        }
    }
    
    /** Thrown in addError indicating that maxVerboseOutput has been exceeded and processing should stop */
    private static class MaxOutputExceededException extends RuntimeException { }
}
