/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.illumina.parser;

import net.sf.picard.illumina.BarcodeUtil;

import java.io.File;
import java.io.FilenameFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import net.sf.picard.PicardException;
import net.sf.picard.util.Log;
import net.sf.samtools.util.AsciiLineReader;

/**
 * Class for specifying options for parsing Illumina basecall output files for a lane, and then creating
 * an IlluminaDataProvider for iterating over the reads in that lane.
 * 
 * @author alecw@broadinstitute.org
 */
public class IlluminaDataProviderFactory {
    private static final Log log = Log.getInstance(IlluminaDataProviderFactory.class);

    public enum BaseCallerVersion { Bustard_1_1, Bustard_1_3, Bustard_1_4, rta, Bustard_1_5, Bustard_1_6 }
    public enum ImageAnalyzerVersion {Firecrest_1_1, Firecrest_1_3, Firecrest_1_4, rta}

    // The following properties must be specified by caller.
    private final File basecallDirectory;
    private final int lane;

    // For all current known Bustard versions, this is not available in the output files therefore
    // it must be set by the caller in the ctor.
    private final Integer barcodeCycle;
    private final Integer barcodeLength;

    // The kinds of data the caller is interested in.  In some cases, more data types
    // than requested are returned, because of the granularity of the data types in files.
    private final Set<IlluminaDataType> dataTypes =
            EnumSet.noneOf(IlluminaDataType.class);

    // The following properties are automatically set, but may be overridden by caller.
    private File rawIntensityDirectory;
    private Boolean pairedEnd;
    private BaseCallerVersion baseCallerVersion;
    private ImageAnalyzerVersion imageAnalyzerVersion;

    // Set to true after computeReadConfiguration() method is called
    private boolean prepared = false;

    // In general the first end must be non-zero length, but occasionally an index-only run is done, i.e. all that
    // is sequenced are the barcodes.
    private boolean allowZeroLengthFirstEnd = false;

    // The final configuration passed to each parser.
    private final ReadConfiguration readConfiguration = new ReadConfiguration();

    /**
     * Prepare to iterate over non-barcoded Illumina Basecall output.
     *
     * @param basecallDirectory Where the Illumina basecall output is.  Some files are found by looking relative to this
     *                          directory, at least by default.
     * @param lane              Which lane to iterate over.
     * @param dataTypes         Which data types to read.  Note that some additional data types may be fetched
     *                          if they reside in the same file type as a requested data type.
     */
    public IlluminaDataProviderFactory(final File basecallDirectory, final int lane,
                                final IlluminaDataType... dataTypes) {
        this.basecallDirectory = basecallDirectory;
        this.lane = lane;
        this.barcodeCycle = null;
        this.barcodeLength = null;
        this.dataTypes.addAll(Arrays.asList(dataTypes));
    }

    /**
     * Prepare to iterate over the Illumina basecall output with barcodes.
     *
     * @param basecallDirectory Where the basecall output is.  Some files are found by looking relative to this
     *                          directory, at least by default.
     * @param lane              Which lane to iterate over.
     * @param barcodeCycle 1-base cycle number where barcode starts
     * @param barcodeLength number of cycles in barcode
     * @param dataTypes         Which data types to read.  Note that some additional data types may be fetched
     *                          if they reside in the same file type as a requested data type.
     */
    public IlluminaDataProviderFactory(final File basecallDirectory, final int lane,
                                       final int barcodeCycle, final int barcodeLength,
                                final IlluminaDataType... dataTypes) {
        this.basecallDirectory = basecallDirectory;
        this.lane = lane;
        if (barcodeCycle < 1) {
            throw new IllegalArgumentException("barcodeCycle is < 1: " + barcodeCycle);
        }
        if (barcodeLength < 1) {
            throw new IllegalArgumentException("barcodeLength is < 1: " + barcodeLength);
        }
        this.barcodeCycle = barcodeCycle;
        this.barcodeLength = barcodeLength;
        this.dataTypes.addAll(Arrays.asList(dataTypes));
    }

    /**
     * Return the list of tiles available for this flowcell and lane.  These are in ascending numerical order.
     * Note: Triggers computeReadConfiguration() if not already called.
     *
     * @return List of all tiles available for this flowcell and lane.
     */
    public List<Integer> getTiles() {
        if (!prepared) {
            computeReadConfiguration();
        }
        final TiledIlluminaFile[] tiledFiles;
        if (baseCallerVersion == BaseCallerVersion.Bustard_1_1) {
            tiledFiles = IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(basecallDirectory, "seq", lane);
        } else {
            tiledFiles = IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, 1);
        }
        final List<Integer> ret = new ArrayList<Integer>(tiledFiles.length);
        for (final TiledIlluminaFile tiledFile : tiledFiles) {
            ret.add(tiledFile.tile);
        }
        return ret;
    }

    /**
     * After calling any desired setters to override default detection, call this method to create
     * an iterator.  Tiles are iterated over in ascending numeric order.
     * Note: Triggers computeReadConfiguration() if not already called.
     *
     * @return An iterator for reading the Illumina basecall output for the lane specified in the ctor.
     */
    public IlluminaDataProvider makeDataProvider() {
        return makeDataProvider(null);
    }
    /**
     * After calling any desired setters to override default detection, call this method to create
     * an iterator.
     * Note: Triggers computeReadConfiguration() if not already called.
     *
     * @param tiles The tiles to iterate over, in order, or null to get all tiles in ascending numerical order.
     * @return An iterator for reading the Illumina basecall output for the lane specified in the ctor.
     */
    public IlluminaDataProvider makeDataProvider(final List<Integer> tiles) {
        if (dataTypes.isEmpty()) {
            throw new PicardException("No data types have been specified for basecall output " + basecallDirectory +
            ", lane " + lane);
        }
        if (!prepared) {
            computeReadConfiguration();
        }
        boolean madeQseqParser = false;
        final List<IlluminaParser> parsers = new ArrayList<IlluminaParser>();
        for (final IlluminaDataType dataType : dataTypes) {
            switch (dataType) {
                case Barcodes:
                    parsers.add(new BarcodeParser(readConfiguration, basecallDirectory, lane, tiles));
                    break;
                case BaseCalls:
                    if (baseCallerVersion == BaseCallerVersion.Bustard_1_1) {
                        parsers.add(new SeqParser(readConfiguration, basecallDirectory, lane, tiles));
                    } else if (!madeQseqParser) {
                        makeQseqParsers(parsers, tiles);
                        madeQseqParser = true;
                    }
                    break;
                case Noise:
                    if (imageAnalyzerVersion == ImageAnalyzerVersion.rta) {
                        parsers.add(new CnfParser(readConfiguration, rawIntensityDirectory, lane, tiles));
                    } else {
                        parsers.add(new NseParser(readConfiguration, rawIntensityDirectory, lane, tiles));
                    }
                    break;
                case PF:
                    if (baseCallerVersion == BaseCallerVersion.Bustard_1_1) {
                        parsers.add(new QhgParser(readConfiguration, basecallDirectory, lane, tiles));
                    } else if (!madeQseqParser) {
                        makeQseqParsers(parsers, tiles);
                        madeQseqParser = true;
                    }
                    break;
                case ProcessedIntensities:
                    parsers.add(new Sig2Parser(readConfiguration, basecallDirectory, lane, tiles));
                    break;
                case QualityScores:
                    if (baseCallerVersion == BaseCallerVersion.Bustard_1_1) {
                        parsers.add(new PrbParser(readConfiguration, basecallDirectory, lane, tiles));
                    } else if (!madeQseqParser) {
                        makeQseqParsers(parsers, tiles);
                        madeQseqParser = true;
                    }
                    break;
                case RawIntensities:
                    if (imageAnalyzerVersion == ImageAnalyzerVersion.rta) {
                        parsers.add(new CifParser(readConfiguration, rawIntensityDirectory, lane, tiles));
                    } else {
                        parsers.add(new IntParser(readConfiguration, rawIntensityDirectory, lane, tiles));
                    }
                    break;
            }
        }
        return new IlluminaDataProvider(readConfiguration.isBarcoded(), readConfiguration.isPairedEnd(), parsers,
                basecallDirectory, lane);
    }


    private void makeQseqParsers(final List<IlluminaParser> parsers, final List<Integer> tiles) {
        if (isBarcodeAwareBaseCaller()) {
            // Barcode is in a separate qseq file.
            int firstEndFileNumber = 1;
            int secondEndFileNumber = 2;
            int barcodeFileNumber = -1;
            if (readConfiguration.isBarcoded()) {
                if (readConfiguration.getBarcodeStart() < readConfiguration.getFirstStart()) {
                    barcodeFileNumber = 1;
                    firstEndFileNumber = 2;
                    secondEndFileNumber = 3;
                } else if (!readConfiguration.isPairedEnd()) {
                    barcodeFileNumber = 2;
                } else if (readConfiguration.getBarcodeStart() < readConfiguration.getSecondStart()) {
                    barcodeFileNumber = 2;
                    secondEndFileNumber = 3;
                } else {
                    barcodeFileNumber = 3;
                }
            }
            parsers.add(new QseqParser(readConfiguration, basecallDirectory, lane, EndType.FIRST, firstEndFileNumber, tiles));
            if (readConfiguration.isPairedEnd()) {
                parsers.add(new QseqParser(readConfiguration, basecallDirectory, lane, EndType.SECOND, secondEndFileNumber, tiles));
            }
            if (readConfiguration.isBarcoded()) {
                parsers.add(new QseqParser(readConfiguration, basecallDirectory, lane, EndType.BARCODE, barcodeFileNumber, tiles));
            }
        } else {
            // Barcode is in one of the regular qseq files.
            parsers.add(new QseqParser(readConfiguration, basecallDirectory, lane, EndType.FIRST, tiles));
            // Handle case in which the barcode is the only thing in read 2, typically because the run was stopped
            // after reading the barcode but before read 2.
            if (pairedEnd || (readConfiguration.isBarcoded() && readConfiguration.getBarcodeRead() == EndType.SECOND)) {
                parsers.add(new QseqParser(readConfiguration, basecallDirectory, lane, EndType.SECOND, tiles));
            }
        }
    }

    /**
     * Check if intensities are available for the given tile.  Note that a call to computeReadConfiguration()
     * is triggered if it has not already been invoked.
     * @param tile For which tile intensities are being sought.
     * @return True if intensities are available for that tile.
     */
    public boolean intensitiesAvailable(final int tile) {
        if (!prepared) {
            computeReadConfiguration();
        }
        if (imageAnalyzerVersion == ImageAnalyzerVersion.rta) {
            return CifParser.cifExists(rawIntensityDirectory, lane, tile);

        } else {
            return IntParser.intExists(rawIntensityDirectory, lane, tile);
        }
    }

    /**
     * This method may be called prior to calling makeDataProvider(), in order to trigger auto-detection
     * of ReadConfiguration.
     */
    public void computeReadConfiguration() {
        if (prepared) {
            throw new IllegalStateException("Already prepared");
        }
        if (baseCallerVersion == null || imageAnalyzerVersion == null) {
            detectPipelineVersion();
        }

        if (rawIntensityDirectory == null) {
            switch (imageAnalyzerVersion) {
                // Currently these are all the same
                case Firecrest_1_1:
                case Firecrest_1_3:
                case Firecrest_1_4:
                case rta:
                    rawIntensityDirectory = basecallDirectory.getParentFile();
                    break;
                default:
                    throw new IllegalStateException("Could not determine image analyzer version for " +
                            basecallDirectory + "; lane " + lane);
            }
        }

        boolean qseqIsBarcodeAware = false;
        switch (baseCallerVersion) {
            case rta:
            case Bustard_1_3:
            case Bustard_1_4:
                computeReadConfigurationFromBarcodeUnawareQseq();
                break;
            case Bustard_1_5:
            case Bustard_1_6:
                computeReadConfigurationFromBarcodeAwareQseq();
                qseqIsBarcodeAware = true;
                break;
            case Bustard_1_1:
                computeReadConfiguationFrom_1_1();
                break;
            default:
                throw new IllegalStateException("Could not determine base caller version for " + basecallDirectory + "; lane " + lane);
        }
        if (!qseqIsBarcodeAware) {
            updateReadConfigurationForBarcode();
        }
        readConfiguration.assertValid(allowZeroLengthFirstEnd);
        prepared = true;
    }


    private void computeReadConfigurationFromBarcodeUnawareQseq() {
        if (pairedEnd == null) {
            pairedEnd = IlluminaFileUtil.endedIlluminaBasecallFilesExist(basecallDirectory, "qseq", lane, 2);
        }
        readConfiguration.setPairedEnd(pairedEnd);
        if (readConfiguration.getFirstLength() == 0) {
            readConfiguration.setFirstStart(1);
            final File firstEnd = IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, 1)[0].file;
            readConfiguration.setFirstEnd(QseqParser.getReadLength(firstEnd));
        }
        if (pairedEnd && readConfiguration.getSecondLength() == 0) {
            final File secondEnd = IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, 2)[0].file;
            final int end2Length = QseqParser.getReadLength(secondEnd);
            readConfiguration.setSecondStart(readConfiguration.getFirstLength() + 1);
            readConfiguration.setSecondEnd(readConfiguration.getFirstLength() + end2Length);
        }
    }


    private void computeReadConfigurationFromBarcodeAwareQseq() {
        final File[] qseqs = new File[3];

        qseqs[0] = IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, 1)[0].file;

        for (int end = 2; end <= 3; ++end) {
            final TiledIlluminaFile[] files = IlluminaFileUtil.getEndedIlluminaBasecallFiles(basecallDirectory, "qseq", lane, end);
            qseqs[end-1] = (files.length > 0? files[0].file: null);
        }

        readConfiguration.setBarcoded(barcodeCycle != null);
        int numQseqs = readConfiguration.isBarcoded()? 2: 1;

        // Determine if paired end
        if (pairedEnd == null) {
            pairedEnd = qseqs[numQseqs] != null;
        }
        if (pairedEnd) {
            ++numQseqs;
        }
        readConfiguration.setPairedEnd(pairedEnd);

        // Read lengths from all qseqs
        final int[] qseqLengths = new int[numQseqs];
        for (int i = 0; i < numQseqs; ++i) {
            qseqLengths[i] = QseqParser.getReadLength(qseqs[i]);
        }

        // Determine order of barcode, first end & second end
        // Which qseq file is the first end
        int firstEndIndex = 0;
        // Which qseq file is the second end.
        int secondEndIndex = 1;
        // Which qseq file contains the barcode
        int barcodeIndex = -1;

        if (readConfiguration.isBarcoded()) {
            // Figure out which of up to 3 files contains the barcode.
            int cyclesSoFar = 1;
            for (barcodeIndex = 0; (barcodeCycle > cyclesSoFar) && (barcodeIndex < qseqLengths.length); ++barcodeIndex) {
                cyclesSoFar += qseqLengths[barcodeIndex];
            }
            if (barcodeCycle != cyclesSoFar) {
                throw new PicardException("Barcode cycle " + barcodeCycle + " does not fall on qseq boundary for barcode-aware qseqs.");
            }
            if (barcodeIndex <= secondEndIndex) {
                ++secondEndIndex;
            }
            if (barcodeIndex == firstEndIndex) {
                ++firstEndIndex;
            }
        }

        // Get the InclusiveRange objects in the order they have been read in this run.
        final ReadConfiguration.InclusiveRange[] ranges = new ReadConfiguration.InclusiveRange[numQseqs];
        ranges[firstEndIndex] = readConfiguration.getFirstRange();
        if (pairedEnd) {
            ranges[secondEndIndex] = readConfiguration.getSecondRange();
        }
        if (readConfiguration.isBarcoded()) {
            ranges[barcodeIndex] = readConfiguration.getBarcodeRange();
        }

        // Fill out the ranges according to where the barcode lies.
        int cyclesSoFar = 0;
        for (int i = 0; i < ranges.length; ++i) {
            ++cyclesSoFar;
            ranges[i].setStart(cyclesSoFar);
            cyclesSoFar += (qseqLengths[i] - 1);
            ranges[i].setEnd(cyclesSoFar);
        }

        if (readConfiguration.isBarcoded()) {
            if (readConfiguration.getBarcodeLength() != barcodeLength) {
                throw new PicardException("Barcode length from qseq file (" + readConfiguration.getBarcodeLength() +
                        ") != barcode length from command line (" + barcodeLength +").");
            }
        }
    }


    private void computeReadConfiguationFrom_1_1() {
        if (rawIntensityDirectory == null) {
            rawIntensityDirectory = basecallDirectory.getParentFile();
        }
        // For Bustard 1.1 we can look for the existence of 1 or 2 matrices off under
        // Matrix directory under the firecrest directory. The file names contain the cycle
        // number of the second cycle in each read!
        final Pattern pattern = Pattern.compile("s_" + lane + "_([0-9]+)_matrix.txt");
        final File[] matrices = new File(this.rawIntensityDirectory, "Matrix").listFiles(new FilenameFilter() {
            public boolean accept(final File parent, final String name) {
                return pattern.matcher(name).matches();
            }
        });

        if (pairedEnd == null) {
            if (matrices.length == 1) {
                this.pairedEnd = false;
            }
            else if (matrices.length == 2) {
                this.pairedEnd = true;
            }
            else {
                throw new IllegalStateException("Cannot determine if Bustard 1.1 output is PE or not.");
            }
        }

        int firstReadLength = -1;
        if (pairedEnd) {
            Matcher matcher = pattern.matcher(matrices[0].getName());
            if (!matcher.matches()) {
                throw new PicardException("That's unpossible");
            }
            final int c1 = Integer.parseInt(matcher.group(1));
            matcher = pattern.matcher(matrices[1].getName());
            if (!matcher.matches()) {
                throw new PicardException("That's unpossible");
            }
            final int c2 = Integer.parseInt(matcher.group(1));

            firstReadLength = Math.max(c1,c2) - Math.min(c1,c2);
        }
        readConfiguration.setPairedEnd(pairedEnd);
        final File seqFile = IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(basecallDirectory, "seq", lane)[0].file;
        final int totalReadLength = SeqParser.getReadLength(seqFile);
        readConfiguration.setFirstStart(1);
        if (this.pairedEnd) {
            readConfiguration.setFirstEnd(firstReadLength);
            readConfiguration.setSecondStart(firstReadLength + 1);
            readConfiguration.setSecondEnd(totalReadLength);
        } else {
            readConfiguration.setFirstEnd(totalReadLength);
        }
    }


    private void updateReadConfigurationForBarcode() {
        if (barcodeCycle == null) {
            readConfiguration.setBarcoded(false);
            return;
        }
        final BarcodeUtil.BarcodePosition barcodePosition =
                BarcodeUtil.findBarcodeEndAndStart(pairedEnd, readConfiguration.getFirstLength(),
                        readConfiguration.getSecondLength(), barcodeCycle, barcodeLength);
        readConfiguration.setBarcoded(true);
        readConfiguration.setBarcodeStart(barcodeCycle);
        readConfiguration.setBarcodeEnd(barcodeCycle + barcodeLength - 1);
        readConfiguration.setBarcodeRead(barcodePosition.barcodeIsInSecondRead ? EndType.SECOND: EndType.FIRST);
        if (barcodePosition.barcodeIsInSecondRead) {
            if (barcodeCycle == readConfiguration.getSecondStart()) {
                // Barcode is at start of second read
                readConfiguration.setSecondStart(barcodeCycle + barcodeLength);
            } else {
                // Barcode is at end of second read
                readConfiguration.setSecondEnd(barcodeCycle - 1);
            }
            if (readConfiguration.getSecondLength() == 0) {
                // PFS-113 Barcode occupies entire second end.  This is probably a run that was terminated prematurely.  
                readConfiguration.setPairedEnd(false);
                setPairedEnd(false);
            }
        } else if (barcodeCycle == readConfiguration.getFirstStart()) {
            // Barcode is at start of first read
            readConfiguration.setFirstStart(barcodeCycle + barcodeLength);
        } else {
            // Barcode is at end of first read
            readConfiguration.setFirstEnd(barcodeCycle - 1);
        }
    }

    /**
     * Add a data type to the list of desired data types specified in the ctor.
     *
     * @param dataType Data type that should be fetched when iterating over reads.
     */
    public void addDataType(final IlluminaDataType dataType) {
        dataTypes.add(dataType);
    }

    public void removeDataType(final IlluminaDataType dataType) {
        dataTypes.remove(dataType);
    }

    /**
     * In general it is not necessary to call this method, as the pairedness can be determined
     * from examining the basecall directory.
     * Must be called before computeReadConfiguration() is triggered, either directly or indirectly.
     */
    public void setPairedEnd(final boolean pairedEnd) {
        this.pairedEnd = pairedEnd;
    }

    public Boolean isPairedEnd() {
        return pairedEnd;
    }

    /**
     * In general it is not necessary to call this method, as the BaseCaller version can be determined
     * from examining the basecall directory.
     * Must be called before computeReadConfiguration() is triggered, either directly or indirectly.
     */
    public void setBaseCallerVersion(final BaseCallerVersion baseCallerVersion) {
        this.baseCallerVersion = baseCallerVersion;
    }

    public BaseCallerVersion getBaseCallerVersion() {
        return baseCallerVersion;
    }

    public ImageAnalyzerVersion getImageAnalyzerVersion() {
        return imageAnalyzerVersion;
    }

    /**
     * In general it is not necessary to call this method, as the ImageAnalyzer version can be determined
     * from examining the basecall directory.
     * Must be called before computeReadConfiguration() is triggered, either directly or indirectly.
     */
    public void setImageAnalyzerVersion(final ImageAnalyzerVersion imageAnalyzerVersion) {
        this.imageAnalyzerVersion = imageAnalyzerVersion;
    }

    /**
     * In general it is not necessary to call this method, as the raw intensity directory can be determined
     * from examining the basecall directory.
     * Must be called before computeReadConfiguration() is triggered, either directly or indirectly.
     */
    public void setRawIntensityDirectory(final File rawIntensityDirectory) {
        this.rawIntensityDirectory = rawIntensityDirectory;
    }

    private boolean isBarcodeAwareBaseCaller() {
        return (baseCallerVersion == BaseCallerVersion.Bustard_1_6) ||
                (baseCallerVersion == BaseCallerVersion.Bustard_1_5);
    }

    public void setAllowZeroLengthFirstEnd(final boolean allowZeroLengthFirstEnd) {
        this.allowZeroLengthFirstEnd = allowZeroLengthFirstEnd;
    }

    /**
     * Emulate the Groovy Elvis operator
     * @return v1 if not null, else v2
     */
    private <T> T elvisOperator(final T v1, final T v2) {
        if (v1 != null) {
            return v1;
        }
        return v2;
    }

    private void detectPipelineVersion() {
        final File solexaBuildVersion = new File(basecallDirectory, ".solexaBuildVersion");
        if (solexaBuildVersion.exists()) {
            final AsciiLineReader reader;
            try {
                reader = new AsciiLineReader(new FileInputStream(solexaBuildVersion));
            } catch (FileNotFoundException e) {
                throw new PicardException("Unexpected FileNotFound: " + solexaBuildVersion, e);
            }
            final String version = reader.readLine();
            log.info("solexaBuildVersion: ", version);
            reader.close();
            if (version.startsWith("1.6")) {
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_6);
                imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.rta);
            }
            else if (version.startsWith("1.5")) {
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_5);
                imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.rta);
            }
            else if (version.startsWith("1.4")) {
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_4);
                imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.rta);
            }
            else if (version.startsWith("1.3")) {
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_3);
                imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.Firecrest_1_3);
            } else if (version.startsWith("1.1")) {
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_1);
                imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.Firecrest_1_1);
            } else {
                log.info("Unrecognized solexaBuildVersion, using fallback to detect version. ", version);
            }
        } else {
            log.info(solexaBuildVersion, " not found, using fallback to detect version.");
        }
        if (new File(basecallDirectory.getParentFile(), "BaseCalls").exists()) {
            log.info("Found BaseCalls directory.");
            if (IlluminaFileUtil.endedIlluminaBasecallFilesExist(basecallDirectory, "qseq", lane, 3)) {
                // 3rd qseq file means barcode-aware base caller.
                log.info("Found barcode qseq file.");
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_5);
            } else {
                log.info("Did not find barcode qseq file.");
                baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_4);
            }
            imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.rta);
        } else if (IlluminaFileUtil.endedIlluminaBasecallFilesExist(basecallDirectory, "qseq", lane, 1)) {
            log.info("Found end-specific qseq file.");
            baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_3);
            imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.Firecrest_1_3);
        }
        else if (IlluminaFileUtil.nonEndedIlluminaBasecallFilesExist(basecallDirectory, "seq", lane)) {
            log.info("Found non-end-specific qseq file.");
            baseCallerVersion = elvisOperator(baseCallerVersion, BaseCallerVersion.Bustard_1_1);
            imageAnalyzerVersion = elvisOperator(imageAnalyzerVersion, ImageAnalyzerVersion.Firecrest_1_1);
        }
        else {
            throw new IllegalStateException("Cannot determine pipeline version for basecall directory: " + basecallDirectory);
        }
        log.info("BaseCallerVersion: " + baseCallerVersion);
        log.info("ImageAnalyzerVersion: " + imageAnalyzerVersion);
    }
}
