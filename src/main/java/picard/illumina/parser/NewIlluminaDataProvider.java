package picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.BarcodeFileReader;
import picard.illumina.parser.readers.CbclReader;
import picard.illumina.parser.readers.LocsFileReader;

import java.io.File;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static picard.illumina.BasecallsConverter.TILE_NUMBER_COMPARATOR;
import static picard.illumina.BasecallsConverter.getTiledFiles;

/**
 * Parse cbcls Illumina Basecall files, and use them the to populate
 * ClusterData objects.
 */
class NewIlluminaDataProvider extends BaseIlluminaDataProvider {
    private final Map<Integer, CbclReader>  tileReaders = new HashMap<>();
    private final Map<Integer, BarcodeFileReader> barcodeFileMap;
    /** The current tile number */
    protected Integer currentTile;
    private final TreeSet<Integer> tileOrder = new TreeSet<>();

    /**
     * Construct a NewIlluminaDataProvider to create a ClusterData iterator over all clusters for a given tile.
     *
     * @param outputMapping     Mapping of reads types to be output.
     * @param basecallDirectory The baseCalls directory of a complete Illumina directory.
     * @param lane              The lane that to provide data for.
     * @param requestedTiles    The list of tiles that data is requested for.
     */
    NewIlluminaDataProvider(final OutputMapping outputMapping,
                            final File basecallDirectory, final int lane, List<Integer> requestedTiles) {
        super(lane, outputMapping);
        requestedTiles.stream().sorted(TILE_NUMBER_COMPARATOR).forEach(tileOrder::add);
        currentTile = tileOrder.first();
        final File laneDir = new File(basecallDirectory, IlluminaFileUtil.longLaneStr(lane));

        final File[] cycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

        //CBCLs
        final List<File> cbcls = Arrays.stream(cycleDirs)
                .flatMap(cycleDir -> Arrays.stream(IOUtil.getFilesMatchingRegexp(cycleDir,
                        "^" + IlluminaFileUtil.longLaneStr(lane) + "_(\\d{1,5}).cbcl$"))).collect(Collectors.toList());

        if (cbcls.size() == 0) {
            throw new PicardException("No CBCL files found.");
        }

        IOUtil.assertFilesAreReadable(cbcls);
        //locs
        final List<AbstractIlluminaPositionFileReader.PositionInfo> locs = new ArrayList<>();
        final File locsFile = new File(basecallDirectory.getParentFile(), AbstractIlluminaPositionFileReader.S_LOCS_FILE);
        IOUtil.assertFileIsReadable(locsFile);
        try (LocsFileReader locsFileReader = new LocsFileReader(locsFile)) {
            while (locsFileReader.hasNext()) {
                locs.add(locsFileReader.next());
            }
        }

        //barcodes
        final Pattern barcodeRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeBarcodeRegex(lane)));
        final File[] barcodeFiles = getTiledFiles(basecallDirectory, barcodeRegex);
        this.barcodeFileMap = new HashMap<>();
        for (File barcodeFile : barcodeFiles) {
            barcodeFileMap.put(fileToTile(barcodeFile.getName()), new BarcodeFileReader(barcodeFile));
        }

        //filter
        final Pattern laneTileRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeLaneTileRegex(".filter", lane)));
        final File[] filterFiles = getTiledFiles(laneDir, laneTileRegex);

        IOUtil.assertFilesAreReadable(Arrays.asList(filterFiles));

        Map<Integer, File> filterFileMap = new HashMap<>();
        for (File filterFile : filterFiles) {
            filterFileMap.put(fileToTile(filterFile.getName()), filterFile);
        }
        for(Integer tile: requestedTiles) {
            this.tileReaders.put(tile, new CbclReader(cbcls, filterFileMap, outputMapping.getOutputReadLengths(), tile, locs, outputMapping.getOutputCycles(), false));
        }
    }

    @Override
    public void close() {
        for(CbclReader reader: this.tileReaders.values()) {
            reader.clear();
            reader.close();
        }
    }


    /*
     * Clear the current set of cycleFileParsers and replace them with the ones for the tile indicated by oneBasedTileNumber
     *
     * @param tile requested tile
     */
    @Override
    public void seekToTile(final Integer tile) {
        currentTile = tile;
    }

    @Override
    public boolean hasNext() {
        return tileReaders.get(currentTile).hasNext() || currentTile < tileOrder.last();
    }

    @Override
    public ClusterData next() {

        if (!tileReaders.get(currentTile).hasNext()) {
            seekToTile(tileOrder.higher(currentTile));
        }
        CbclData cbclData = tileReaders.get(currentTile).next();

        if (cbclData == null) return null;

        final ClusterData cluster = new ClusterData(outputReadTypes);
        cluster.setLane(lane);
        cluster.setTile(cbclData.getTile());
        BarcodeFileReader barcodeReader = barcodeFileMap.get(cbclData.getTile());
        if (barcodeReader != null) {
            cluster.setMatchedBarcode(barcodeReader.next());
        }

        addReadData(cluster, numReads, cbclData);

        return cluster;
    }

}
