package picard.illumina.parser;

import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.BarcodeFileReader;
import picard.illumina.parser.readers.CbclReader;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parse cbcls Illumina Basecall files, and use them the to populate
 * ClusterData objects.
 */
public class NewIlluminaDataProvider extends BaseIlluminaDataProvider {
    private final CbclReader reader;
    private BarcodeFileReader barcodeReader = null;

    /**
     * Construct a NewIlluminaDataProvider to create a ClusterData iterator over all clusters for a given tile.
     *
     * @param cbcls              A list of cbcls to use when creating this data provider.
     * @param filterFiles        A list of the pf filter files to use when creating this data provider.
     */
    NewIlluminaDataProvider(final List<File> cbcls, final List<AbstractIlluminaPositionFileReader.PositionInfo> locs,
                            final File[] filterFiles, final int lane, final int tileNum,
                            final OutputMapping outputMapping, final File barcodeFile) {
        super(lane, outputMapping);

        Map<Integer, File> filterFileMap = new HashMap<>();
        for (File filterFile : filterFiles) {
            filterFileMap.put(fileToTile(filterFile.getName()), filterFile);
        }
        this.reader = new CbclReader(cbcls, filterFileMap, outputMapping.getOutputReadLengths(), tileNum, locs, outputMapping.getOutputCycles(), false);
        if (barcodeFile != null) {
            this.barcodeReader = new BarcodeFileReader(barcodeFile);
        }
    }

    @Override
    public void close() {
        reader.clear();
        reader.close();
    }

    @Override
    public void seekToTile(int seekAfterFirstRead) {

    }

    @Override
    public boolean hasNext() {
        return reader.hasNext();
    }

    @Override
    public ClusterData next() {
        CbclData cbclData = reader.next();

        if (cbclData == null) return null;

        final ClusterData cluster = new ClusterData(outputReadTypes);
        cluster.setLane(lane);
        cluster.setTile(cbclData.getTile());
        if (barcodeReader != null) {
            cluster.setMatchedBarcode(barcodeReader.next());
        }

        addReadData(cluster, numReads, cbclData);

        return cluster;
    }

    public static Integer fileToTile(final String fileName) {
        final Matcher matcher = Pattern.compile("^s_\\d+_(\\d{1,5}).+").matcher(fileName);
        if (!matcher.matches()) {
            return null;
        }
        return Integer.parseInt(matcher.group(1));
    }
}
