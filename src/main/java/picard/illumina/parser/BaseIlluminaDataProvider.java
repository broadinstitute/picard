package picard.illumina.parser;

import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.StreamSupport;

/**
 * Parse various formats and versions of Illumina Basecall files, and use them the to populate
 * ClusterData objects.
 */
public abstract class BaseIlluminaDataProvider implements Iterator<ClusterData>, Iterable<ClusterData>, AutoCloseable {

    public static final Pattern FILE_NAME_PATTERN = Pattern.compile("^s_\\d+_(\\d{1,5}).+");
    protected final int lane;
    /**
     * Calculated once, outputReadTypes describes the type of read data for each ReadData that will be found in output ClusterData objects
     */
    final ReadType[] outputReadTypes;
    /**
     * Number of reads in each ClusterData
     */
    final int numReads;

    public BaseIlluminaDataProvider(final int lane, final OutputMapping outputMapping) {
        numReads = outputMapping.numOutputReads();
        this.lane = lane;
        this.outputReadTypes = StreamSupport.stream(outputMapping.getOutputDescriptors().spliterator(), false)
                .map(rd -> rd.type).toArray(ReadType[]::new);
    }

    @Override
    public Iterator<ClusterData> iterator() {
        return this;
    }

    public abstract void close();

    /*
     * Methods for that transfer data from the IlluminaData objects to the current cluster
     */
    protected void addData(final ClusterData clusterData, final PositionalData posData) {
        clusterData.setX(posData.getXCoordinate());
        clusterData.setY(posData.getYCoordinate());
    }

    protected void addData(final ClusterData clusterData, final PfData pfData) {
        clusterData.setPf(pfData.isPf());
    }

    protected void addData(final ClusterData clusterData, final BarcodeData barcodeData) {
        clusterData.setMatchedBarcode(barcodeData.getBarcode());
    }

    protected void addReadData(final ClusterData clusterData, final int numReads, final BaseData baseData) {
        final byte[][] bases = baseData.getBases();
        for (int i = 0; i < numReads; i++) {
            clusterData.getRead(i).setBases(bases[i]);
        }
    }

    protected void addReadData(final ClusterData clusterData, final int numReads, final QualityData qualityData) {
        final byte[][] qualities = qualityData.getQualities();
        for (int i = 0; i < numReads; i++) {
            clusterData.getRead(i).setQualities(qualities[i]);
        }
    }

    protected void addReadData(final ClusterData clusterData, final int numReads, final CbclData cbclData) {
        final byte[][] bases = cbclData.getBases();
        for (int i = 0; i < numReads; i++) {
            clusterData.getRead(i).setBases(bases[i]);
        }
        final byte[][] qualities = cbclData.getQualities();
        for (int i = 0; i < numReads; i++) {
            clusterData.getRead(i).setQualities(qualities[i]);
        }
        clusterData.setPf(cbclData.isPf());
        clusterData.setX(cbclData.getPositionInfo().xQseqCoord);
        clusterData.setY(cbclData.getPositionInfo().yQseqCoord);
    }

    abstract void seekToTile(Integer seekAfterFirstRead);

    public static Integer fileToTile(final String fileName) {
        final Matcher matcher = FILE_NAME_PATTERN.matcher(fileName);
        if (!matcher.matches()) {
            return null;
        }
        return Integer.parseInt(matcher.group(1));
    }
}
