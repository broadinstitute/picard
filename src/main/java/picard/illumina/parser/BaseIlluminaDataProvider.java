package picard.illumina.parser;

import java.util.Iterator;

public abstract class BaseIlluminaDataProvider implements Iterator<ClusterData>, Iterable<ClusterData> {

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
        this.outputReadTypes = new ReadType[numReads];
        final int[] i = {0};
        outputMapping.getOutputDescriptors().forEach(rd -> outputReadTypes[i[0]++] = rd.type);
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

    abstract void seekToTile(int seekAfterFirstRead);
}
