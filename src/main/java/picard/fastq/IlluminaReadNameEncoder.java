package picard.fastq;

import picard.illumina.parser.ClusterData;

/**
 * A read name encoder following the encoding initially produced by picard fastq writers.
 * 
 * <a href="http://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers">Illumina sequence identifiers</a>
 * almost describes the format used here, except instead of an instrument name, we write the run barcode.
 *
 * @author mccowan
 */
public class IlluminaReadNameEncoder implements ReadNameEncoder {
    final String runBarcode;
    public IlluminaReadNameEncoder(final String runBarcode) {
        this.runBarcode = runBarcode;
    }
    
    @Override
    public String generateReadName(final ClusterData cluster, final Integer pairNumber) {
        return generateShortName(cluster) + generatePairNumberSuffix(pairNumber);
    }

    @Override
    public String generateShortName(ClusterData cluster) {
        return runBarcode + ":" + cluster.getLane() + ":" + cluster.getTile() + ":" + cluster.getX() + ":" + cluster.getY();
    }

    private static String generatePairNumberSuffix(final Integer pairNumber) {
        if (pairNumber == null)
            return ""; 
        else
            return "/" + pairNumber;
    }
}
