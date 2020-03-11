package picard.fastq;

import htsjdk.samtools.util.StringUtil;
import picard.illumina.parser.ClusterData;

/**
 * A read name encoder conforming to the standard described by Illumina Casava 1.8.
 * 
 * @see <a href="http://biowulf.nih.gov/apps/CASAVA1_8_Changes.pdf">Casava 1.8 update</a>
 * @author mccowan
 */
public class Casava18ReadNameEncoder implements ReadNameEncoder {
    final static int CONTROL_FIELD_VALUE = 0;
    final String runId, instrumentName, flowcellId;
    
    static enum IsFilteredLabel {
        Y, N;
        static IsFilteredLabel get(final boolean passesFilter) {
            return passesFilter ? N : Y;
        }
    }
    
    public Casava18ReadNameEncoder(final String instrumentName, final String runId, final String flowcellId) {
        this.runId = runId;
        this.instrumentName = instrumentName;
        this.flowcellId = flowcellId;
    }

    @Override
    public String generateReadName(final ClusterData cluster, final Integer pairNumber) {
        return new StringBuilder().append(instrumentName).append(":")
                .append(runId).append(":")
                .append(flowcellId).append(":")
                .append(cluster.getLane()).append(":")
                .append(cluster.getTile()).append(":")
                .append(cluster.getX()).append(":")
                .append(cluster.getY()).append(" ")
                .append(StringUtil.asEmptyIfNull(pairNumber)).append(":")
                .append(IsFilteredLabel.get(cluster.isPf())).append(":")
                .append(CONTROL_FIELD_VALUE).append(":")
                .append(StringUtil.asEmptyIfNull(cluster.getMatchedBarcode())).toString();
    }
}
