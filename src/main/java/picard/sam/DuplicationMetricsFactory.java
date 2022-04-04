package picard.sam;

public class DuplicationMetricsFactory {

    // create a DuplicationMetrics for a specific read group
    public static DuplicationMetrics createMetrics(final boolean flowMetrics) {

        // create based on the presence of flow order
        if ( !flowMetrics ) {
            return new DuplicationMetrics();
        } else {
            return new FlowBasedDuplicationMetrics();
        }
    }

    public static DuplicationMetrics createMetrics() {
        return new DuplicationMetrics();
    }
}
