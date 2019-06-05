package picard.sam.SamErrorMetric;


import htsjdk.samtools.util.Log;

public class BooleanSamErrorReadFilterCriterion extends SamErrorReadFilterCriterion<Boolean> {
    public BooleanSamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final Boolean value) {
        super(comparator, value);
    }

    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    @Override
    public void checkCriterion(final Object stratus) {
        Boolean boolValue = (Boolean)stratus;
        switch (comparator) {
            case Equal: if(boolValue == value) { setSatisifed(); } break;
            case NotEqual: if(boolValue != value) { setSatisifed(); } break;
            default: throw new IllegalArgumentException("Invalid operator for this criterion type.");
        }
    }
}