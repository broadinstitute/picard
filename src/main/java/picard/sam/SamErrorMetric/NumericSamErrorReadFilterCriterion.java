package picard.sam.SamErrorMetric;


import htsjdk.samtools.util.Log;

public class NumericSamErrorReadFilterCriterion extends SamErrorReadFilterCriterion<Integer> {
    public NumericSamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final Integer value) {
        super(comparator, value);
    }

    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    @Override
    public void checkCriterion(final Object stratus) {
        if(stratus == null)
            return;
        Integer intValue = (Integer)stratus;
        switch (comparator) {
            case Equal: if (intValue.equals(value)) { setSatisifed(); } break;
            case NotEqual: if (!intValue.equals(value)) { setSatisifed(); } break;
            case Smaller: if (intValue < value) { setSatisifed(); } break;
            case SmallerOrEqual: if (intValue <= value) { setSatisifed(); } break;
            case Greater: if (intValue > value) { setSatisifed(); } break;
            case GreaterOrEqual: if (intValue >= value) { setSatisifed(); } break;
            default: throw new IllegalArgumentException("Invalid operator for this criterion type.");
        }
    }
}