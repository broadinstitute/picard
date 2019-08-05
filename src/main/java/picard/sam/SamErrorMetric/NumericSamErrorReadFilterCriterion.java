package picard.sam.SamErrorMetric;

public class NumericSamErrorReadFilterCriterion extends SamErrorReadFilterCriterion<Integer> {
    public NumericSamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final Integer value) {
        super(comparator, value);
    }

    @Override
    public void checkCriterion(final Integer stratus) {
        switch (comparator) {
            case Equal: if (stratus.equals(value)) { setSatisifed(); } break;
            case NotEqual: if (!stratus.equals(value)) { setSatisifed(); } break;
            case Smaller: if (stratus < value) { setSatisifed(); } break;
            case SmallerOrEqual: if (stratus <= value) { setSatisifed(); } break;
            case Greater: if (stratus > value) { setSatisifed(); } break;
            case GreaterOrEqual: if (stratus >= value) { setSatisifed(); } break;
            default: throw new IllegalArgumentException("Invalid operator for this criterion type.");
        }
    }
}