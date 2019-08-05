package picard.sam.SamErrorMetric;

public class BooleanSamErrorReadFilterCriterion extends SamErrorReadFilterCriterion<Boolean> {
    public BooleanSamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final Boolean value) {
        super(comparator, value);
    }

    @Override
    public void checkCriterion(final Boolean stratus) {
        switch (comparator) {
            case Equal: if(stratus == value) { setSatisifed(); } break;
            case NotEqual: if(stratus != value) { setSatisifed(); } break;
            default: throw new IllegalArgumentException("Invalid operator for this criterion type.");
        }
    }
}