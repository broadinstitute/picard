package picard.sam.SamErrorMetric;

/**
 * Class holding a criterion to be used in the SamErrorReadFilter class
 * @param <T>
 */
public abstract class SamErrorReadFilterCriterion<T extends Comparable<T>> {
    public final SamErrorReadFilter.Comparator comparator;
    public final T value;
    private boolean satisifed = false;

    /**
     * Creates a criterion object
     * @param comparator Comparison operation to perform
     * @param value Value to compare with to satisfy the criterion
     */
    public SamErrorReadFilterCriterion(final SamErrorReadFilter.Comparator comparator, final T value) {
        this.comparator = comparator;
        this.value = value;
    }

    /**
     * Method to be called in order to check whether this criterion is satisfied
     * @param stratus Object to check whether it satisfies the criterion
     */
    public abstract void checkCriterion(T stratus);

    /**
     * Returns whether or not the criterion is satisfied
     */
    public boolean isSatisifed() {
        return satisifed;
    }

    protected void setSatisifed() {
        satisifed = true;
    }

    /**
     * Resets this criterion. This needs to be called whenever a new RecordAndOffset is considered for filtering.
     */
    public void reset() {
        satisifed = false;
    }
}
