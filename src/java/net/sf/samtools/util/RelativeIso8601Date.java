package net.sf.samtools.util;

import java.util.Date;

/**
 * Like {@link Iso8601Date}, but also comes in a "lazy now" flavor.
 * 
 * When "lazy now" mode is enabled, this instance's date value is undefined until the first time it is queried, at which time it is set to 
 * {@link System#currentTimeMillis()}.  This value is returned on subsequent queries, so it is consistent.
 * 
 * The "lazy state" is conveyed via {@link #toString()}.  A "lazy now" instance will answer {@link #toString()} with
 * {@link #LAZY_NOW_LABEL} if the time has not yet been queried/set, or a {@link Iso8601Date}-formatted date of the query time if it
 * has been queried.  This characteristic is useful for serialization and persistence purposes.
 * 
 * Consumers can create "lazy now" instances via the {@link #generateLazyNowInstance()} factory method or by passing {@link #LAZY_NOW_LABEL} to 
 * {@link #RelativeIso8601Date(String)}.
 * 
 * @author mccowan
 */
public class RelativeIso8601Date extends Iso8601Date {

    public static final String LAZY_NOW_LABEL = "NOW";
    
    /** Flag that indicates this instance is lazy and has not yet been queried (and so its value should be updated at the next query). */
    private boolean doSetTimeNextQuery;

    /** Returns a "lazy now" instance. */
    public static RelativeIso8601Date generateLazyNowInstance() {
        return new RelativeIso8601Date(LAZY_NOW_LABEL);
    }
    
    public RelativeIso8601Date(final Date date) {
        super(date);
        doSetTimeNextQuery = false;
    }

    public RelativeIso8601Date(final String dateStr) {
        /** 
         * We must pass a date parsable {@link Iso8601Date#Iso8601Date(String)}; we will never actually read the passed value, so it doesn't
         * matter what it is.
         */
        super(LAZY_NOW_LABEL.equals(dateStr) ? new Iso8601Date(new Date()).toString() : dateStr);
        doSetTimeNextQuery = LAZY_NOW_LABEL.equals(dateStr);
    }

    /** Updates the time stored by this instance if it's a "lazy now" instance and has never been stored. */
    private synchronized void conditionallyUpdateTime() {
        if (!doSetTimeNextQuery) {
            super.setTime(System.currentTimeMillis());
            doSetTimeNextQuery = false;
        }
    }

    /**
     * Returns a {@link String} representation of this date.
     * @return An {@link Iso8601Date}-formatted string, or the value of {@link #LAZY_NOW_LABEL} if this is a "lazy now" instance.
     */
    @Override
    public String toString() {
        return doSetTimeNextQuery ? LAZY_NOW_LABEL : super.toString();
    }
    
    @Override
    public long getTime() {
        conditionallyUpdateTime();
        return super.getTime();
    }
    
    @Override
    public boolean after(final Date when) {
        conditionallyUpdateTime();
        return super.after(when);
    }

    @Override
    public boolean before(final Date when) {
        conditionallyUpdateTime();
        return super.before(when);
    }

    @Override
    public Object clone() {
        conditionallyUpdateTime();
        return super.clone();
    }

    @Override
    public int compareTo(final Date anotherDate) {
        conditionallyUpdateTime();
        return super.compareTo(anotherDate);
    }

    @Override
    public boolean equals(final Object obj) {
        conditionallyUpdateTime();
        return super.equals(obj);
    }

    @Override
    public int getDate() {
        conditionallyUpdateTime();
        return super.getDate();
    }

    @Override
    public int getDay() {
        conditionallyUpdateTime();
        return super.getDay();
    }

    @Override
    public int getHours() {
        conditionallyUpdateTime();
        return super.getHours();
    }

    @Override
    public int getMinutes() {
        conditionallyUpdateTime();
        return super.getMinutes();
    }

    @Override
    public int getMonth() {
        conditionallyUpdateTime();
        return super.getMonth();
    }

    @Override
    public int getSeconds() {
        conditionallyUpdateTime();
        return super.getSeconds();
    }

    @Override
    public int getTimezoneOffset() {
        conditionallyUpdateTime();
        return super.getTimezoneOffset();
    }

    @Override
    public int getYear() {
        conditionallyUpdateTime();
        return super.getYear();
    }

    @Override
    public int hashCode() {
        conditionallyUpdateTime();
        return super.hashCode();
    }

    @Override
    public void setDate(final int date) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setHours(final int hours) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setMinutes(final int minutes) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setMonth(final int month) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setSeconds(final int seconds) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setTime(final long time) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setYear(final int year) {
        throw new UnsupportedOperationException();
    }
}
