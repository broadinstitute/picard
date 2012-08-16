package net.sf.picard.util;

import net.sf.samtools.SAMRecord;

import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Little progress logging class to facilitate consistent output of useful information when progressing
 * through a stream of SAM records.
 *
 * @author Tim Fennell
 */
public class ProgressLogger {
    private final Log log;
    private final int n;
    private final String verb;
    private final long startTime = System.currentTimeMillis();
    
    private final NumberFormat fmt = new DecimalFormat("#,###");
    private final NumberFormat timeFmt = new DecimalFormat("00");
    
    private long processed = 0;
    private long lastStartTime = startTime;

    /**
     * Construct a progress logger.
     * @param log the Log object to write outputs to
     * @param n the frequency with which to output (i.e. every N records)
     * @param verb the verb to log, e.g. "Processed, Read, Written".
     */
    public ProgressLogger(final Log log, final int n, final String verb) {
        this.log = log;
        this.n = n;
        this.verb = verb;
    }

    /**
     * Construct a progress logger with the desired log and frequency and the verb "Processed".
     * @param log the Log object to write outputs to
     * @param n the frequency with which to output (i.e. every N records)
     */
    public ProgressLogger(final Log log, final int n) { this(log, n, "Processed"); }

    /**
     * Construct a progress logger with the desired log, the verb "Processed" and a period of 1m records.
     * @param log the Log object to write outputs to
     */
    public ProgressLogger(final Log log) { this(log, 1000000); }

    /**
     * Records that a given record has been processed and triggers logging if necessary.
     * @return boolean true if logging was triggered, false otherwise
     */
    public synchronized boolean record(final SAMRecord rec) {
        if (++this.processed % this.n == 0) {
            final long now = System.currentTimeMillis();
            final long lastPeriodSeconds = (now - this.lastStartTime) / 1000;
            this.lastStartTime = now;

            final long seconds = (System.currentTimeMillis() - startTime) / 1000;
            final String elapsed   = formatElapseTime(seconds);
            final String period    = pad(fmt.format(lastPeriodSeconds), 4);
            final String processed = pad(fmt.format(this.processed), 13);
            
            final String readInfo;
            if (rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                readInfo = "*/*";
            } else {
                readInfo = rec.getReferenceName() + ":" + fmt.format(rec.getAlignmentStart());
            }
            
            log.info(this.verb, " ", processed, " records.  Elapsed time: ", elapsed, "s.  Time for last ", fmt.format(this.n),
                     ": ", period, "s.  Last read position: ", readInfo);
            return true;
        }
        else {
            return false;
        }
    }
    
    /** Records multiple SAMRecords and triggers logging if necessary. */
    public boolean record(final SAMRecord... recs) {
        boolean triggered = false;
        for (final SAMRecord rec : recs) triggered = record(rec) || triggered;
        return triggered;
    }
    
    /** Returns the count of records processed. */
    public long getCount() { return this.processed; }

    /** Returns the number of seconds since progress tracking began. */
    public long getElapsedSeconds() { return (System.currentTimeMillis() - this.startTime) / 1000; }
    
    /** Left pads a string until it is at least the given length. */
    private String pad (String in, final int length) {
        while (in.length() < length) {
            in = " " + in;
        }
        
        return in;
    }
    
    /** Formats a number of seconds into hours:minutes:seconds. */
    private String formatElapseTime(final long seconds) {
        final long s = seconds % 60;
        final long allMinutes = seconds / 60;
        final long m = allMinutes % 60;
        final long h = allMinutes / 60;

        return timeFmt.format(h) + ":" + timeFmt.format(m) + ":" + timeFmt.format(s);
    }
}
