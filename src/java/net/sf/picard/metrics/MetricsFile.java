/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package net.sf.picard.metrics;

import net.sf.picard.PicardException;
import net.sf.picard.util.FormatUtil;
import net.sf.picard.util.Histogram;
import net.sf.samtools.util.StringUtil;

import java.io.*;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

/**
 * Contains a set of metrics that can be written to a file and parsed back
 * again. The set of metrics is composed of zero or more instances of a class,
 * BEAN, that extends {@link MetricBase} (all instances must be of the same type)
 * and may optionally include one or more histograms that share the same key set.
 *
 * @author Tim Fennell
 */
public class MetricsFile<BEAN extends MetricBase, HKEY extends Comparable> {
    public static final String MAJOR_HEADER_PREFIX = "## ";
    public static final String MINOR_HEADER_PREFIX = "# ";
    public static final String SEPARATOR = "\t";
    public static final String HISTO_HEADER = "## HISTOGRAM\t";
    public static final String METRIC_HEADER = "## METRICS CLASS\t";

    private final List<Header> headers = new ArrayList<Header>();
    private final List<BEAN> metrics = new ArrayList<BEAN>();
    private final List<Histogram<HKEY>> histograms = new ArrayList<Histogram<HKEY>>();

    /** Adds a header to the collection of metrics. */
    public void addHeader(Header h) { this.headers.add(h); }

    /** Returns the list of headers. */
    public List<Header> getHeaders() { return Collections.unmodifiableList(this.headers); }

    /** Adds a bean to the collection of metrics. */
    public void addMetric(BEAN bean) { this.metrics.add(bean); }

    /** Returns the list of headers. */
    public List<BEAN> getMetrics() { return Collections.unmodifiableList(this.metrics); }

    /** Returns the histogram contained in the metrics file if any. */
    public Histogram<HKEY> getHistogram() {
        if (histograms.size() > 0) return this.histograms.get(0);
        else return null;
    }

    /** Sets the histogram contained in the metrics file. */
    public void setHistogram(Histogram<HKEY> histogram) {
        if (this.histograms.isEmpty()) this.histograms.add(histogram);
        else this.histograms.set(0, histogram);
    }

    /** Adds a histogram to the list of histograms in the metrics file. */
    public void addHistogram(Histogram<HKEY> histogram) {
        this.histograms.add(histogram);
    }

    /** Returns the list of headers with the specified type. */
    public List<Header> getHeaders(Class<? extends Header> type) {
        List<Header> tmp = new ArrayList<Header>();
        for (Header h : this.headers) {
            if (h.getClass().equals(type)) {
                tmp.add(h);
            }
        }

        return tmp;
    }

    /**
     * Writes out the metrics file to the supplied file. The file is written out
     * headers first, metrics second and histogram third.
     *
     * @param f a File into which to write the metrics
     */
    public void write(File f) {
        FileWriter w = null;
        try {
            w = new FileWriter(f);
            write(w);
        }
        catch (IOException ioe) {
            throw new PicardException("Could not write metrics to file: " + f.getAbsolutePath(), ioe);
        }
        finally {
            if (w != null) {
                try {
                    w.close();
                } catch (IOException e) {
                }
            }
        }
    }

    /**
     * Writes out the metrics file to the supplied writer. The file is written out
     * headers first, metrics second and histogram third.
     *
     * @param w a Writer into which to write the metrics
     */
    public void write(Writer w) {
        try {
            FormatUtil formatter = new FormatUtil();
            BufferedWriter out = new BufferedWriter(w);
            printHeaders(out);
            out.newLine();

            printBeanMetrics(out, formatter);
            out.newLine();

            printHistogram(out, formatter);
            out.newLine();
            out.flush();
        }
        catch (IOException ioe) {
            throw new PicardException("Could not write metrics file.", ioe);
        }
    }

    /** Prints the headers into the provided PrintWriter. */
    private void printHeaders(BufferedWriter out) throws IOException {
        for (Header h : this.headers) {
            out.append(MAJOR_HEADER_PREFIX);
            out.append(h.getClass().getName());
            out.newLine();
            out.append(MINOR_HEADER_PREFIX);
            out.append(h.toString());
            out.newLine();
        }
    }

    /** Prints each of the metrics entries into the provided PrintWriter. */
    private void printBeanMetrics(BufferedWriter out, FormatUtil formatter) throws IOException {
        if (this.metrics.isEmpty()) {
            return;
        }

        // Write out a header row with the type of the metric class
        out.append(METRIC_HEADER + getBeanType().getName());
        out.newLine();

        // Write out the column headers
        Field[] fields = getBeanType().getFields();
        final int fieldCount = fields.length;

        for (int i=0; i<fieldCount; ++i) {
            out.append(fields[i].getName());
            if (i < fieldCount - 1) {
                out.append(MetricsFile.SEPARATOR);
            }
            else {
                out.newLine();
            }
        }

        // Write out each of the data rows
        for (BEAN bean : this.metrics) {
            for (int i=0; i<fieldCount; ++i) {
                try {
                    Object value = fields[i].get(bean);
                    out.append(StringUtil.assertCharactersNotInString(formatter.format(value), '\t', '\n'));

                    if (i < fieldCount - 1) {
                        out.append(MetricsFile.SEPARATOR);
                    }
                    else {
                        out.newLine();
                    }
                }
                catch (IllegalAccessException iae) {
                    throw new PicardException("Could not read property " + fields[i].getName()
                            + " from class of type " + bean.getClass());
                }
            }
        }

        out.flush();
    }

    /** Prints the histogram if one is present. */
    private void printHistogram(BufferedWriter out, FormatUtil formatter) throws IOException {
        if (this.histograms.isEmpty()) {
            return;
        }

        // Build a combined key set
        java.util.Set<HKEY> keys = new TreeSet<HKEY>();
        for (Histogram<HKEY> histo : histograms) {
            keys.addAll(histo.keySet());
        }

        // Add a header for the histogram key type
        out.append(HISTO_HEADER + this.histograms.get(0).keySet().iterator().next().getClass().getName());
        out.newLine();

        // Output a header row
        out.append(StringUtil.assertCharactersNotInString(this.histograms.get(0).getBinLabel(), '\t', '\n'));
        for (Histogram<HKEY> histo : this.histograms) {
            out.append(SEPARATOR);
            out.append(StringUtil.assertCharactersNotInString(histo.getValueLabel(), '\t', '\n'));
        }
        out.newLine();

        for (HKEY key : keys) {
            out.append(key.toString());

            for (Histogram<HKEY> histo : this.histograms) {
                Histogram<HKEY>.Bin bin = histo.get(key);
                final double value = (bin == null ? 0 : bin.getValue());

                out.append(SEPARATOR);
                out.append(formatter.format(value));
            }

            out.newLine();
        }
    }

    /** Gets the type of the metrics bean being used. */
    private Class<?> getBeanType() {
        if (this.metrics == null || this.metrics.isEmpty()) {
            return null;
        } else {
            return this.metrics.get(0).getClass();
        }
    }

    /** Reads the Metrics in from the given reader. */
    public void read(Reader r) {
        BufferedReader in = new BufferedReader(r);
        FormatUtil formatter = new FormatUtil();
        String line = null;

        try {
            // First read the headers
            Header header = null;
            boolean inHeader = true;
            while ((line = in.readLine()) != null && inHeader) {
                line = line.trim();
                // A blank line signals the end of the headers, otherwise parse out
                // the header types and values and build the headers.
                if ("".equals(line)) {
                    inHeader = false;
                }
                else if (line.startsWith(MAJOR_HEADER_PREFIX)) {
                    if (header != null) {
                        throw new IllegalStateException("Consecutive header class lines encountered.");
                    }
                    
                    String className = line.substring(MAJOR_HEADER_PREFIX.length()).trim();
                    try {
                        header = (Header) loadClass(className).newInstance();
                    }
                    catch (Exception e) {
                        throw new PicardException("Error load and/or instantiating an instance of " + className, e);
                    }
                }
                else if (line.startsWith(MINOR_HEADER_PREFIX)) {
                    if (header == null) {
                        throw new IllegalStateException("Header class must precede header value:" + line);
                    }
                    header.parse(line.substring(MINOR_HEADER_PREFIX.length()));
                    this.headers.add(header);
                    header = null;
                }
                else {
                    throw new PicardException("Illegal state. Found following string in metrics file header: " + line);
                }
            }

            if (line == null) {
                throw new PicardException("No lines in metrics file after header.");
            }
            // Then read the metrics if there are any
            while (!line.startsWith(MAJOR_HEADER_PREFIX)) {
                line = in.readLine().trim();
            }
            if (line.startsWith(METRIC_HEADER)) {
                // Get the metric class from the header
                String className = line.split(SEPARATOR)[1];
                Class<?> type = null;
                try {
                    type = loadClass(className);
                }
                catch (ClassNotFoundException cnfe) {
                    throw new PicardException("Could not locate class with name " + className, cnfe);
                }

                // Read the next line with the column headers
                String[] fieldNames = in.readLine().split(SEPARATOR);
                Field[] fields = new Field[fieldNames.length];
                for (int i=0; i<fieldNames.length; ++i) {
                    try {
                        fields[i] = type.getField(fieldNames[i]);
                    }
                    catch (Exception e) {
                        throw new PicardException("Could not get field with name " + fieldNames[i] +
                            " from class " + type.getName());
                    }
                }

                // Now read the values
                while ((line = in.readLine()) != null) {
                    line = line.trim();
                    if ("".equals(line)) {
                        break;
                    }
                    else {
                        String[] values = line.split(SEPARATOR);
                        BEAN bean = null;

                        try { bean = (BEAN) type.newInstance(); }
                        catch (Exception e) { throw new PicardException("Error instantiating a " + type.getName(), e); }


                        for (int i=0; i<fields.length; ++i) {
                            Object value = null;
                            if (values[i] != null && values[i].length() > 0) {
                                value = formatter.parseObject(values[i], fields[i].getType());
                            }

                            try { fields[i].set(bean, value); }
                            catch (Exception e) {
                                throw new PicardException("Error setting field " + fields[i].getName() +
                                        " on class of type " + type.getName(), e);
                            }
                        }

                        this.metrics.add(bean);
                    }
                }
            }

            // Then read the histograms if any are present
            while (line != null && !line.startsWith(MAJOR_HEADER_PREFIX)) {
                line = in.readLine();
            }
            if (line != null && line.startsWith(HISTO_HEADER)) {
                // Get the key type of the histogram
                String keyClassName = line.split(SEPARATOR)[1].trim();
                Class<?> keyClass = null;

                try { keyClass = loadClass(keyClassName); }
                catch (ClassNotFoundException cnfe) { throw new PicardException("Could not load class with name " + keyClassName); }

                // Read the next line with the bin and value labels
                String[] labels = in.readLine().split(SEPARATOR);
                for (int i=1; i<labels.length; ++i) {
                    this.histograms.add(new Histogram<HKEY>(labels[0], labels[i]));
                }

                // Read the entries in the histograms
                while ((line = in.readLine()) != null && !"".equals(line)) {
                    String[] fields = line.trim().split(SEPARATOR);
                    HKEY key = (HKEY) formatter.parseObject(fields[0], keyClass);

                    for (int i=1; i<fields.length; ++i) {
                        double value = formatter.parseDouble(fields[i]);
                        this.histograms.get(i-1).increment(key, value);
                    }
                }
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Could not read metrics from reader.", ioe);
        }
    }

    /** Attempts to load a class, taking into account that some classes have "migrated" from the broad to sf. */
    private Class<?> loadClass(String className) throws ClassNotFoundException {
        try {
            return Class.forName(className);
        }
        catch (ClassNotFoundException cnfe) {
            if (className.startsWith("edu.mit.broad.picard")) {
                return loadClass(className.replace("edu.mit.broad.picard", "net.sf.picard"));
            }
            else {
                throw cnfe;
            }
        }
    }

    /** Checks that the headers, metrics and histogram are all equal. */
    @Override
    public boolean equals(Object o) {
        if (o == null) {
            return false;
        }
        if (getClass() != o.getClass()) {
            return false;
        }
        MetricsFile that = (MetricsFile) o;

        if (!this.headers.equals(that.headers)) {
            return false;
        }
        if (!this.metrics.equals(that.metrics)) {
            return false;
        }
        if (!this.histograms.equals(that.histograms)) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = headers.hashCode();
        result = 31 * result + metrics.hashCode();
        return result;
    }
}
