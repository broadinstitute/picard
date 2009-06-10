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

/**
 * Contains a set of metrics that can be written to a file and parsed back
 * again. The set of metrics is composed of zero or more instances of a class,
 * BEAN, that extends {@link MetricBase} (all instances must be of the same type)
 * and may optionally include a histogram of data.
 *
 * @author Tim Fennell
 */
public class MetricsFile<BEAN extends MetricBase, HKEY extends Comparable> {
    public static final String MAJOR_HEADER_PREFIX = "## ";
    public static final String MINOR_HEADER_PREFIX = "# ";
    public static final String SEPARATOR = "\t";
    public static final String HISTO_HEADER = "## HISTOGRAM\t";
    public static final String METRIC_HEADER = "## METRICS CLASS\t";

    private List<Header> headers = new ArrayList<Header>();
    private List<BEAN> metrics = new ArrayList<BEAN>();
    private Histogram<HKEY> histogram;

    /** Adds a header to the collection of metrics. */
    public void addHeader(Header h) { this.headers.add(h); }

    /** Returns the list of headers. */
    public List<Header> getHeaders() { return Collections.unmodifiableList(this.headers); }

    /** Adds a bean to the collection of metrics. */
    public void addMetric(BEAN bean) { this.metrics.add(bean); }

    /** Returns the list of headers. */
    public List<BEAN> getMetrics() { return Collections.unmodifiableList(this.metrics); }

    /** Returns the histogram contained in the metrics file if any. */
    public Histogram<HKEY> getHistogram() { return histogram; }

    /** Sets the histogram contained in the metrics file. */
    public void setHistogram(Histogram<HKEY> histogram) { this.histogram = histogram; }

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
        if (this.histogram == null || this.histogram.isEmpty()) {
            return;
        }

        // Add a header for the histogram key type
        out.append(HISTO_HEADER + this.histogram.keySet().iterator().next().getClass().getName());
        out.newLine();
        
        if (this.histogram != null) {
            out.append(StringUtil.assertCharactersNotInString(this.histogram.getBinLabel(), '\t', '\n'));
            out.append(SEPARATOR);
            out.append(StringUtil.assertCharactersNotInString(this.histogram.getValueLabel(), '\t', '\n'));
            out.newLine();
            
            for (Histogram<HKEY>.Bin bin : this.histogram.values()) {
                out.append(StringUtil.assertCharactersNotInString(formatter.format(bin.getId()), '\t', '\n'));
                out.append(MetricsFile.SEPARATOR);
                out.append(formatter.format(bin.getValue()));
                out.newLine();
            }
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
                        header = (Header) Class.forName(className).newInstance();
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
                    type = Class.forName(className);
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

            // Then read the histogram if it is present
            while (line != null && !line.startsWith(MAJOR_HEADER_PREFIX)) {
                line = in.readLine();
            }
            if (line != null && line.startsWith(HISTO_HEADER)) {
                // Get the key type of the histogram
                String keyClassName = line.split(SEPARATOR)[1].trim();
                Class<?> keyClass = null;

                try { keyClass = Class.forName(keyClassName); }
                catch (ClassNotFoundException cnfe) { throw new PicardException("Could not load class with name " + keyClassName); }

                // Read the next line with the bin and value labels
                String[] labels = in.readLine().split(SEPARATOR);
                this.histogram = new Histogram(labels[0], labels[1]);

                // Read the entries in the histogram
                while ((line = in.readLine()) != null && !"".equals(line)) {
                    String[] fields = line.trim().split(SEPARATOR);
                    HKEY key = (HKEY) formatter.parseObject(fields[0], keyClass);
                    double value = formatter.parseDouble(fields[1]);
                    this.histogram.increment(key, value);
                }
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Could not read metrics from reader.", ioe);
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
        if (this.histogram == null && that.histogram == null) {
            return true;
        } else if (this.histogram != null) {
            return this.histogram.equals(that.histogram);
        } else if (that.histogram != null) {
            return that.histogram.equals(this.histogram);
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = headers.hashCode();
        result = 31 * result + metrics.hashCode();
        result = 31 * result + (histogram != null ? histogram.hashCode() : 0);
        return result;
    }
}
