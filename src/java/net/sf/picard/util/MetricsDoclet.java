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

package net.sf.picard.util;

import com.sun.javadoc.*;

import net.sf.picard.metrics.MetricBase;

import java.util.TreeMap;
import java.util.SortedMap;
import java.io.PrintStream;
import java.io.File;
import java.io.FileNotFoundException;

/**
 * Doclet for use with JavaDoc that will find all classes extending MetricBase and
 * output information about the metrics definitions that go along with the classes.
 *
 * Takes a single parameter (-f file) to tell it where to output the resulting
 * documentation file in HTML format.
 *
 * @author Tim Fennell
 */
public class MetricsDoclet {
    /**
     * Entry point called by the javadoc command line tool. Loops over all the
     * classes identifying metrics classes and then produces some basic information
     * about each in a single HTML file.
     *
     * @param root the root of the javadoc object hierarchy
     * @return true if completed successfully, false otherwise
     */
    public static boolean start(final RootDoc root) {
        // Build a set of metrics classes sorted by name
        final SortedMap<String,ClassDoc> metricsClasses = new TreeMap<String,ClassDoc>();
        for (final ClassDoc doc : root.classes()) {
            if (isMetricsClass(doc)) {
                System.out.println("Processing " + doc.qualifiedTypeName());
                metricsClasses.put(doc.typeName(), doc);
            }
        }

        // Get a print stream to write to
        final PrintStream out = getOutput(root);
        if (out == null) return false;

        // Write the headings
        out.println("<head>");
        out.println("  <title>Picard Metrics Definitions</title>");
        out.println("  <style>");
        out.println("    body { font-family: Arial; background-color: white; font-size: 10pt; }");
        out.println("    h2 { color: red; }");
        out.println("    .class_description { font-style: italic; }");
        out.println("    .metric_column_def { padding-bottom: 0px; margin-bottom: 8px; }");
        out.println("    .metric_column_def span { font-weight: bold; color: #222222; font-style: italic; }");
        out.println("  </style>");
        out.println("</head>");
        out.println("<h1>Picard Metrics Definitions</h1>");

        // Write out the TOC
        out.println("<h2>Table Of Contents</h2>");
        out.println("<ol>");
        for (final ClassDoc doc : metricsClasses.values()) {
            out.println("<li><a href=\"#" + doc.name() + "\">" + doc.name() + "</a>: " +
                        firstSentence(doc) + "</li>");
        }
        out.println("</ol>");

        // Now print out each class
        for (final ClassDoc doc : metricsClasses.values()) {
            out.println("<a name=\"" + doc.name() + "\"");
            out.println("<h2>" + doc.name() + "</h2>");
            out.println("<p class=\"class_description\">" + doc.commentText() + "</p>");
            out.println("<h3>Column Definitions</h3>");

            for (final FieldDoc field : doc.fields()) {
                if (field.isPublic() && !field.isStatic()) {
                    out.append("<div class=\"metric_column_def\"><span>" + field.name() + ": </span>");
                    out.append(field.commentText());
                    out.println("</div>");
                }
            }
        }

        out.close();
        return true;
    }

    /**
     * Checks to see if the class extends MetricBase using only the JavaDoc
     * metadata provided about the class.
     *
     * @param doc the ClassDoc representing the class to be tested
     * @return true if the class is a metrics class, false otherwise
     */
    protected static boolean isMetricsClass(ClassDoc doc) {
        final String metricBaseFqn = MetricBase.class.getName();
        if (!doc.isClass()) return false;
        if (doc.qualifiedTypeName().contains("personal")) return false;

        do {
            doc = doc.superclass();
            if (doc != null && metricBaseFqn.equals(doc.qualifiedTypeName())) return true;
        }
        while (doc != null);

        return false;
    }

    /**
     * Gets the file output parameter from the RootDoc and then opens an
     * PrintStream to write to the file.
     */
    protected static PrintStream getOutput(final RootDoc root) {
        for (final String[] arg : root.options()) {
            if (arg[0].equals("-f") && arg.length == 2) {
                try {
                    return new PrintStream(new File(arg[1]));
                }
                catch (FileNotFoundException fnfe) {
                    root.printError("Could not open destination file: " + arg[1]);
                    fnfe.printStackTrace();
                    return null;
                }
            }
        }

        root.printError("Destination file parameter -f not supplied.");
        return null;
    }

    /**
     * Required method by the javadoc caller that returns the expected number of elements
     * for doclet specific command line arguments.
     */
    public static int optionLength(final String option) {
        if(option.equals("-f")) {
	        return 2;
        }
        return 0;
    }

    /**
     * Takes a Doc object and uses the firstSentenceTags() to recreate the first sentence
     * text.
     */
    protected static String firstSentence(final Doc doc) {
        final Tag[] tags = doc.firstSentenceTags();
        final StringBuilder builder = new StringBuilder(128);
        for (final Tag tag : tags) {
            builder.append(tag.text());
        }

        return builder.toString();
    }

}
