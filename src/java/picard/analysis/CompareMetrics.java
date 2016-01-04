/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.PositionalArguments;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.io.FileReader;
import java.util.List;

/**
 * Compare two metrics files.
 */
@CommandLineProgramProperties(
        usage = CompareMetrics.USAGE,
        usageShort = CompareMetrics.USAGE,
        programGroup = Metrics.class
)
public class CompareMetrics extends CommandLineProgram {

    static final String USAGE = "Compare two metrics files";

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> metricsFiles;

    private static final Log log = Log.getInstance(CompareMetrics.class);

    @Override
    protected int doWork() {
        IOUtil.assertFilesAreReadable(metricsFiles);
        final MetricsFile<?, ?> metricsA = new MetricsFile();
        final MetricsFile<?, ?> metricsB = new MetricsFile();
        try {
            metricsA.read(new FileReader(metricsFiles.get(0)));
            metricsB.read(new FileReader(metricsFiles.get(1)));
            final boolean areEqual = metricsA.areMetricsEqual(metricsB) && metricsA.areHistogramsEqual(metricsB);
            final String status = areEqual ? "EQUAL" : "NOT EQUAL";
            log.info("Files " + metricsFiles.get(0) + " and " + metricsFiles.get(1) + "are " + status);
        } catch (final Exception e) {
            throw new PicardException(e.getMessage());
        }
        return 0;
    }

    public static void main(String[] argv) {
        new CompareMetrics().instanceMainWithExit(argv);
    }
}
