
/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

package picard.fingerprint;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Fingerprinting;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Program to check that all read groups within the set of BAM files appear to come from the same
 * individual.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "DEPRECATED: USE CrosscheckFingerprints. Checks if all read groups within a set of BAM files appear to come from the same individual",
        usageShort = "DEPRECATED: USE CrosscheckFingerprints. Checks if all read groups appear to come from the same individual.",
        programGroup = Fingerprinting.class
)
/**
 @deprecated 6/6/2017
 Use  {@link CrosscheckFingerprints} instead.
 */
@Deprecated
public class CrosscheckReadGroupFingerprints extends CrosscheckFingerprints {


    @Option(doc = "Instead of producing the normal comparison of read-groups, roll fingerprints up to the sample level " +
            "and print out a sample x sample matrix with LOD scores.")
    public boolean CROSSCHECK_SAMPLES = false;

    @Option(doc = "Instead of producing the normal comparison of read-groups, roll fingerprints up to the library level " +
            "and print out a library x library matrix with LOD scores.")
    public boolean CROSSCHECK_LIBRARIES = false;

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();

        if (CROSSCHECK_BY != CrosscheckMetric.DataType.READGROUP) {
            errors.add("When calling CrosscheckReadGroupFingerprints, please refrain from supplying a CROSSCHECK_BY argument. " +
                            "(Found value " + CROSSCHECK_BY + "\n" +
                    "Use CrosscheckFingerprints if you would like to do that.");
        }

        if (MATRIX_OUTPUT != null) {
            errors.add("When calling CrosscheckReadGroupFingerprints, please refrain from supplying a MATRIX_OUTPUT argument.\n" +
                    "(Found value " + MATRIX_OUTPUT + "\n" +
                    "Use CrosscheckFingerprints if you would like to do that.");
        }

        if (!errors.isEmpty()) {
            return errors.toArray(new String[errors.size()]);
        } else {
            return super.customCommandLineValidation();
        }
    }

    @Override
    protected int doWork() {

        if (CROSSCHECK_LIBRARIES) {
            CROSSCHECK_BY = CrosscheckMetric.DataType.LIBRARY;
        } else if (CROSSCHECK_SAMPLES) {
            CROSSCHECK_BY = CrosscheckMetric.DataType.SAMPLE;
        } else {
            CROSSCHECK_BY = CrosscheckMetric.DataType.READGROUP;
            MATRIX_OUTPUT = OUTPUT;
            OUTPUT = new File("/dev/null");
        }

        return super.doWork();
    }
}
     