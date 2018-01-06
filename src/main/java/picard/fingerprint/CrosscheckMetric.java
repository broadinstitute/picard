/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.metrics.MetricBase;

/**
 * A class to hold the result of crosschecking fingerprints.
 * The same metric will be used for crosschecking Readgroups, libraries, samples, or files.
 *
 * @author Yossi Farjoun
 */
public class CrosscheckMetric extends MetricBase {

    // An enum representing whether the result of the fingerprinting was expected and whether it was a match.
    public enum FingerprintResult {
        EXPECTED_MATCH(true, true),
        EXPECTED_MISMATCH(true, false),
        UNEXPECTED_MATCH(false, true),
        UNEXPECTED_MISMATCH(false, false),
        INCONCLUSIVE(null, null);

        private final Boolean isExpected;
        private final Boolean isMatch;

        FingerprintResult(Boolean isExpected, Boolean isMatch) {
            this.isExpected = isExpected;
            this.isMatch = isMatch;
        }

        public Boolean isExpected() {
            return isExpected;
        }

        public Boolean isMatch() {
            return isMatch;
        }
    }

    public enum DataType {
        FILE,
        SAMPLE,
        LIBRARY,
        READGROUP
    }

    public String LEFT_GROUP_VALUE;
    public String RIGHT_GROUP_VALUE;

    // The overall result of the match
    public FingerprintResult RESULT;
    // The data type that was being compared
    public DataType DATA_TYPE;

    // The resulting LOD score comparing LEFT and RIGHT data
    public Double LOD_SCORE;
    // The resulting LOD score comparing LEFT as tumor and RIGHT as normal
    public Double LOD_SCORE_TUMOR_NORMAL;
    // The resulting LOD score comparing LEFT as normal and RIGHT as tumor
    public Double LOD_SCORE_NORMAL_TUMOR;

    // The LEFT run barcode (PU field) expected to look like : D047KACXX110901.1.ACCAACTG
    public String LEFT_RUN_BARCODE;
    // The LEFT lane
    public Integer LEFT_LANE;
    // The LEFT molecular (sample) barcode
    public String LEFT_MOLECULAR_BARCODE_SEQUENCE;
    // The LEFT library identifier
    public String LEFT_LIBRARY;
    // The LEFT sample identifier
    public String LEFT_SAMPLE;
    // The LEFT file from which the fingerprint was obtained
    public String LEFT_FILE;

    // The RIGHT run barcode (PU field) expected to look like : D047KACXX110901.1.ACCAACTG
    public String RIGHT_RUN_BARCODE;
    // The LEFT lane
    public Integer RIGHT_LANE;
    // The LEFT molecular (sample) barcode
    public String RIGHT_MOLECULAR_BARCODE_SEQUENCE;
    // The LEFT library identifier
    public String RIGHT_LIBRARY;
    // The LEFT sample identifier
    public String RIGHT_SAMPLE;
    // The LEFT file from which the fingerprint was obtained
    public String RIGHT_FILE;
}
