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

package picard.vcf;

import picard.PicardException;
import picard.vcf.GenotypeConcordanceStates.*;

/**
 * The default scheme is derived from the GA4GH Benchmarking Work Group's proposed evaluation scheme. This scheme has been edited to count MISSING
 * sites in the truth set differently for use with truth sets where the HOM_REF sites are not included in the data set. The MISSING truth sites
 * are called like they are HOM_REF sites.
 */

public class GA4GHSchemeWithMissingAsHomRef extends GenotypeConcordanceScheme{
    @Override
    protected void initiateScheme() {
        /**          ROW STATE            MISSING       HOM_REF       HET_REF_VAR1       HET_VAR1_VAR2        HOM_VAR1        NO_CALL        LOW_GQ        LOW_DP        VC_FILTERED   GT_FILTERED   IS_MIXED    **/
        addRow(CallState.MISSING,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_REF,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR1,    FP_TN,        FP_TN,        TP_TN,             TP_FN,               TP_FN,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR2,    NA,           NA,           FP_TN_FN,          NA,                  FP_FN,          NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HET_REF_VAR3,    NA,           NA,           NA,                FP_FN,               NA,             NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HET_VAR1_VAR2,   FP_ONLY,      FP_ONLY,      TP_FP,             TP_ONLY,             TP_FP_FN,       EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_VAR1_VAR3,   NA,           NA,           NA,                TP_FP_FN,            NA,             NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HET_VAR3_VAR4,   FP_ONLY,      FP_ONLY,      FP_FN,             FP_FN,               FP_FN,          NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HOM_VAR1,        FP_ONLY,      FP_ONLY,      TP_FP,             TP_FN,               TP_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_VAR2,        NA,           NA,           FP_FN,             TP_FN,               FP_FN,          NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HOM_VAR3,        NA,           NA,           NA,                FP_FN,               NA,             NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.NO_CALL,         EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.VC_FILTERED,     TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.GT_FILTERED,     TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.LOW_GQ,          TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.LOW_DP,          TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.IS_MIXED,        EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
    }
}
