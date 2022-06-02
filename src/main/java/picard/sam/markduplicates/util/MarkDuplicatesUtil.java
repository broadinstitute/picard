/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Log;
import htsjdk.utils.ValidationUtils;

public class MarkDuplicatesUtil {
    //singleton.
    private MarkDuplicatesUtil() {
    }

    private static short MIN_INFORMATIVE_MAPPING_QUALITY = -1;
    private static final Log log = Log.getInstance(MarkDuplicatesUtil.class);

    private static boolean warnedForMissingMQ = false;

    private static boolean readAlignedForMarkDuplicates(final SAMRecord rec, final short minInformativeMappingQuality) {
        return !rec.getReadUnmappedFlag() && rec.getMappingQuality() >= minInformativeMappingQuality;
    }

    private static boolean mateAlignedForMarkDuplicates(final SAMRecord rec, final short minInformativeMappingQuality) {
        final Short mateMappingQuality = rec.getShortAttribute(SAMTag.MQ.name());
        if (null == mateMappingQuality) { // cannot rule out, so consider the pair
            if (!warnedForMissingMQ) {
                log.warn("Missing mate mapping quality (MQ) tag in record: " + rec.getSAMString()+ " Future such warnings wll be supressed.");
                warnedForMissingMQ = true;
            }
            return !rec.getMateUnmappedFlag();
        }
        return !rec.getMateUnmappedFlag() && mateMappingQuality >= minInformativeMappingQuality;
    }

    public static boolean pairedForMarkDuplicates(final SAMRecord rec, final short minInformativeMappingQuality) {
        return rec.getReadPairedFlag() &&
                readAlignedForMarkDuplicates(rec, minInformativeMappingQuality) &&
                mateAlignedForMarkDuplicates(rec, minInformativeMappingQuality);
    }

    public static boolean pairedForMarkDuplicates(final SAMRecord rec) {
        ValidationUtils.validateArg(MIN_INFORMATIVE_MAPPING_QUALITY >= 0, "MIN_INFORMATIVE_MAPPING_QUALITY wasn't initialized. Must be >=0.");

        return pairedForMarkDuplicates(rec, MIN_INFORMATIVE_MAPPING_QUALITY);
    }

    public static void setMinInformativeMappingQuality(final short minInformativeMappingQuality) {
        MIN_INFORMATIVE_MAPPING_QUALITY = minInformativeMappingQuality;
    }
}
