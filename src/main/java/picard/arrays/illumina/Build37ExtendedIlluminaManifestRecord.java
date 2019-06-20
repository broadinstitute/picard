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

package picard.arrays.illumina;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang.StringUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * A class to represent a record (line) from an Extended Illumina Manifest [Assay] entry
 */
public class Build37ExtendedIlluminaManifestRecord extends IlluminaManifestRecord {
    protected enum Flag {
        ILLUMINA_FLAGGED,   // Illumina flagged
        LIFTOVER_FAILED,
        UNSUPPORTED_GENOME_BUILD,
        SEQUENCE_MISMATCH,      // mismatch between what the manifest claims is the reference vs. the actual reference.
        INDEL_SEQ_MISMATCH,     // Unable to reconcile indel situation
        INDEL_EXTENSION_ERROR,  // extension base conflict.
        DUPE,
        PASS,
    }

    private String b37Chr;
    private Integer b37Pos;
    private String snpRefAllele;
    private String snpAlleleA;
    private String snpAlleleB;
    private String rsId;
    private Flag flag = Flag.PASS;

    private Allele A;
    private Allele B;
    private Allele ref;
    private final Log log = Log.getInstance(Build37ExtendedIlluminaManifestRecord.class);

    /**
     * This constructor is used to read records from an already created Build37ExtendedIlluminaManifestRecord file.
     * It does not work to set the Extended-specific fields
     */
    Build37ExtendedIlluminaManifestRecord(final Map<String, Integer> columnNameToIndex, final String[] line, final int index) {
        super(columnNameToIndex, line, index);

        final int end = line.length;
        flag = Flag.valueOf(line[end - 1]);

        if (!isBad()) {
            b37Chr = line[end - 7];
            b37Pos = parseIntOrNull(line[end - 6]);
            snpRefAllele = line[end - 5];
            snpAlleleA = line[end - 4];
            snpAlleleB = line[end - 3];
            rsId = line[end - 2];

            A = Allele.create(snpAlleleA, snpAlleleA.equals(snpRefAllele));
            B = Allele.create(snpAlleleB, snpAlleleB.equals(snpRefAllele));
            ref = Allele.create(snpRefAllele, true);
        } else {
            b37Chr = "0";
            b37Pos = 0;
            snpRefAllele = "";
            snpAlleleA = "";
            snpAlleleB = "";
            rsId = "";

            A = Allele.NO_CALL;
            B = Allele.NO_CALL;
            ref = Allele.NO_CALL;
        }
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }

    public Allele getRefAllele() {
        return ref;
    }

    public String getB37Chr() {
        return b37Chr;
    }

    public Integer getB37Pos() {
        return b37Pos;
    }

    public String getRsId() { return rsId; }

    public Boolean isBad() {
        return flag != Flag.DUPE && flag != Flag.PASS;
    }

    public Boolean isDupe() {
        return flag == Flag.DUPE;
    }

    public Flag getFlag() {
        return flag;
    }

    @Override
    public String getLine() {
        final String originalLine = super.getLine();

        final List<String> extensions = new ArrayList<>();
        extensions.add(b37Chr);
        extensions.add(b37Pos != null ? b37Pos.toString() : null);
        extensions.add(snpRefAllele);
        extensions.add(snpAlleleA);
        extensions.add(snpAlleleB);
        extensions.add(rsId);
        extensions.add(flag.name());

        return originalLine + "," + StringUtils.join(extensions, ",");
    }

}