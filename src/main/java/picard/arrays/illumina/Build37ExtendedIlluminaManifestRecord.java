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

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.StringUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * A class to represent a record (line) from an Extended Illumina Manifest [Assay] entry
 */
public class Build37ExtendedIlluminaManifestRecord extends IlluminaManifestRecord {
    protected enum Flag {
        /** The record in the manifest passes validation and will be used in VCF generation */
        PASS,

        /** The record passes validation but is a duplicate (another assay at the same locus with the same alleles) */
        DUPE,

        /** Flagged by Illumina as a bad assay */
        ILLUMINA_FLAGGED,

        LIFTOVER_FAILED,

        UNSUPPORTED_GENOME_BUILD,

        /** Probe sequence not found in reference. */
        PROBE_SEQUENCE_MISMATCH,

        /** Manifest contained no allele B probe sequence */
        MISSING_ALLELE_B_PROBESEQ,

        /** Source sequence is invalid (contains invalid character). */
        SOURCE_SEQUENCE_INVALID,

        /** Neither insertion nor deletion sequence found in reference. */
        INDEL_NOT_FOUND,

        /** Both insertion and deletion sequence found in reference. */
        INDEL_CONFLICT,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        PROBE_SEQUENCE_STRAND_INVALID,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        SOURCE_SEQUENCE_MISMATCH,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        SOURCE_SEQUENCE_STRAND_INVALID,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        SEQUENCE_MISMATCH,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        INDEL_SEQ_MISMATCH,

        /** @deprecated - but used in existing extended manifest files. */
        @Deprecated
        INDEL_EXTENSION_ERROR
    }

    String b37Chr;
    Integer b37Pos;
    String snpRefAllele;
    String snpAlleleA;
    String snpAlleleB;
    String rsId;
    Flag flag = Flag.PASS;

    Allele aAllele = null;
    Allele bAllele = null;
    Allele refAllele = null;

    // The refStrand if provided in the Illumina manifest, otherwise calculated
    Strand referenceStrand = null;

    /**
     * This constructor is used to read records from an already created Build37ExtendedIlluminaManifestRecord file.
     */
    Build37ExtendedIlluminaManifestRecord(final Map<String, Integer> columnNameToIndex, final String[] line, final int index) {
        super(columnNameToIndex, line, index);

        final int end = line.length;
        flag = Flag.valueOf(line[end - 1]);

        if (!isFail()) {
            b37Chr = line[end - 7];
            b37Pos = parseIntOrNull(line[end - 6]);
            snpRefAllele = line[end - 5];
            snpAlleleA = line[end - 4];
            snpAlleleB = line[end - 3];
            rsId = line[end - 2];
        } else {
            b37Chr = "0";
            b37Pos = 0;
            snpRefAllele = "";
            snpAlleleA = "";
            snpAlleleB = "";
            rsId = "";
        }
    }

    Build37ExtendedIlluminaManifestRecord(final IlluminaManifestRecord record,
                                   final Flag flag,
                                   final String b37Chr,
                                   final Integer b37Pos,
                                   final String snpRefAllele,
                                   final String snpAlleleA,
                                   final String snpAlleleB,
                                   final String rsId) {
        super(record);
        this.flag = flag;
        this.b37Chr = b37Chr;
        this.b37Pos = b37Pos;
        this.snpRefAllele = snpRefAllele;
        this.snpAlleleA = snpAlleleA;
        this.snpAlleleB = snpAlleleB;
        this.rsId = rsId;
    }

    public Allele getAlleleA() {
        if (aAllele == null) {
            aAllele = Allele.NO_CALL;
            if (!isFail() && !StringUtils.isEmpty(snpAlleleA)) {
                aAllele = Allele.create(snpAlleleA, snpAlleleA.equals(snpRefAllele));
            }
        }
        return aAllele;
    }

    public Allele getAlleleB() {
        if (bAllele == null) {
            bAllele = Allele.NO_CALL;
            if (!isFail() && !StringUtils.isEmpty(snpAlleleB)) {
                bAllele = Allele.create(snpAlleleB, snpAlleleB.equals(snpRefAllele));
            }
        }
        return bAllele;
    }

    public Allele getRefAllele() {
        if (refAllele == null) {
            refAllele = Allele.NO_CALL;
            if (!isFail() && !StringUtils.isEmpty(snpRefAllele)) {
                refAllele = Allele.create(snpRefAllele, true);
            }
        }
        return refAllele;
    }

    public Strand getReferenceStrand() { return referenceStrand; }

    public String getB37Chr() {
        return b37Chr;
    }

    public Integer getB37Pos() {
        return b37Pos;
    }

    public String getSnpRefAllele() {
        return snpRefAllele;
    }

    public String getSnpAlleleA() {
        return snpAlleleA;
    }

    public String getSnpAlleleB() {
        return snpAlleleB;
    }

    public String getRsId() { return rsId; }

    public Boolean isFail() {
        return flag != Flag.DUPE && flag != Flag.PASS;
    }

    public Boolean isDupe() {
        return flag == Flag.DUPE;
    }

    public Flag getFlag() {
        return flag;
    }

    public void setRsId(String rsId) {
        this.rsId = rsId;
    }

    public void setDupe(boolean isDupe) {
        if (!isFail()) {
            if (isDupe) {
                flag = Flag.DUPE;
            }
        }
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

        return originalLine + "," + String.join(",", extensions);
    }
}
