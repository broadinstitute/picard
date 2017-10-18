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

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * class to hold the details of a element of fingerprinting PU tag
 *
 * @author Yossi Farjoun
 */

public class FingerprintIdDetails {
    String platformUnit;
    String runBarcode;
    Integer runLane;
    String molecularBarcode;
    String library;
    String file;
    String sample;

    String group;// not used for equals or hash since it's expected to be made up of one of the other fields..

    static final String multipleValuesString = "<MULTIPLE_VALUES>";

    public FingerprintIdDetails() {}

    // If platformUnit is not populated, or "improperly" formatted (or missing), then fields will be initialized as
    // flowcellBarcode="?", lane=-1, molecularBarcode="?"
    public FingerprintIdDetails(final String platformUnit, final String file) {
        getPlatformUnitDetails(platformUnit);
        this.platformUnit = platformUnit;
        this.file = file;
    }

    public FingerprintIdDetails(final SAMReadGroupRecord rg, final String file) {
        this(rg.getPlatformUnit(), file);
        this.sample = rg.getSample();
        this.library = rg.getLibrary();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final FingerprintIdDetails that = (FingerprintIdDetails) o;

        if (platformUnit != null ? !platformUnit.equals(that.platformUnit) : that.platformUnit != null) return false;
        if (runBarcode != null ? !runBarcode.equals(that.runBarcode) : that.runBarcode != null) return false;
        if (runLane != null ? !runLane.equals(that.runLane) : that.runLane != null) return false;
        if (molecularBarcode != null ? !molecularBarcode.equals(that.molecularBarcode) : that.molecularBarcode != null)
            return false;
        if (library != null ? !library.equals(that.library) : that.library != null) return false;
        if (file != null ? !file.equals(that.file) : that.file != null) return false;

        return sample != null ? sample.equals(that.sample) : that.sample == null;
    }

    @Override
    public int hashCode() {
        int result = platformUnit != null ? platformUnit.hashCode() : 0;
        result = 31 * result + (runBarcode != null ? runBarcode.hashCode() : 0);
        result = 31 * result + (runLane != null ? runLane.hashCode() : 0);
        result = 31 * result + (molecularBarcode != null ? molecularBarcode.hashCode() : 0);
        result = 31 * result + (library != null ? library.hashCode() : 0);
        result = 31 * result + (file != null ? file.hashCode() : 0);
        result = 31 * result + (sample != null ? sample.hashCode() : 0);
        return result;
    }

    public FingerprintIdDetails merge(final FingerprintIdDetails other) {

        platformUnit     = equalValueOrElse(platformUnit,     other.platformUnit,     multipleValuesString);
        runBarcode       = equalValueOrElse(runBarcode,       other.runBarcode,       multipleValuesString);
        runLane          = equalValueOrElse(runLane,          other.runLane,          Integer.MIN_VALUE);
        library          = equalValueOrElse(library,          other.library,          multipleValuesString);
        file             = equalValueOrElse(file,             other.file,             multipleValuesString);
        sample           = equalValueOrElse(sample,           other.sample,           multipleValuesString);
        molecularBarcode = equalValueOrElse(molecularBarcode, other.molecularBarcode, multipleValuesString);

        return this;
    }

    static private <T> T equalValueOrElse(final T lhs, final T rhs, final T orElse) {
        if (rhs == null) return lhs;
        if (lhs == null) return rhs;

        return lhs.equals(rhs) ? lhs : orElse;
    }

    public String getPlatformUnit() {
        return platformUnit;
    }

    public String getSample() {
        return sample;
    }

    /**
     * Fills the relevant fields from the platformUnit string.
     *
     * @param puString platform Unit tag (from @RG) under consideration
     */
    private void getPlatformUnitDetails(final String puString) {

        this.runBarcode = "?";
        this.runLane = -1;
        this.molecularBarcode = "?";

        if (puString == null) return;

        final String[] tmp = puString.split("\\."); // Expect to look like: D047KACXX110901.1.ACCAACTG
        if ((tmp.length == 3) || (tmp.length == 2)) {
            this.runBarcode = tmp[0];
            this.molecularBarcode = (tmp.length == 3) ? tmp[2] : "";  // In older BAMS there may be no molecular barcode sequence
            try {
                this.runLane = Integer.parseInt(tmp[1]);
            } catch (final NumberFormatException e) {
                //no-op. return with whatever was parsed
            }
        }
    }
}
