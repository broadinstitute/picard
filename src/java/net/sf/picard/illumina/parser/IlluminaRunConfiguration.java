/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Describes the intended logical output configuration of an Illumina run.
 * (e.g. If the input data consists of 80 base
 * clusters and we provide a run configuration of "36T8B36T" then those bases should be split into 3 reads:
 *     read one should be 36 cycles of template,
 *     read two should be 8 cycles of barcode,
 *     read three should be another 36 cycle template read.)
 *  Note: In future releases, IlluminaRunConfigurations will be specified by clients of IlluminaDataProvider(currently
 *  configuration are detected by IlluminaDataProviderFactory via the structure of QSeq files). When working with
 *  QSeq formats, while the individual reads need not fall on QSeq end file boundaries the total number of cycles
 *  should equal the total number of cycles found in all end files for one tile.  (e.g. if I have 80 reads and
 *  3 end files per tile, those end files should have a total of 80 reads in them regardless of how many reads
 *  appear in each individual file)
 *
 *  @author jburke@broadinstitute.org
 */
public class IlluminaRunConfiguration {
    public final List<ReadDescriptor> descriptors;
    public final int totalCycles;
    public final int numBarcodes;
    public final int numTemplates;
    public final int numSkips;
    public final int numDescriptors;

    /** index into descriptors of templates */
    public final int [] templateIndices;

    /** index into descriptors of all barcodes */
    public final int [] barcodeIndices;

    /** index into descriptors of all skips */
    public final int [] skipIndices;

    /** Characters representing valid ReadTypes */
    private static final String ValidTypeChars;

    /** ValidTypeChars except characters are separated by commas */
    private static final String ValidTypeCharsWSep;

    static {
        String validTypes = "";
        String vtWSep = "";

        boolean written = false;
        for(ReadType rt : ReadType.values()) {
            if(written) {
                vtWSep += ",";
            }

            validTypes += rt.name();
            vtWSep += rt.name();
        }

        ValidTypeChars = validTypes;
        ValidTypeCharsWSep = vtWSep;
    }

    private static final String RunConfigMsg = "Run configuration must be formatted as follows: " +
            "<number of bases><type><number of bases><type>...<number of bases> where number of bases is a " +
            "positive (NON-ZERO) integer and type is one of the following characters " + ValidTypeCharsWSep +
            " (e.g. 76T8B68T would denote a paired-end run with a 76 base first end an 8 base barcode followed by a 68 base second end).";
    private static final Pattern FullPattern = java.util.regex.Pattern.compile("^((\\d+[" + ValidTypeChars + "]{1}))+$");
    private static final Pattern SubPattern = java.util.regex.Pattern.compile("(\\d+)([" + ValidTypeChars + "]{1})");

    /**
     * Copies collection into descriptors (making descriptors unmodifiable) and then calculates relevant statistics about descriptors.
     * @param collection A collection of ReadDescriptors that describes this IlluminaRunConfiguration
     */
    public IlluminaRunConfiguration(final List<ReadDescriptor> collection) {
        if(collection.size() == 0) { //If this changes, change hashcode
            throw new IllegalArgumentException("IlluminaRunConfiguration does not support 0 length clusters!");
        }

        this.descriptors = Collections.unmodifiableList(collection);
        int cycles = 0;

        final List<Integer> barcodeIndicesList  = new ArrayList<Integer>();
        final List<Integer> templateIndicesList = new ArrayList<Integer>();
        final List<Integer> skipIndicesList     = new ArrayList<Integer>();

        int descIndex = 0;
        for(final ReadDescriptor desc : descriptors) {
            if(desc.length == 0 || desc.length < 0) {
                throw new IllegalArgumentException("IlluminaRunConfiguration only supports ReadDescriptor lengths > 0, found(" + desc.length + ")");
            }

            cycles += desc.length;
            switch(desc.type) {
                case B:
                    barcodeIndicesList.add(descIndex);
                    break;
                case T:
                    templateIndicesList.add(descIndex);
                    break;
                case S:
                    skipIndicesList.add(descIndex);
                    break;

                default:
                    throw new IllegalArgumentException("Unsupported ReadType (" + desc.type + ") encountered by IlluminaRunConfiugration!");
            }
            ++descIndex;
        }

        this.totalCycles    = cycles;
        this.numBarcodes    = barcodeIndicesList.size();
        this.numTemplates   = templateIndicesList.size();
        this.numSkips       = skipIndicesList.size();
        this.numDescriptors = descriptors.size();

        this.barcodeIndices = new int[numBarcodes];
        for(int i = 0; i < numBarcodes; i++) {
            this.barcodeIndices[i] = barcodeIndicesList.get(i);
        }

        this.templateIndices = new int[numTemplates];
        for(int i = 0; i < numTemplates; i++) {
            this.templateIndices[i] = templateIndicesList.get(i);
        }

        this.skipIndices = new int[numSkips];
        for(int i = 0; i < numSkips; i++) {
            this.skipIndices[i] = skipIndicesList.get(i);
        }
    }

    /**
     * Converts configStr into a List<ReadDescriptor> and calls the primary constructor using this List as it's argument.
     * @param configStr A string of the format <number of bases><type><number of bases><type>...<number of bases><type> describing
     * this run configuration
     */
    public IlluminaRunConfiguration(final String configStr) {
        this(configStrToDescriptors(configStr));
    }

    /**
     * Converts this object into a String using rules complementary to the single string constructor above.
     * @return A string of the form <number of bases><type><number of bases><type>...<number of bases><type> with one
     * <number of bases><type> per ReadDescriptor in descriptors.
     */
    @Override
    public String toString() {
        String out = "";
        for(final ReadDescriptor rd : descriptors) {
            out += rd.toString();
        }
        return out;
    }

    /**
     * Converts configStr into a List<ReadDescriptor>
     * @param configStr A string of the format <number of bases><type><number of bases><type>...<number of bases><type> describing
     * a run configuration
     * @return A List<ReadDescriptor> corresponding to the input string
     */
    private final static List<ReadDescriptor> configStrToDescriptors(final String configStr) {
        final Matcher fullMatcher = FullPattern.matcher(configStr);
        if(!fullMatcher.matches()) {
            throw new IllegalArgumentException(configStr + " cannot be parsed as an Illumina Run Configuration! " + RunConfigMsg);
        }


        final Matcher subMatcher = SubPattern.matcher(configStr);
        final List<ReadDescriptor> descriptors = new ArrayList<ReadDescriptor>();
        while(subMatcher.find()) {
            final ReadDescriptor rd =  new ReadDescriptor(Integer.parseInt(subMatcher.group(1)), ReadType.valueOf(subMatcher.group(2)));
            descriptors.add(rd);
        }

        return descriptors;
    }

    @Override
    public boolean equals(final Object thatObj) {
        if(this == thatObj) return true;
        if(this.getClass() != thatObj.getClass()) return false;

        final IlluminaRunConfiguration that = (IlluminaRunConfiguration) thatObj;
        if(this.descriptors.size() != that.descriptors.size()) {
            return false;
        }

        for(int i = 0; i < this.descriptors.size(); i++) {
            if(!this.descriptors.get(i).equals(that.descriptors.get(i))) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int hashCode() {
        int res = descriptors.get(0).hashCode();
        for(int i = 1; i < descriptors.size(); i++) {
            res *= descriptors.get(i).hashCode();
        }

        return res;
    }
}
