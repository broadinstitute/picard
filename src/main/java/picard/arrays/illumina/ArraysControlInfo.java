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

/**
 * A simple class to store names and counts for the the Control Information fields that are stored in an Illumina GTC file.
 */
public class ArraysControlInfo {

    private final String control;
    private final String category;
    private final int red;
    private final int green;

    public ArraysControlInfo(String control, String category, int red, int green) {
        this.control = control;
        this.category = category;
        this.red = red;
        this.green = green;
    }

    public static ArraysControlInfo[] CONTROL_INFO = new ArraysControlInfo[]{
            new ArraysControlInfo("DNP(High)", "Staining", 0, 0),
            new ArraysControlInfo("DNP(Bgnd)", "Staining", 0, 0),
            new ArraysControlInfo("Biotin(High)", "Staining", 0, 0),
            new ArraysControlInfo("Biotin(Bgnd)", "Staining", 0, 0),
            new ArraysControlInfo("Extension(A)", "Extension", 0, 0),
            new ArraysControlInfo("Extension(T)", "Extension", 0, 0),
            new ArraysControlInfo("Extension(C)", "Extension", 0, 0),
            new ArraysControlInfo("Extension(G)", "Extension", 0, 0),
            new ArraysControlInfo("TargetRemoval", "TargetRemoval", 0, 0),
            new ArraysControlInfo("Hyb(High)", "Hybridization", 0, 0),
            new ArraysControlInfo("Hyb(Medium)", "Hybridization", 0, 0),
            new ArraysControlInfo("Hyb(Low)", "Hybridization", 0, 0),
            new ArraysControlInfo("String(PM)", "Stringency", 0, 0),
            new ArraysControlInfo("String(MM)", "Stringency", 0, 0),
            new ArraysControlInfo("NSB(Bgnd)Red", "Non-SpecificBinding", 0, 0),
            new ArraysControlInfo("NSB(Bgnd)Purple", "Non-SpecificBinding", 0, 0),
            new ArraysControlInfo("NSB(Bgnd)Blue", "Non-SpecificBinding", 0, 0),
            new ArraysControlInfo("NSB(Bgnd)Green", "Non-SpecificBinding", 0, 0),
            new ArraysControlInfo("NP(A)", "Non-Polymorphic", 0, 0),
            new ArraysControlInfo("NP(T)", "Non-Polymorphic", 0, 0),
            new ArraysControlInfo("NP(C)", "Non-Polymorphic", 0, 0),
            new ArraysControlInfo("NP(G)", "Non-Polymorphic", 0, 0),
            new ArraysControlInfo("Restore", "Restoration", 0, 0)};

    public String getCategory() {
        return category;
    }

    public String getControl() {
        return control;
    }

    @Override
    public String toString() {
        return control + "|" + category;
    }

    public int getGreen() {
        return green;
    }

    public int getRed() {
        return red;
    }

    public String fullString() {
        return toString() + "|" + red + "|" + green;
    }
}

