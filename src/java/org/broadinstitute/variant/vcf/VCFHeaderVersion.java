/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.vcf;

import org.broad.tribble.TribbleException;

/**
 * information that identifies each header version
 */
public enum VCFHeaderVersion {
    VCF3_2("VCRv3.2","format"),
    VCF3_3("VCFv3.3","fileformat"),
    VCF4_0("VCFv4.0","fileformat"),
    VCF4_1("VCFv4.1","fileformat");

    private final String versionString;
    private final String formatString;

    /**
     * create the enum, privately, using:
     * @param vString the version string
     * @param fString the format string
     */
    VCFHeaderVersion(String vString, String fString) {
        this.versionString = vString;
        this.formatString = fString;
    }

    /**
     * get the header version
     * @param version the version string
     * @return a VCFHeaderVersion object
     */
    public static VCFHeaderVersion toHeaderVersion(String version) {
        version = clean(version);
        for (VCFHeaderVersion hv : VCFHeaderVersion.values())
            if (hv.versionString.equals(version))
                    return hv;
        return null;
    }

    /**
     * are we a valid version string of some type
     * @param version the version string
     * @return true if we're valid of some type, false otherwise
     */
    public static boolean isVersionString(String version){
        return toHeaderVersion(version) != null;
    }

    /**
     * are we a valid format string for some type
     * @param format the format string
     * @return true if we're valid of some type, false otherwise
     */
    public static boolean isFormatString(String format){
        format = clean(format);
        for (VCFHeaderVersion hv : VCFHeaderVersion.values())
            if (hv.formatString.equals(format))
                    return true;
        return false;
    }

    public static VCFHeaderVersion getHeaderVersion(String versionLine) {
        String[] lineFields = versionLine.split("=");
        if ( lineFields.length != 2 || !isFormatString(lineFields[0].substring(2)) )
            throw new TribbleException.InvalidHeader(versionLine + " is not a valid VCF version line");

        if ( !isVersionString(lineFields[1]) )
            throw new TribbleException.InvalidHeader(lineFields[1] + " is not a supported version");

        return toHeaderVersion(lineFields[1]);
    }

    /**
     * Utility function to clean up a VCF header string
     * 
     * @param s string
     * @return  trimmed version of s
     */
    private static String clean(String s) {
        return s.trim();
    }


    public String getVersionString() {
        return versionString;
    }

    public String getFormatString() {
        return formatString;
    }
}
