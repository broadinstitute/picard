/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble;


/**
 * @author Aaron
 *
 * The base Tribble exception; this allows external libraries to catch any exception Tribble generates
 * 
 */
public class TribbleException extends RuntimeException {
    // what file or input source we are working from
    String source;

    public TribbleException(String msg) {
        super(msg);
    }

    public TribbleException(String message, Throwable throwable) {
        super(message, throwable);
    }

    /**
     * set the source for the file; where we got lines from
     * @param source the source location, usually a file though it could be a http link or other source
     */
    public void setSource(String source) {
        this.source = source;
    }

    /**
     * override the default message with ours, which attaches the source file in question
     * @return a string with our internal error, along with the causitive source file (or other input source)
     */
    public String getMessage() {
        String ret = super.getMessage();
        if ( source != null )
            ret = ret + ", for input source: " + source;
        return ret;
    }

    // //////////////////////////////////////////////////////////////////////
    // other more specific exceptions generated in Tribble
    // //////////////////////////////////////////////////////////////////////


    // //////////////////////////////////////////////////////////////////////
    // Codec exception
    // //////////////////////////////////////////////////////////////////////
    // if the line to decode is incorrect
    public static class InvalidDecodeLine extends TribbleException {
        public InvalidDecodeLine(String message, String line) { super (message + ", line = " + line); }

        public InvalidDecodeLine(String message, int lineNo) { super (message + ", at line number " + lineNo); }
    }

    public static class InvalidHeader extends TribbleException {
        public InvalidHeader(String message) { super ("Your input file has a malformed header: " + message); }
    }

    // capture other internal codec exceptions
    public static class InternalCodecException extends TribbleException {
        public InternalCodecException(String message) { super (message); }
    }

    // //////////////////////////////////////////////////////////////////////
    // Index exceptions
    // //////////////////////////////////////////////////////////////////////
    public static class UnableToCreateCorrectIndexType extends TribbleException {
        public UnableToCreateCorrectIndexType(String message, Exception e) {
            super(message,e);
        }
        public UnableToCreateCorrectIndexType(String message) {
            super(message);
        }
    }

    // //////////////////////////////////////////////////////////////////////
    // Source exceptions
    // //////////////////////////////////////////////////////////////////////
    public static class FeatureFileDoesntExist extends TribbleException {
        public FeatureFileDoesntExist(String message, String file) {
            super(message);
            setSource(file);
        }
    }

    public static class MalformedFeatureFile extends TribbleException {
        public MalformedFeatureFile(String message, String f, Exception e) {
            super(message,e);
            setSource(f);
        }
        public MalformedFeatureFile(String message, String f) {
            super(message);
            setSource(f);
        }
    }

    public static class UnableToReadIndexFile extends TribbleException {
        public UnableToReadIndexFile(String message, String f, Exception e) {
            super(message,e);
            setSource(f);
        }
    }

    public static class TabixReaderFailure extends TribbleException {
        public TabixReaderFailure(String message, String f, Exception e) {
            super(message,e);
            setSource(f);
        }

        public TabixReaderFailure(String message, String f) {
            super(message);
            setSource(f);
        }
    }
}
