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

package picard.cmdline.argumentcollections;

import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.nio.file.Path;

/**
 * Argument collection for references that are required (and not common).
 */
public class RequiredReferenceArgumentCollection implements ReferenceArgumentCollection {
    private final static Log log = Log.getInstance(RequiredReferenceArgumentCollection.class);

    @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.", common = false)
    public PicardHtsPath REFERENCE_SEQUENCE;

    /**
     * @return The reference provided by the user as a File. May be null.
     */
    public File getReferenceFile() {
        return ReferenceArgumentCollection.getFileSafe(REFERENCE_SEQUENCE, log);
    }

    /**
     * @return The reference provided by the user as an nio Path. May be null.
     */
    @Override
    public Path getReferencePath() { return REFERENCE_SEQUENCE == null ? null: REFERENCE_SEQUENCE.toPath(); }

    /**
     * @return The reference provided by the user as a PicardHtsPath. May be null.
     */
    @Override
    public PicardHtsPath getHtsPath() { return REFERENCE_SEQUENCE; }

}
