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

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.nio.file.Path;

/**
 * Picard default argument collection for an optional reference.
 */
public class OptionalReferenceArgumentCollection implements ReferenceArgumentCollection {
    private final static Log log = Log.getInstance(OptionalReferenceArgumentCollection.class);

    @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.", common = true, optional = true)
    public PicardHtsPath REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA == null ? null : new PicardHtsPath(Defaults.REFERENCE_FASTA.getAbsolutePath());

    /**
     * @return The reference provided by the user, if any, or the default defined by {@code htsjdk.samtools.Defaults.REFERENCE_FASTA}
     */
    @Override
    public File getReferenceFile() {
        return ReferenceArgumentCollection.getFileSafe(REFERENCE_SEQUENCE, log);
    }

    /**
     * @return The reference provided by the user, if any, or the default, if any, as an nio Path.
     */
    @Override
    public Path getReferencePath() { return REFERENCE_SEQUENCE == null ? null: REFERENCE_SEQUENCE.toPath(); }

    /**
     * @return The reference provided by the user, if any, or the default, if any, as a PicardHtsPath.
     */
    public PicardHtsPath getREFERENCE_SEQUENCE() { return REFERENCE_SEQUENCE; }

}
