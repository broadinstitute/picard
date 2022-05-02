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

import java.io.File;
import java.nio.file.Path;

import htsjdk.samtools.util.Log;
import picard.nio.PicardHtsPath;

/**
 * Base interface for a reference argument collection.
 */
public interface ReferenceArgumentCollection {
    /**
     * @return The reference provided by the user, if any, or the default, if any, as a File. May be null.
     */
    File getReferenceFile();

    /**
     * @return The reference provided by the user, if any, or the default, if any, as an nio Path. May be null.
     */
    default Path getReferencePath() { return getReferenceFile() == null ? null : getReferenceFile().toPath(); }

    /**
     * @return The reference provided by the user, if any, or the default, if any, as a PicardHtsPath. May be null.
     */
    default PicardHtsPath getHtsPath() {
        return getReferenceFile() == null ?
                null :
                new PicardHtsPath(getReferenceFile().getAbsolutePath()); }

    /**
     * @return A "safe" way to obtain a File object for any reference path.
     *
     * For files that reside on a local file system, this returns a valid File object. Files that reside on
     * a non-local file system can't be accessed via a File object, so return a placeholder File object that
     * will defer failure when/if the file is actually accessed. This allows code paths that blindly propagate
     * the value returned by calls to getReferenceFile to not get an NPE, and to fail gracefully downstream
     * with an error message that includes the reference file specifier.
     */
    static File getFileSafe(final PicardHtsPath picardPath, final Log log) {
        if (picardPath == null) {
            return null;
        } else if (picardPath.getScheme().equals(PicardHtsPath.FILE_SCHEME)) {
            // file on a local file system
            return picardPath == null ? null : picardPath.toPath().toFile();
        } else {
            log.warn(String.format(
                    "The reference specified by %s cannot be used as a local file object",
                    picardPath.getRawInputString()));
            // toPath().toFile() would throw here, so use the File constructor to create a placeholder
            // object and defer failure until downstream code attempts to actually use the File object
            return new File(picardPath.getRawInputString());
        }
    }
}
