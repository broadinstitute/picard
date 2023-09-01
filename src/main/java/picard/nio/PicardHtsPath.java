/*
 * The MIT License
 *
 * Copyright (c) 2022 The Broad Institute
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

package picard.nio;

import htsjdk.io.HtsPath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.utils.ValidationUtils;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * A Subclass of {@link HtsPath} with conversion to {@link Path} making use of {@link IOUtil}
 */
public class PicardHtsPath extends HtsPath {
    /**
     * Create a PicardHtsPath from a raw input path string.
     * <p>
     * If the raw input string already contains a scheme (including a "file" scheme), assume its already
     * properly escape/encoded. If no scheme component is present, assume it references a raw path on the
     * local file system, so try to get a Path first, and then retrieve the URI from the resulting Path.
     * This ensures that input strings that are local file references without a scheme component and contain
     * embedded characters are valid in file names, but which would otherwise be interpreted as excluded
     * URI characters (such as the URI fragment delimiter "#") are properly escape/encoded.
     *
     * @param rawInputString a string specifying an input path. May not be null.
     */
    public PicardHtsPath(final String rawInputString) {
        super(rawInputString);
    }

    /**
     * Create a PicardHtsPath from an existing {@link HtsPath} or subclass.
     *
     * @param htsPath an existing PathSpecifier. May not be null.
     */
    public PicardHtsPath(final HtsPath htsPath) {
        super(htsPath);
    }

    /**
     * Create a PicardHtsPath from a {@link File} reference. Uses the URI string of {@code file}.
     * @param file the file reference to create this object from
     */
    public PicardHtsPath(final File file){
        this(file.toURI().toString());
    }

    /**
     * Create a {@link List<PicardHtsPath>} from path representations.
     * @param paths URIs or local paths. May not be null but may be empty.
     * @return the converted {@link List}
     */
    public static List<PicardHtsPath> fromPaths(Collection<String> paths) {
        Objects.requireNonNull(paths);
        return paths.stream().map(PicardHtsPath::new).collect(Collectors.toList());
    }

    /**
     * Resolve the URI of this object to a {@link Path} object.
     *
     * @return the resulting {@link Path}
     * @throws RuntimeException if an I/O error occurs when creating the file system
     */
    @Override
    public Path toPath() {
        try {
            return IOUtil.getPath(super.getURIString());
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Construct a {@link PicardHtsPath} from a {@link Path}
     * @param path may NOT be null
     * @return a new object representing path
     */
    public static PicardHtsPath fromPath(final Path path){
        Objects.requireNonNull(path);
        return new PicardHtsPath(new HtsPath(path.toUri().toString()));
    }

    /**
     * Create a {@link List<Path>} from {@link PicardHtsPath}s
     * @param picardHtsPaths may NOT be null
     * @return Path representations of the input picardHtsPaths
     */
    public static List<Path> toPaths(final Collection<PicardHtsPath> picardHtsPaths){
        Objects.requireNonNull(picardHtsPaths);
        return picardHtsPaths.stream().map(PicardHtsPath::toPath).collect(Collectors.toList());
    }

    /**
     * Instead of a static method in PicardHtsPath, it could very well be a static method in PicardIOUtils.
     * Or, it should be an instance variable in HtsPath, that returns a new HtsPath object...although we probalby want
     * it to return a PicardHtsPath...similar to fromPath()
     *
     * Examples:
     *     - (test_na12878.bam, .bai) -> test_na12878.bai (append = false)
     *     - (test_na12878.bam, .bai) -> test_na12878.bam.md5 (append = true)
     *
     * @param path the original path
     * @param append whether to append (true) or replace (false) the new extension
     * @param newExtension the extension including the dot e.g. ".txt"
     * @return a new PicardHtsPath object pointed to a file with
     */
    public static PicardHtsPath replaceExtension(final PicardHtsPath path, final String newExtension, final boolean append){
        ValidationUtils.validateArg(newExtension.startsWith("."), "newExtension must start with a dot '.'");

        if (append){
            return new PicardHtsPath(path.getURIString() + newExtension);
        } else {
            final Optional<String> oldExtension = path.getExtension();

            if (oldExtension.isEmpty()){
                throw new PicardException("The extension cannot be identified for the path: " + path.getURIString());
            }

            final String oldFileName = path.toPath().getFileName().toString();
            return PicardHtsPath.fromPath(path.toPath().resolveSibling(oldFileName.replaceAll(oldExtension.get() + "$", newExtension)));
        }
    }

    /**
     *
     * Returns the extension of a filename including the preceding dot '.'
     *
     * e.g. /Users/jsoto/error.log -> ".log"
     *      /Users/jsoto/stderr -> ""
     *
     * As a first pass, we will use the 'lastIndexOf' implementation.
     *
     *
     * **/
//    public static String getExtension(final String uRIString) {
//        final int lastIndexOfPeriod = uRIString.lastIndexOf(".");
//
//        // This is probably not correct --- if say someone is on a Windows machine and giving it a gcloud path
//        final int lastIndexOfSeparator = uRIString.lastIndexOf(FileSystems.getDefault().getSeparator());
//
//        // How about this?
//        String separator;
//        try {
//            separator = FileSystems.getFileSystem(new URI(uRIString)).getSeparator();
//        } catch (URISyntaxException e){
//            // No problem: just use the default
//            separator = FileSystems.getDefault().getSeparator();
//        }
//
//        String a = File.separator;
//
//        String c = System.getProperty("file.separator");
//
//
//        // What if there is no period...
//        String extension = uRIString.substring(uRIString.lastIndexOf("."));
//        return extension;
//    }

    /**
     * Wrapper for Path.resolve()
     */
    public static PicardHtsPath resolve(final PicardHtsPath absPath, final String relativePath){
        return PicardHtsPath.fromPath(absPath.toPath().resolve(relativePath));
    }
}
