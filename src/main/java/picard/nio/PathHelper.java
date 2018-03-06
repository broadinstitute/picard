/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

/**
 * A class who's purpose is to initialize the various plugins that provide Path support.
 * <p>
 * At the moment only google-GCS support is available, but it should be easy to follow the pattern and add support for
 * other Path providers.
 */
public class PathHelper {

    public enum PathProviders {
        GCS {
            @Override
            public void initialize() {
                if (this.isAvailable()) {
                    GoogleStorageUtils.initialize();
                }
            }

            @Override
            public boolean isAvailable() {
                try {
                    new GoogleStorageUtils();
                    return true;
                } catch (NoClassDefFoundError e) {
                    return false;
                }
            }
        };

        public abstract void initialize();

        public abstract boolean isAvailable();
    }

    // making sure this is a singleton class
    private PathHelper() {
    }

    /**
     * calls PathProviders::initialize() for all PathProviders. Needs only be called once,
     * is currently called in PicardCommandlineProgram so that all classes derived from that class need not worry
     * about it.
     */
    public static void initilizeAll() {
        for (PathProviders pathProviders : PathProviders.values()) {
            pathProviders.initialize();
        }
    }
}
