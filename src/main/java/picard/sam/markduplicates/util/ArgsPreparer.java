/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

package picard.sam.markduplicates.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class ArgsPreparer {

    private final String[] args;
    private final Set<String> removedArgs;
    private final ArrayList<String> addedArgs;

    public ArgsPreparer(String[] args) {
        this.args = args;
        this.removedArgs = new HashSet<>();
        this.addedArgs = new ArrayList<>();
    }

    public ArgsPreparer exclude(String arg) {
        removedArgs.add(arg);
        return this;
    }

    public ArgsPreparer add(String argKey, String argValue) {
        addedArgs.add(argKey + "=" + argValue);
        return this;
    }

    public String[] toArray() {
        return Arrays.stream(args)
                .filter(arg -> !removedArgs.contains(arg.split("=")[0]))
                .toArray(length -> {
                    String[] args = new String[length + addedArgs.size()];

                    for (int i = 0; i < addedArgs.size(); i++) {
                        args[length + i] = addedArgs.get(i);
                    }

                    return args;
                });
    }
}
