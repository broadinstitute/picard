package picard.cmdline;

import htsjdk.samtools.util.Log;

import java.util.Arrays;
import java.util.ServiceLoader;

/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
public class PicardCommandLine {
    private static final Log log = Log.getInstance(PicardCommandLine.class);

    public int instanceMain(final String[] args) {
        final CommandLineProgram program = extractCommandLineProgram(args);
        // we can lop off the first two arguments but it requires an array copy or alternatively we could update CLP to remove them
        // in the constructor do the former in this implementation.
        final String[] mainArgs = Arrays.copyOfRange(args, 2, args.length);
        return program.instanceMain(mainArgs);
    }

    public static void main(final String[] args) {
        System.exit(new PicardCommandLine().instanceMain(args));
    }

    private static CommandLineProgram extractCommandLineProgram(final String[] args) {
        final ServiceLoader<CommandLineProgram> loader = ServiceLoader.load(CommandLineProgram.class);
        if (args.length < 2) {
            printUsage(loader);
            System.exit(1);
        } else {
            if (!args[0].equals("-T") || args[0].equals("-h")) {
                printUsage(loader);
                System.exit(1);
            } else {

                for (final CommandLineProgram program : loader) {
                    if (program.getClass().getSimpleName().equals(args[1])) {
                        return program;
                    }
                }
                printUsage(loader);
                System.exit(1);
            }
        }
        return null;
    }

    private static void printUsage(final ServiceLoader<CommandLineProgram> loader) {
        final StringBuilder builder = new StringBuilder();
        builder.append("USAGE: PicardCommandLine -T <program name> [-h]\n\n");
        builder.append("Available Programs:\n");
        builder.append("--------------------------------------------------------------------------------------\n");
        for (final CommandLineProgram program : loader) {
            builder.append(String.format(" %-45s%s\n", program.getClass().getSimpleName(), program.getShortUsage()));
        }
        System.err.println(builder.toString());
    }
}
