/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.samtools;


import net.sf.samtools.util.CloseableIterator;
import java.io.*;


/**
 * Command line utility for manipulating SAM/BAM files.
 */
public class SAMTools
{
    private String mCommand = null;
    private File mInputFile = null;


    public static void main(final String[] args)
        throws Exception {
        final int status = new SAMTools().run(args);
        if (status != 0) {
            System.exit(status);
        }
    }

    private SAMTools() {
    }

    private void usage() {
        System.out.println();
        System.out.println("SAMTools version 0.1.0");
        System.out.println("Tools for manipulating SAM/BAM files");
        System.out.println();
        System.out.println("Usage: SAMTools <command> <options...>");
        System.out.println();
        System.out.println("Commands:");
        System.out.println("  help");
        System.out.println("  view        <file>");
        System.out.println();
    }

    private boolean parseArguments(final String[] args) {
        if (args.length == 0) {
            usage();
            return true;
        }
        final String command = args[0];
        final int argpos = 1;
        final int argcount = args.length - argpos;
        if (command.equals("help")) {
            usage();
            return true;
        } else if (command.equals("view")) {
            if (argcount != 1) {
                usage();
                return false;
            }
            mInputFile = new File(args[1]);
            if (!mInputFile.exists()) {
                System.out.println("Input file not found: " + mInputFile);
                return false;
            }
        } else {
            System.out.println("Unrecognized command: " + command);
            System.out.println();
            usage();
            return false;
        }
        mCommand = command;
        return true;
    }

    private int run(final String[] args)
        throws Exception {
        if (!parseArguments(args)) {
            return 1;
        }
        if (mCommand == null) {
            return 0;
        }
        if (mCommand.equals("view")) {
            return runView();
        }
        return 1;
    }

    private int runView() {
        final SAMFileReader reader = new SAMFileReader(mInputFile);
        final CloseableIterator<SAMRecord> iterator = reader.iterator();
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            System.out.println(record.format());
        }
        iterator.close();
        return 0;
    }
}
