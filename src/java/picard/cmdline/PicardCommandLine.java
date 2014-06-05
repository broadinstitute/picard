package picard.cmdline;

import htsjdk.samtools.util.Log;
import org.reflections.Reflections;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Set;

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

    public static void main(final String[] args) {
        final Reflections reflections = new Reflections();
        final Set<Class<? extends CommandLineProgram>> allClps =
                reflections.getSubTypesOf(CommandLineProgram.class);
        try {
            final String commandLineProgram = extractCommandLineProgram(args, allClps);
            final Class<?> clpClass = Class.forName(commandLineProgram);
            if (CommandLineProgram.class.isAssignableFrom(clpClass)) {
                final Method mainMethod = clpClass.getMethod("main", String[].class);
                // we can lop off the first two arguments but it requires an array copy or alternatively we could update CLP to remove them
                // in the constructor do the former in this implementation.
                final String[] mainArgs = Arrays.copyOfRange(args, 2, args.length);
                mainMethod.invoke(null, (Object) mainArgs);
            } else {
                log.error("Program with the name " + commandLineProgram +
                        " does not extend CommandLineProgram. Use -h to print out a list of available programs and their descriptions");
                System.exit(1);
            }
        } catch (InstantiationException | ClassNotFoundException | NoSuchMethodException | InvocationTargetException | IllegalAccessException e) {
            log.error("Could not find command line program with the name " + args[1] +
                    ". Use -h to print out a list of available programs and their descriptions");
            System.exit(1);
        }
    }

    private static String extractCommandLineProgram(final String[] args, final Set<Class<? extends CommandLineProgram>> allClps)
            throws InstantiationException, IllegalAccessException {
        if (args.length < 2) {
            printUsage(allClps);
            System.exit(1);
        } else {
            if (!args[0].equals("-t") || args[0].equals("-h")) {
                printUsage(allClps);
                System.exit(1);
            } else {
                for(Class<? extends CommandLineProgram> programClass : allClps){
                    if(programClass.getSimpleName().equals(args[1])){
                        return programClass.getName();
                    }
                }
                printUsage(allClps);
                System.exit(1);
            }
        }
        return "";
    }

    private static void printUsage(final Set<Class<? extends CommandLineProgram>> allClps) throws IllegalAccessException, InstantiationException {
        final StringBuilder builder = new StringBuilder();
        builder.append("USAGE: PicardCommandLine -t <program name>\n\n\tAvailable Programs: \n\t");
        for (final Class<? extends CommandLineProgram> programClass : allClps) {
            final CommandLineProgram program = programClass.newInstance();
            builder.append(program.getCommandLine());
            builder.append(program.getStandardUsagePreamble());
            builder.append("\n\t");
        }
        System.err.println(builder.toString());
    }

}
