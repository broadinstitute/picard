package picard.cmdline;

import htsjdk.samtools.util.Log;
import org.reflections.Reflections;
import org.reflections.scanners.FieldAnnotationsScanner;
import org.reflections.util.FilterBuilder;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;

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
        try {
            final String commandLineProgram = extractCommandLineProgram(args);
            final Class<?> clpClass = Class.forName(commandLineProgram);
            if (CommandLineProgram.class.isAssignableFrom(clpClass)) {
                final CommandLineProgram program = (CommandLineProgram) clpClass.newInstance();
                final Method mainMethod = clpClass.getMethod("instanceMain", String[].class);
                // we can lop off the first two arguments but it requires an array copy or alternatively we could update CLP to remove them
                // in the constructor do the former in this implementation.
                final String[] mainArgs = Arrays.copyOfRange(args, 2, args.length);
                return (Integer) mainMethod.invoke(program, (Object) mainArgs);
            } else {
                log.error("Program with the name " + commandLineProgram +
                        " does not extend CommandLineProgram. Use -h to print out a list of available programs and their descriptions");
                return 1;
            }
        } catch (final IllegalAccessException e) {
            log.error("Could not access class  " + args[1]);
            return 1;
        } catch (final ClassNotFoundException e) {
            log.error("Could not find command line program with the name " + args[1] +
                    ". Use -h to print out a list of available programs and their descriptions");
            return 1;
        } catch (final NoSuchMethodException e) {
            log.error("Could not find method instanceMain on class " + args[1]);
            return 1;
        } catch (final InstantiationException e) {
            log.error("Could not find instantiate class with the name " + args[1]);
            return 1;
        } catch (final InvocationTargetException e) {
            // Runtime exceptions cause InvocationTargetExceptions with the RuntimeException as the target. We don't want to gobble these up
            // so if we have one be sure to rethrow it.
            if (e.getTargetException() != null) {
                if (RuntimeException.class.isAssignableFrom(e.getTargetException().getClass())) {
                    throw (RuntimeException) e.getTargetException();
                }
            }
            log.error("Could not invoke instanceMain on class " + args[1]);
            return 1;
        }
    }

    public static void main(final String[] args) {
        System.exit(new PicardCommandLine().instanceMain(args));
    }

    private static String extractCommandLineProgram(final String[] args)
            throws InstantiationException, IllegalAccessException {
        if (args.length < 2) {
            printUsage();
            System.exit(1);
        } else {
            if (!args[0].equals("-t") || args[0].equals("-h")) {
                printUsage();
                System.exit(1);
            } else {
                final Set<Class<? extends CommandLineProgram>> allClps = getAllClpClasses();

                for (final Class<? extends CommandLineProgram> programClass : allClps) {
                    if (programClass.getSimpleName().equals(args[1])) {
                        return programClass.getName();
                    }
                }
                printUsage();
                System.exit(1);
            }
        }
        return "";
    }

    private static Set<Class<? extends CommandLineProgram>> getAllClpClasses() {
        Reflections.log = null;
        final Reflections reflections = new Reflections("picard", "edu.mit.broad.picard",
                new FilterBuilder().excludePackage("edu.mit.broad.picard.personal"));
        return reflections.getSubTypesOf(CommandLineProgram.class);
    }

    private static void printUsage() throws IllegalAccessException, InstantiationException {
        final Set<Class<? extends CommandLineProgram>> allClps = getSortedClps();
        final StringBuilder builder = new StringBuilder();
        builder.append("USAGE: PicardCommandLine -t <program name>\n\n");
        builder.append(allClps.size());
        builder.append(" Available Programs: \n");
        builder.append("--------------------------------------------------------------------------------------\n");
        for (final Class<? extends CommandLineProgram> programClass : allClps) {
            if (!Modifier.isAbstract(programClass.getModifiers())) {
                final Reflections reflections = new Reflections(programClass.getName(), new FieldAnnotationsScanner());
                final CommandLineProgram program = programClass.newInstance();
                builder.append(programClass.getSimpleName());
                builder.append("\n\n");
                final Set<Field> usages = reflections.getFieldsAnnotatedWith(Usage.class);
                if (usages.size() == 1) {
                    builder.append(usages.iterator().next().get(program));
                }
                builder.append("\n--------------------------------------------------------------------------------------\n");
            }
        }
        System.err.println(builder.toString());
    }

    private static Set<Class<? extends CommandLineProgram>> getSortedClps() {
        final Set<Class<? extends CommandLineProgram>> allClps = getAllClpClasses();
        final Set<Class<? extends CommandLineProgram>> sortedClps = new TreeSet<Class<? extends CommandLineProgram>>(new Comparator<Class>() {
            @Override
            public int compare(final Class class1, final Class class2) {
                return class1.getSimpleName().compareTo(class2.getSimpleName());
            }
        });
        sortedClps.addAll(allClps);
        return sortedClps;
    }

}
