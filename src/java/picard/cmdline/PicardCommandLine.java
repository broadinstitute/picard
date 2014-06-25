package picard.cmdline;

import htsjdk.samtools.util.Log;
import picard.sam.SamToFastq;

import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.ServiceLoader;
import java.util.Set;
import java.util.TreeMap;
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

    /** Provides ANSI colors for the terminal output **/
    private final static String KNRM = "\u001B[0m"; // reset
    private final static String KBLD = "\u001B[1m"; // Bold
    private final static String KRED = "\u001B[31m";
    private final static String KGRN = "\u001B[32m";
    private final static String KYEL = "\u001B[33m";
    private final static String KBLU = "\u001B[34m";
    private final static String KMAG = "\u001B[35m";
    private final static String KCYN = "\u001B[36m";
    private final static String KWHT = "\u001B[37m";
    private final static String KBLDRED = "\u001B[1m\u001B[31m";

    /** similarity floor for matching in help **/
    private final static int HELP_SIMILARITY_FLOOR = 7;
    private final static int MINIMUM_SUBSTRING_LENGTH = 5;

    /** Override this if you wish to search for command line programs in different packages **/
    protected static String[] packageList = {"picard"};

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
        /** Get the set of classes that are our command line programs **/
        final ClassFinder classFinder = new ClassFinder();
        for (final String pkg : packageList) {
            classFinder.find(pkg, CommandLineProgram.class);
        }
        final Set<Class<?>> classes = new HashSet<Class<?>>();
        for (final Class clazz : classFinder.getClasses()) {
            // No interfaces, synthetic, primitive, local, or abstract classes.
            if (!clazz.isInterface() && !clazz.isSynthetic() && !clazz.isPrimitive() && !clazz.isLocalClass() && !Modifier.isAbstract(clazz.getModifiers())) {
                classes.add(clazz);
            }
        }

        if (args.length < 2) {
            printUsage(classes);
            System.exit(1);
        } else {
            if (!args[0].equals("-T") || args[0].equals("-h")) {
                printUsage(classes);
                System.exit(1);
            } else {
                for (final Class clazz : classes) {
                    if (clazz.getSimpleName().equals(args[1])) {
                        try {
                            return (CommandLineProgram)clazz.newInstance();
                        } catch (InstantiationException e) {
                            throw new RuntimeException(e);
                        } catch (IllegalAccessException e) {
                            throw new RuntimeException(e);
                        }
                    }
                }
                printUsage(classes);
                printUnknown(classes, args[1]);
                System.exit(1);
            }
        }
        return null;
    }

    private static void printUsage(final Set<Class<?>> classes) {
        final StringBuilder builder = new StringBuilder();
        builder.append(KBLDRED + "USAGE: PicardCommandLine -T " + KGRN + "<program name>" + KBLDRED + " [-h]\n\n" + KNRM);
        builder.append(KBLDRED + "Available Programs:\n" + KNRM);

        // Maps a subclass of CommandLineProgramGroup to a list of CommandLineProgramProperties annotations that have that subclass as their programGroup
        final Map<Class, List<CommandLineProgramProperties>> classToPropertiesMap = new TreeMap<Class, List<CommandLineProgramProperties>>(CommandLineProgramGroup.comparator);
        // Maps a CommandLineProgramProperties property to a specific CommandLineProgramGroup subclass
        final Map<CommandLineProgramProperties, Class> propertiesToClassMap = new HashMap<CommandLineProgramProperties, Class>();
        for (final Class clazz : classes) {
            final CommandLineProgramProperties property = (CommandLineProgramProperties)clazz.getAnnotation(CommandLineProgramProperties.class);
            if (null == property) {
                throw new RuntimeException("No CommandLineProgramProperties annotation for: " + clazz.getSimpleName());
            }
            List<CommandLineProgramProperties> programs = classToPropertiesMap.get(property.programGroup());
            if (null == programs) {
                classToPropertiesMap.put(property.programGroup(), programs = new ArrayList<CommandLineProgramProperties>());
            }
            programs.add(property);
            propertiesToClassMap.put(property, clazz);
        }

        // Print out the programs in each group
        for (final Map.Entry<Class, List<CommandLineProgramProperties>> entry : classToPropertiesMap.entrySet()) {
            Class propertyClass = entry.getKey();
            CommandLineProgramGroup programGroup;
            try {
                programGroup = (CommandLineProgramGroup)propertyClass.newInstance();
            } catch (InstantiationException e) {
                throw new RuntimeException(e);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
            builder.append("--------------------------------------------------------------------------------------\n");
            builder.append(String.format("%s%-48s %-45s%s\n", KRED, programGroup.getName() + ":", programGroup.getDescription(), KNRM));
            for (final CommandLineProgramProperties property : entry.getValue()) {
                builder.append(String.format("%s    %-45s%s%s%s\n", KGRN, propertiesToClassMap.get(property).getSimpleName(), KCYN, property.usageShort(), KNRM));
            }
            builder.append(String.format("\n"));
        }
        builder.append(KWHT + "--------------------------------------------------------------------------------------\n" + KNRM);
        System.err.println(builder.toString());
    }

    public static void printUnknown(final Set<Class<?>> classes, final String command) {
        Map<Class, Integer> distances = new HashMap<Class, Integer>();

        int bestDistance = Integer.MAX_VALUE;
        int bestN = 0;

        // Score against all classes
        for (final Class clazz : classes) {
            final String name = clazz.getSimpleName();
            final int distance;
            if (name.equals(command)) {
                throw new RuntimeException("Command matches: " + command);
            }
            if (name.startsWith(command) || (MINIMUM_SUBSTRING_LENGTH <= command.length() && name.contains(command))) {
                distance = 0;
            }
            else {
                distance = levenshteinDistance(command, name, 0, 2, 1, 4);
            }
            distances.put(clazz, distance);

            if (distance < bestDistance) {
                bestDistance = distance;
                bestN = 1;
            }
            else if (distance == bestDistance) {
                bestN++;
            }
        }

        // Upper bound on the similarity score
        if (0 == bestDistance && bestN == classes.size()) {
            bestDistance = HELP_SIMILARITY_FLOOR + 1;
        }

        // Output similar matches
        // TODO
        System.err.println(String.format("'%s' is not a valid command. See PicardCommandLine --help for more information.", command));
        if (bestDistance < HELP_SIMILARITY_FLOOR) {
            System.err.println(String.format("Did you mean %s?", (bestN < 2) ? "this" : "one of these"));
            for (final Class clazz : classes) {
                if (bestDistance == distances.get(clazz)) {
                    System.err.println(String.format("        %s", clazz.getSimpleName()));
                }
            }
        }
    }

    /*
     * This s from GIT!
     *  This function implements the Damerau-Levenshtein algorithm to
     * calculate a distance between strings.
     *
     * Basically, it says how many letters need to be swapped, substituted,
     * deleted from, or added to string1, at least, to get string2.
     *
     * The idea is to build a distance matrix for the substrings of both
     * strings.  To avoid a large space complexity, only the last three rows
     * are kept in memory (if swaps had the same or higher cost as one deletion
     * plus one insertion, only two rows would be needed).
     *
     * At any stage, "i + 1" denotes the length of the current substring of
     * string1 that the distance is calculated for.
     *
     * row2 holds the current row, row1 the previous row (i.e. for the substring
     * of string1 of length "i"), and row0 the row before that.
     *
     * In other words, at the start of the big loop, row2[j + 1] contains the
     * Damerau-Levenshtein distance between the substring of string1 of length
     * "i" and the substring of string2 of length "j + 1".
     *
     * All the big loop does is determine the partial minimum-cost paths.
     *
     * It does so by calculating the costs of the path ending in characters
     * i (in string1) and j (in string2), respectively, given that the last
     * operation is a substitution, a swap, a deletion, or an insertion.
     *
     * This implementation allows the costs to be weighted:
     *
     * Note that this algorithm calculates a distance _iff_ d == a.
     */
    private static int levenshteinDistance(final String string1, final String string2, int swap, int substitution, int insertion, int deletion) {
        int i, j;

        int[] row0 = new int[(string2.length() + 1)];
        int[] row1 = new int[(string2.length() + 1)];
        int[] row2 = new int[(string2.length() + 1)];
        int[] dummy;


        byte[] str1 = string1.getBytes();
        byte[] str2 = string2.getBytes();

        for (j = 0; j < str2.length; j++) {
            row1[j] = j * insertion;
        }
        for (i = 0; i < str1.length; i++) {
            row2[0] = (i + 1) * deletion;
            for (j = 0; j < str2.length; j++) {
                /* substitution */
                row2[j + 1] = row1[j];
                if (str1[i] != str2[j]) {
                    row2[j + 1] += substitution;
                }
                /* swap */
                if (i > 0 && j > 0 && str1[i - 1] == str2[j] &&
                        str1[i] == str2[j - 1] &&
                        row2[j + 1] > row0[j - 1] + swap) {
                    row2[j + 1] = row0[j - 1] + swap;
                }
                /* deletion */
                if (row2[j + 1] > row1[j + 1] + deletion) {
                    row2[j + 1] = row1[j + 1] + deletion;
                }
                /* insertion */
                if (row2[j + 1] > row2[j] + insertion) {
                    row2[j + 1] = row2[j] + insertion;
                }
            }

            dummy = row0;
            row0 = row1;
            row1 = row2;
            row2 = dummy;
        }

        i = row1[str2.length];

        return i;
    }
}
