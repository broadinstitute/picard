package picard.cmdline;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;

import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

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


/**
 * This is the main class of Picard and is the way of executing individual command line programs.
 *
 * CommandLinePrograms are listed in a single command line interface based on the java package specified to instanceMain.
 *
 * If you want your own single command line program, extend this class and give instanceMain a new list of java packages in which to
 * search for classes that extend CommandLineProgram.
 *
 */
public class PicardCommandLine {
    private static final Log log = Log.getInstance(PicardCommandLine.class);
    private static String initializeColor(final String color) {
        if (CommandLineDefaults.COLOR_STATUS) return color;
        else return "";
    }

    /** Provides ANSI colors for the terminal output **/
    private final static String KNRM = initializeColor("\u001B[0m"); // reset
    private final static String KBLD = initializeColor("\u001B[1m"); // Bold
    private final static String KRED = initializeColor("\u001B[31m");
    private final static String KGRN = initializeColor("\u001B[32m");
    private final static String KYEL = initializeColor("\u001B[33m");
    private final static String KBLU = initializeColor("\u001B[34m");
    private final static String KMAG = initializeColor("\u001B[35m");
    private final static String KCYN = initializeColor("\u001B[36m");
    private final static String KWHT = initializeColor("\u001B[37m");
    private final static String KBLDRED = initializeColor("\u001B[1m\u001B[31m");

    /** The name of this unified command line program **/
    private final static String COMMAND_LINE_NAME = PicardCommandLine.class.getSimpleName();

    /** Prefixes for class that annotated by @ExperimentalFeature and @BetaFeature **/
    private final static String BETA_PREFIX = "**BETA - WORK IN PROGRESS** ";
    private final static String EXPERIMENTAL_PREFIX = "**EXPERIMENTAL - USE AT YOUR OWN RISK** ";

    /** The packages we wish to include in our command line **/
    protected static List<String> getPackageList() {
        final List<String> packageList = new ArrayList<String>();
        packageList.add("picard");
        return packageList;
    }

    /**
     * The main method.
     *
     * Give a list of java packages in which to search for classes that extend CommandLineProgram.  Those will be included
     * on the command line.
     * **/
    protected int instanceMain(final String[] args, final List<String> packageList, final String commandLineName) {
        final CommandLineProgram program = extractCommandLineProgram(args, packageList, commandLineName);
        if (null == program) return 1; // no program found!
        // we can lop off the first two arguments but it requires an array copy or alternatively we could update CLP to remove them
        // in the constructor do the former in this implementation.
        final String[] mainArgs = Arrays.copyOfRange(args, 1, args.length);
        return program.instanceMain(mainArgs);
    }

    /** For testing **/
    protected int instanceMain(final String[] args) {
        return instanceMain(args, getPackageList(), COMMAND_LINE_NAME);
    }

    /** Override this if you want to include different java packages to search for classes that extend CommandLineProgram. **/
    public static void main(final String[] args) {
        System.exit(new PicardCommandLine().instanceMain(args, getPackageList(), COMMAND_LINE_NAME));
    }

    /** Returns the command line program specified, or prints the usage and exits with exit code 1 **/
    private static CommandLineProgram extractCommandLineProgram(final String[] args, final List<String> packageList, final String commandLineName) {
        final Map<String, Class<?>> simpleNameToClass = new HashMap<>();
        final List<String> missingAnnotationClasses = new ArrayList<>();
        processAllCommandLinePrograms(
                packageList,
                (Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    // Check for missing annotations
                    if (null == clProperties) {
                        missingAnnotationClasses.add(clazz.getSimpleName());
                    }
                    else if (!clProperties.omitFromCommandLine()) { /** We should check for missing annotations later **/
                        if (simpleNameToClass.containsKey(clazz.getSimpleName())) {
                            throw new RuntimeException("Simple class name collision: " + clazz.getSimpleName());
                        }
                        simpleNameToClass.put(clazz.getSimpleName(), clazz);
                    }
                }
        );
        if (!missingAnnotationClasses.isEmpty()) {
            throw new RuntimeException("The following classes are missing the required CommandLineProgramProperties annotation: " +
                    missingAnnotationClasses.stream().collect(Collectors.joining((", "))));
        }

        final Set<Class<?>> classes = new HashSet<>();
        classes.addAll(simpleNameToClass.values());

        if (args.length < 1) {
            printUsage(classes, commandLineName);
        } else {
            if (args[0].equals("-h")) {
                printUsage(classes, commandLineName);
            } else if(args[0].equals("--list-commands")) {
                printCommandList(classes);
            } else {
                if (simpleNameToClass.containsKey(args[0])) {
                    final Class clazz = simpleNameToClass.get(args[0]);
                    try {
                        return (CommandLineProgram)clazz.newInstance();
                    } catch (final InstantiationException e) {
                        throw new RuntimeException(e);
                    } catch (final IllegalAccessException e) {
                        throw new RuntimeException(e);
                    }
                }
                printUsage(classes, commandLineName);
                printUnknown(classes, args[0]);
            }
        }
        return null;
    }

    /**
     * Process each {@code CommandLineProgram}-derived class given a list of packages.
     * @param packageList list of packages to search
     * @param clpClassProcessor function to process each CommandLineProgram class found in {@code packageList} (note
     *                          that the {@code CommandLineProgramProperties} argument may be null)
     */
    public static void processAllCommandLinePrograms(
            final List<String> packageList,
            final BiConsumer<Class<CommandLineProgram>, CommandLineProgramProperties> clpClassProcessor) {
        final ClassFinder classFinder = new ClassFinder();
        packageList.forEach(pkg -> classFinder.find(pkg, CommandLineProgram.class));

        for (final Class clazz : classFinder.getClasses()) {
            // No interfaces, synthetic, primitive, local, or abstract classes.
            if (!clazz.isInterface() && !clazz.isSynthetic() && !clazz.isPrimitive() && !clazz.isLocalClass()
                    && !Modifier.isAbstract(clazz.getModifiers())) {
                clpClassProcessor.accept(clazz, PicardCommandLine.getProgramProperty(clazz));
            }
        }
    }

    public static CommandLineProgramProperties getProgramProperty(Class clazz) {
        return (CommandLineProgramProperties)clazz.getAnnotation(CommandLineProgramProperties.class);
    }

    private static class SimpleNameComparator implements Comparator<Class> {
        @Override
        public int compare(final Class aClass, final Class bClass) {
            return aClass.getSimpleName().compareTo(bClass.getSimpleName());
        }
    }

    private static void printCommandList(final Set<Class<?>> classes) {
        printUsage(classes, null, true, true);
    }

    private static void printUsage(final Set<Class<?>> classes, final String commandLineName) {
        printUsage(classes, commandLineName, false, false);
    }

    private static void printUsage(final Set<Class<?>> classes, final String commandLineName, boolean commandListOnly, boolean toStdout) {
        final StringBuilder builder = new StringBuilder();
        if (!commandListOnly) {
            builder.append(KBLDRED + "USAGE: " + commandLineName + " " + KGRN + "<program name>" + KBLDRED + " [-h]\n\n" + KNRM);
            builder.append(KBLDRED + "Available Programs:\n" + KNRM);
        }

        /** Group CommandLinePrograms by CommandLineProgramGroup **/
        final Map<Class<? extends CommandLineProgramGroup>, CommandLineProgramGroup> programGroupClassToProgramGroupInstance = new HashMap<Class<? extends CommandLineProgramGroup>, CommandLineProgramGroup>();
        final Map<CommandLineProgramGroup, List<Class>> programsByGroup = new TreeMap<CommandLineProgramGroup, List<Class>>(CommandLineProgramGroup.comparator);
        final Map<Class, CommandLineProgramProperties> programsToProperty = new HashMap<Class, CommandLineProgramProperties>();
        for (final Class clazz : classes) {
            // Get the command line property for this command line program
            final CommandLineProgramProperties property = getProgramProperty(clazz);
            if (null == property) {
                throw new RuntimeException(String.format("The class '%s' is missing the required CommandLineProgramProperties annotation.", clazz.getSimpleName()));
            }
            programsToProperty.put(clazz, property);
            // Get the command line program group for the command line property
            // NB: we want to minimize the number of times we make a new instance, hence programGroupClassToProgramGroupInstance
            CommandLineProgramGroup programGroup = programGroupClassToProgramGroupInstance.get(property.programGroup());
            if (null == programGroup) {
                try {
                    programGroup = property.programGroup().newInstance();
                } catch (final InstantiationException e) {
                    throw new RuntimeException(e);
                } catch (final IllegalAccessException e) {
                    throw new RuntimeException(e);
                }
                programGroupClassToProgramGroupInstance.put(property.programGroup(), programGroup);
            }
            List<Class> programs = programsByGroup.get(programGroup);
            if (null == programs) {
                programsByGroup.put(programGroup, programs = new ArrayList<Class>());
            }
            programs.add(clazz);
        }

        /** Print out the programs in each group **/
        for (final Map.Entry<CommandLineProgramGroup, List<Class>> entry : programsByGroup.entrySet()) {
            final CommandLineProgramGroup programGroup = entry.getKey();

            if (!commandListOnly) {
                builder.append(KWHT + "--------------------------------------------------------------------------------------\n" + KNRM);
                builder.append(String.format("%s%-48s %-45s%s\n", KRED, programGroup.getName() + ":", programGroup.getDescription(), KNRM));
            }

            final List<Class> sortedClasses = new ArrayList<Class>();
            sortedClasses.addAll(entry.getValue());
            Collections.sort(sortedClasses, new SimpleNameComparator());

            for (final Class clazz : sortedClasses) {
                final CommandLineProgramProperties property = programsToProperty.get(clazz);
                if (null == property) {
                    throw new RuntimeException(String.format("Unexpected error: did not find the CommandLineProgramProperties annotation for '%s'", clazz.getSimpleName()));
                }
                if (!commandListOnly) {
                    builder.append(String.format(
                            clazz.getSimpleName().length() >= 45
                                ? "%s    %s    %s%s%s%s%s\n"
                                : "%s    %-45s%s%s%s%s%s\n",
                            KGRN, clazz.getSimpleName(),
                            KRED, getToolSummaryPrefix(clazz),
                            KCYN, property.oneLineSummary(),
                            KNRM));
                }
                else {
                    builder.append(clazz.getSimpleName() + "\n");
                }
            }
            if (!commandListOnly) builder.append(String.format("\n"));
        }
        if (!commandListOnly) builder.append(KWHT + "--------------------------------------------------------------------------------------\n\n" + KNRM);
        if (toStdout) {
            System.out.print(builder.toString());
        }
        else {
            System.err.print(builder.toString());
        }
    }

    private static String getToolSummaryPrefix(Class<?> clazz) {
        if (clazz.getAnnotation(ExperimentalFeature.class) != null) {
            return EXPERIMENTAL_PREFIX;
        }

        if (clazz.getAnnotation(BetaFeature.class) != null) {
            return BETA_PREFIX;
        }

        return "";
    }

    /** similarity floor for matching in printUnknown **/
    private final static int HELP_SIMILARITY_FLOOR = 7;
    private final static int MINIMUM_SUBSTRING_LENGTH = 5;

    /** When a command does not match any known command, searches for similar commands, using the same method as GIT **/
    public static void printUnknown(final Set<Class<?>> classes, final String command) {
        final Map<Class, Integer> distances = new HashMap<Class, Integer>();

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
                distance = StringUtil.levenshteinDistance(command, name, 0, 2, 1, 4);
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
        System.err.println(String.format("'%s' is not a valid command. See PicardCommandLine -h for more information.", command));
        if (bestDistance < HELP_SIMILARITY_FLOOR) {
            System.err.println(String.format("Did you mean %s?", (bestN < 2) ? "this" : "one of these"));
            for (final Class clazz : classes) {
                if (bestDistance == distances.get(clazz)) {
                    System.err.println(String.format("        %s", clazz.getSimpleName()));
                }
            }
        }
    }
}
