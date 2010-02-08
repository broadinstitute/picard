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
package net.sf.picard.cmdline;

import java.io.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.*;

import net.sf.samtools.util.StringUtil;
import net.sf.samtools.util.CloserUtil;
import net.sf.picard.PicardException;

/**
 * Annotation-driven utility for parsing command-line arguments, checking for errors, and producing usage message.
 *
 * This class supports options of the form KEY=VALUE, plus positional arguments.  Positional arguments must not contain
 * an equal sign lest they be mistaken for a KEY=VALUE pair.
 *
 * The caller must supply an object that both defines the command line and has the parsed options set into it.
 * For each possible KEY=VALUE option, there must be a public data member annotated with @Option.  The KEY name is
 * the name of the data member.  An abbreviated name may also be specified with the shortName attribute of @Option.
 * If the data member is a List<T>, then the option may be specified multiple times.  The type of the data member,
 * or the type of the List element must either have a ctor T(String), or must be an Enum.  List options must
 * be initialized by the caller with some kind of list.  Any other option that is non-null is assumed to have the given
 * value as a default.  If an option has no default value, and does not have the optional attribute of @Option set,
 * is required.  For List options, minimum and maximum number of elements may be specified in the @Option annotation.
 *
 * A single List data member may be annotated with the @PositionalArguments.  This behaves similarly to a Option
 * with List data member: the caller must initialize the data member, the type must be constructable from String, and
 * min and max number of elements may be specified.  If no @PositionalArguments annotation appears in the object,
 * then it is an error for the command line to contain positional arguments.
 *
 * A single String public data member may be annotated with @Usage.  This string, if present, is used to
 * construct the usage message.  Details about the possible options are automatically appended to this string.
 * If @Usage does not appear, a boilerplate usage message is used.
 */
public class CommandLineParser {
    // For formatting option section of usage message.
    private static final int OPTION_COLUMN_WIDTH = 30;
    private static final int DESCRIPTION_COLUMN_WIDTH = 90;

    private static final Boolean[] TRUE_FALSE_VALUES = {Boolean.TRUE, Boolean.FALSE};

    // Use these if no @Usage annotation
    private static final String defaultUsagePreamble = "Usage: program [options...]\n";
    private static final String defaultUsagePreambleWithPositionalArguments =
            "Usage: program [options...] [positional-arguments...]\n";
    private static final String OPTIONS_FILE = "OPTIONS_FILE";

    /**
     * A typical command line program will call this to get the beginning of the usage message,
     * and then append a description of the program, like this:
     *
     * \@Usage(programVersion=PROGRAM_VERSION)
     * public String USAGE = CommandLineParser.getStandardUsagePreamble(getClass()) + "Frobnicates the freebozzle."
     */
    public static String getStandardUsagePreamble(final Class mainClass) {
        return "USAGE: " + mainClass.getName() + " [options]\n\n";
    }

    // This is the object that the caller has provided that contains annotations,
    // and into which the values will be assigned.
    private final Object callerOptions;

    private String usagePreamble;
    // null if no @PositionalArguments annotation
    private Field positionalArguments;
    private int minPositionalArguments;
    private int maxPositionalArguments;

    // List of all the data members with @Option annotation
    private final List<OptionDefinition> optionDefinitions = new ArrayList<OptionDefinition>();

    // Maps long name, and short name, if present, to an option definition that is
    // also in the optionDefinitions list.
    private final Map<String, OptionDefinition> optionMap = new HashMap<String, OptionDefinition>();

    // For printing error messages when parsing command line.
    private PrintStream messageStream;

    // In case implementation wants to get at arg for some reason.
    private String[] argv;

    private String programVersion = "";

    // The command line used to launch this program, including non-null default options that
    // weren't explicitly specified. This is used for logging and debugging.
    private String commandLine;

    /**
     * This attribute is here just to facilitate printing usage for OPTIONS_FILE
     */
    public File IGNORE_THIS_PROPERTY;

    /**
     * Prepare for parsing command line arguments, by validating annotations.
     * @param callerOptions This object contains annotations that define the acceptable command-line options,
     * and ultimately will receive the settings when a command line is parsed.
     */
    public CommandLineParser(final Object callerOptions) {
        this.callerOptions = callerOptions;

        for (final Field field : this.callerOptions.getClass().getFields()) {
            if (field.getAnnotation(PositionalArguments.class) != null) {
                handlePositionalArgumentAnnotation(field);
            }
            if (field.getAnnotation(Usage.class) != null) {
                handleUsageAnnotation(field);
            }
            if (field.getAnnotation(Option.class) != null) {
                handleOptionAnnotation(field);
            }
        }

        if (usagePreamble == null) {
            if (positionalArguments == null) {
                usagePreamble = defaultUsagePreamble;
            } else {
                usagePreamble = defaultUsagePreambleWithPositionalArguments;
            }
        }
    }

    /**
     * Print a usage message based on the options object passed to the ctor.
     * @param stream Where to write the usage message.
     */
    public void usage(final PrintStream stream) {
        stream.print(usagePreamble);
        if (!optionDefinitions.isEmpty()) {
            stream.println("\nOptions:\n");
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                printOptionUsage(stream, optionDefinition);
            }
        }
        final Field fileField;
        try {
            fileField = getClass().getField("IGNORE_THIS_PROPERTY");
        } catch (NoSuchFieldException e) {
            throw new PicardException("Should never happen", e);
        }
        final OptionDefinition optionsFileOptionDefinition =
                new OptionDefinition(fileField, OPTIONS_FILE, "",
                        "File of OPTION_NAME=value pairs.  No positional parameters allowed.  Unlike command-line options, " +
                "unrecognized options are ignored.  " + "A single-valued option set in an options file may be overridden " +
                "by a subsequent command-line option.  " +
                "A line starting with '#' is considered a comment.", false, true, 0, Integer.MAX_VALUE, null, new String[0]);
        printOptionUsage(stream, optionsFileOptionDefinition);
    }

    public void htmlUsage(final PrintStream stream, final String programName) {
        stream.println("<a name=\"" + programName + "\"/>");
        stream.println("<h3>" + programName + "</h3>");
        stream.println("<p>" + usagePreamble + "</p>");
        stream.println("<table>");
        stream.println("<tr><th>Option</th><th>Description</th></tr>");
        for (final OptionDefinition optionDefinition : optionDefinitions) {
            printHtmlOptionUsage(stream, optionDefinition);
        }
        stream.println("</table>");
        stream.println("<br/>");
    }

    /**
     * Parse command-line options, and store values in callerOptions object passed to ctor.
     * @param messageStream Where to write error messages.
     * @param args Command line tokens.
     * @return true if command line is valid.
     */
    public boolean parseOptions(final PrintStream messageStream, final String[] args) {
        this.argv = args;
        this.messageStream = messageStream;
        for (int i = 0; i < args.length; ++i) {
            final String arg = args[i];
            if (arg.equals("-h") || arg.equals("--help")) {
                usage(messageStream);
                return false;
            }
            final String[] pair = arg.split("=", 2);
            if (pair.length == 2 && pair[1].length() == 0) {

                if (i < args.length - 1) {
                    pair[1] = args[++i];
                }
            }
            if (pair.length == 2) {
                if (pair[0].equals(OPTIONS_FILE)) {
                    if (!parseOptionsFile(pair[1])) {
                        messageStream.println();
                        usage(messageStream);
                        return false;
                    }
                } else {
                    if (!parseOption(pair[0], pair[1], false)) {
                        messageStream.println();
                        usage(messageStream);
                        return false;
                    }
                }
            } else if (!parsePositionalArgument(arg)) {
                messageStream.println();
                usage(messageStream);
                return false;
            }
        }
        if (!checkNumArguments()) {
            messageStream.println();
            usage(messageStream);
            return false;
        }

        return true;
    }

    /**
     * After command line has been parsed, make sure that all required options have values, and that
     * lists with minimum # of elements have sufficient.
     * @return true if valid
     */
    private boolean checkNumArguments() {
        //Also, since we're iterating over all options and args, use this opportunity to recreate the commandLineString
        final StringBuffer commandLineString = new StringBuffer();
        commandLineString.append( callerOptions.getClass().getName() );
        try {
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                final StringBuilder mutextOptionNames = new StringBuilder();
                for (final String mutexOption : optionDefinition.mutuallyExclusive) {
                    final OptionDefinition mutextOptionDef = optionMap.get(mutexOption);
                    if (mutextOptionDef != null && mutextOptionDef.hasBeenSet) {
                        mutextOptionNames.append(" ").append(mutextOptionDef.name);
                    }
                }
                if (optionDefinition.hasBeenSet && mutextOptionNames.length() > 0) {
                    messageStream.println("ERROR: Option '" + optionDefinition.name +
                            "' cannot be used in conjunction with option(s)" +
                            mutextOptionNames.toString());
                    return false;
                }
                if (optionDefinition.isCollection) {
                    final Collection c = (Collection)optionDefinition.field.get(callerOptions);
                    if (c.size() < optionDefinition.minElements) {
                        messageStream.println("ERROR: Option '" + optionDefinition.name + "' must be specified at least " +
                        optionDefinition.minElements + " times.");
                        return false;
                    }
                } else if (!optionDefinition.optional && !optionDefinition.hasBeenSet && mutextOptionNames.length() == 0) {
                    messageStream.print("ERROR: Option '" + optionDefinition.name + "' is required");
                    if (optionDefinition.mutuallyExclusive.isEmpty()) {
                        messageStream.println(".");
                    } else {
                        messageStream.println(" unless any of " + optionDefinition.mutuallyExclusive + " are specified.");
                    }
                    return false;
                }

            }
            if (positionalArguments != null) {
                final Collection c = (Collection)positionalArguments.get(callerOptions);
                if (c.size() < minPositionalArguments) {
                    messageStream.println("ERROR: At least " + minPositionalArguments +
                            " positional arguments must be specified.");
                    return false;
                }
                for( Object posArg : c ) {
                    commandLineString.append(" " + posArg.toString());
                }
            }
            //first, append args that were explicitly set
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                if(optionDefinition.hasBeenSet) {
                    commandLineString.append(" " + optionDefinition.name + "=" + optionDefinition.field.get(callerOptions));
                }
            }
            commandLineString.append("   "); //separator to tell the 2 apart
            //next, append args that weren't explicitly set, but have a default value
            for (final OptionDefinition optionDefinition : optionDefinitions) {
                if(!optionDefinition.hasBeenSet && !optionDefinition.defaultValue.equals("null")) {
                    commandLineString.append(" " + optionDefinition.name + "=" + optionDefinition.defaultValue);
                }
            }
            this.commandLine = commandLineString.toString();
            return true;
        } catch (IllegalAccessException e) {
            // Should never happen because lack of publicness has already been checked.
            throw new RuntimeException(e);
        }


    }

    private boolean parsePositionalArgument(final String stringValue) {
        if (positionalArguments == null) {
            messageStream.println("ERROR: Invalid argument '" + stringValue + "'.");
            return false;
        }
        final Object value;
        try {
            value = constructFromString(getUnderlyingType(positionalArguments), stringValue);
        } catch (CommandLineParseException e) {
            messageStream.println("ERROR: " + e.getMessage());
            return false;
        }
        final Collection c;
        try {
            c = (Collection)positionalArguments.get(callerOptions);
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
        if (c.size() >= maxPositionalArguments) {
            messageStream.println("ERROR: No more than " + maxPositionalArguments +
                    " positional arguments may be specified on the command line.");
            return false;
        }
        c.add(value);
        return true;
    }

    private boolean parseOption(String key, final String stringValue, final boolean optionsFile) {
        key = key.toUpperCase();
        final OptionDefinition optionDefinition = optionMap.get(key);
        if (optionDefinition == null) {
            if (optionsFile) {
                // Silently ignore unrecognized option from options file
                return true;
            }
            messageStream.println("ERROR: Unrecognized option: " + key);
            return false;
        }
        if (!optionDefinition.isCollection) {
            if (optionDefinition.hasBeenSet && !optionDefinition.hasBeenSetFromOptionsFile) {
                messageStream.println("ERROR: Option '" + key + "' cannot be specified more than once.");
                return false;
            }
        }
        final Object value;
        try {
            if(stringValue.equals("null")) {
                //"null" is a special value that allows the user to override any default
                //value set for this arg. It can only be used for optional args. When
                //used for a list arg, it will clear the list.
                if(optionDefinition.optional) {
                    value = null;
                } else {
                    messageStream.println("ERROR: non-null value must be provided for '" + key + "'.");
                    return false;
                }
            } else {
                value = constructFromString(getUnderlyingType(optionDefinition.field), stringValue);
            }

        } catch (CommandLineParseException e) {
            messageStream.println("ERROR: " + e.getMessage());
            return false;
        }
        try {
            if (optionDefinition.isCollection) {
                final Collection c = (Collection)optionDefinition.field.get(callerOptions);
                if(value == null) {
                    //user specified this arg=null which is interpreted as empty list
                    c.clear();
                } else if (c.size() >= optionDefinition.maxElements) {
                    messageStream.println("ERROR: Option '" + key + "' cannot be used more than " +
                            optionDefinition.maxElements + " times.");
                    return false;
                } else {
                    c.add(value);
                }
            } else {
                optionDefinition.field.set(callerOptions, value);
                optionDefinition.hasBeenSet = true;
                optionDefinition.hasBeenSetFromOptionsFile = optionsFile;
            }
        } catch (IllegalAccessException e) {
            // Should never happen because we only iterate through public fields.
            throw new RuntimeException(e);
        }
        return true;
    }

    /**
     * Parsing of options from file is looser than normal.  Any unrecognized options are
     * ignored, and a single-valued option that is set in a file may be overridden by a
     * subsequent appearance of that option.
     * A line that starts with '#' is ignored.
     * @param optionsFile
     * @return false if a fatal error occurred
     */
    private boolean parseOptionsFile(final String optionsFile) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(optionsFile));
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                final String[] pair = line.split("=", 2);
                if (pair.length == 2) {
                    if (!parseOption(pair[0], pair[1], true)) {
                        messageStream.println();
                        usage(messageStream);
                        return false;
                    }
                } else {
                    messageStream.println("Strange line in OPTIONS_FILE " + optionsFile + ": " + line);
                    usage(messageStream);
                    return false;
                }
            }
            reader.close();
            return true;

        } catch (IOException e) {
            throw new PicardException("I/O error loading OPTIONS_FILE=" + optionsFile, e);
        } finally {
            CloserUtil.close(reader);
        }
    }

    private void printHtmlOptionUsage(final PrintStream stream, final OptionDefinition optionDefinition) {
        final String type = getUnderlyingType(optionDefinition.field).getSimpleName();
        String optionLabel = optionDefinition.name + "=" + type;
        stream.println("<tr><td>" + optionLabel + "</td><td>" +
                makeOptionDescription(optionDefinition) + "</td></tr>");
    }

    private void printOptionUsage(final PrintStream stream, final OptionDefinition optionDefinition) {
        final String type = getUnderlyingType(optionDefinition.field).getSimpleName();
        String optionLabel = optionDefinition.name + "=" + type;
        stream.print(optionLabel);
        if (optionDefinition.shortName.length() > 0) {
            stream.println();
        }
        if (optionDefinition.shortName.length() > 0) {
            optionLabel = optionDefinition.shortName + "=" + type;
            stream.print(optionLabel);
        }
        int numSpaces = OPTION_COLUMN_WIDTH - optionLabel.length();
        if (optionLabel.length() > OPTION_COLUMN_WIDTH) {
            stream.println();
            numSpaces = OPTION_COLUMN_WIDTH;
        }
        printSpaces(stream, numSpaces);
        final String optionDescription = makeOptionDescription(optionDefinition);
        final String wrappedDescription = StringUtil.wordWrap(optionDescription, DESCRIPTION_COLUMN_WIDTH);
        final String[] descriptionLines = wrappedDescription.split("\n");
        for (int i = 0; i < descriptionLines.length; ++i) {
            if (i > 0) {
                printSpaces(stream, OPTION_COLUMN_WIDTH);
            }
            stream.println(descriptionLines[i]);
        }
        stream.println();
    }

    private String makeOptionDescription(final OptionDefinition optionDefinition) {
        final StringBuilder sb = new StringBuilder();
        if (optionDefinition.doc.length() > 0) {
            sb.append(optionDefinition.doc);
            sb.append(" ");
        }
        if (optionDefinition.optional && !optionDefinition.isCollection) {
            sb.append("Default value: ");
            sb.append(optionDefinition.defaultValue);
            sb.append(". ");
            if(!optionDefinition.defaultValue.equals("null")) {
                sb.append("This option can be set to 'null' to clear the default value. ");
            }
        } else if (!optionDefinition.isCollection){
            sb.append("Required. ");
        }
        Object[] enumConstants = getUnderlyingType(optionDefinition.field).getEnumConstants();
        if (enumConstants == null && getUnderlyingType(optionDefinition.field) == Boolean.class) {
            enumConstants = TRUE_FALSE_VALUES;
        }
        if (enumConstants != null) {
            sb.append("Possible values: {");
            for (int i = 0; i < enumConstants.length; ++i) {
                if (i > 0) {
                    sb.append(", ");
                }
                sb.append(enumConstants[i].toString());
            }
            sb.append("} ");
        }
        if (optionDefinition.isCollection) {
            if (optionDefinition.minElements == 0) {
                if (optionDefinition.maxElements == Integer.MAX_VALUE) {
                    sb.append("This option may be specified 0 or more times. ");
                } else {
                    sb.append("This option must be specified no more than " + optionDefinition.maxElements + " times. ");
                }
            } else  if (optionDefinition.maxElements == Integer.MAX_VALUE) {
                sb.append("This option must be specified at least " + optionDefinition.minElements + " times. ");
            } else {
                sb.append("This option may be specified between " + optionDefinition.minElements +
                " and " + optionDefinition.maxElements + " times. ");
            }

            if(!optionDefinition.defaultValue.equals("null")) {
                sb.append("This option can be set to 'null' to clear the default list. ");
            }

        }
        if (!optionDefinition.mutuallyExclusive.isEmpty()) {
            sb.append(" Cannot be used in conjuction with option(s)");
            for (final String option : optionDefinition.mutuallyExclusive) {
                final OptionDefinition mutextOptionDefinition = optionMap.get(option);
                sb.append(" ").append(mutextOptionDefinition.name);
                if (mutextOptionDefinition.shortName.length() > 0) {
                    sb.append(" (").append(mutextOptionDefinition.shortName).append(")");
                }
            }
        }
        return sb.toString();
    }

    private void printSpaces(final PrintStream stream, final int numSpaces) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numSpaces; ++i) {
            sb.append(" ");
        }
        stream.print(sb);
    }

    private void handleOptionAnnotation(final Field field) {
        try {
            final Option optionAnnotation = field.getAnnotation(Option.class);
            final boolean isCollection = isCollectionField(field);
            if (isCollection) {
                if (optionAnnotation.maxElements() == 0) {
                    throw new CommandLineParserDefinitionException("@Option member " + field.getName() +
                    "has maxElements = 0");
                }
                if (optionAnnotation.minElements() > optionAnnotation.maxElements()) {
                    throw new CommandLineParserDefinitionException("In @Option member " + field.getName() +
                            ", minElements cannot be > maxElements");
                }
                if (field.get(callerOptions) == null) {
                    createCollection(field, callerOptions, "@Option");
                }
            }
            if (!canBeMadeFromString(getUnderlyingType(field))) {
                throw new CommandLineParserDefinitionException("@Option member " + field.getName() +
                " must have a String ctor or be an enum");
            }

            final OptionDefinition optionDefinition = new OptionDefinition(field,
                    field.getName(),
                    optionAnnotation.shortName(),
                    optionAnnotation.doc(), optionAnnotation.optional() || (field.get(callerOptions) != null),
                    isCollection, optionAnnotation.minElements(),
                    optionAnnotation.maxElements(), field.get(callerOptions),
                    optionAnnotation.mutex());

            for (final String option : optionAnnotation.mutex()) {
                final OptionDefinition mutextOptionDef = optionMap.get(option);
                if (mutextOptionDef != null) {
                    mutextOptionDef.mutuallyExclusive.add(field.getName());
                }
            }
            if (optionMap.containsKey(optionDefinition.name)) {
                throw new CommandLineParserDefinitionException(optionDefinition.name + " has already been used");
            }
            optionMap.put(optionDefinition.name, optionDefinition);
            if (optionDefinition.shortName.length() > 0) {
                if (optionMap.containsKey(optionDefinition.shortName)) {
                    throw new CommandLineParserDefinitionException(optionDefinition.shortName + " has already been used");
                }
                optionMap.put(optionDefinition.shortName, optionDefinition);
            }
            optionDefinitions.add(optionDefinition);
        } catch (IllegalAccessException e) {
            throw new CommandLineParserDefinitionException(field.getName() +
                    " must have public visibility to have @Option annotation");
        }
    }

    private void handleUsageAnnotation(final Field field) {
        if (usagePreamble != null) {
            throw new CommandLineParserDefinitionException
                    ("@Usage cannot be used more than once in an option class.");
        }
        try {
            usagePreamble = (String)field.get(callerOptions);
            final Usage usageAnnotation = field.getAnnotation(Usage.class);
            if (usageAnnotation.programVersion().length() > 0) {
                this.programVersion = usageAnnotation.programVersion();
                usagePreamble += "Version: " + usageAnnotation.programVersion() + "\n";
            }
        } catch (IllegalAccessException e) {
            throw new CommandLineParserDefinitionException("@Usage data member must be public");
        } catch (ClassCastException e) {
            throw new CommandLineParserDefinitionException
                    ("@Usage can only be applied to a String data member.");
        }
    }

    private void handlePositionalArgumentAnnotation(final Field field) {
        if (positionalArguments != null) {
            throw new CommandLineParserDefinitionException
                    ("@PositionalArguments cannot be used more than once in an option class.");
        }
        positionalArguments = field;
        if (!isCollectionField(field)) {
            throw new CommandLineParserDefinitionException("@PositionalArguments must be applied to a Collection");
        }

        if (!canBeMadeFromString(getUnderlyingType(field))) {
            throw new CommandLineParserDefinitionException("@PositionalParameters member " + field.getName() +
            "does not have a String ctor");
        }

        final PositionalArguments positionalArgumentsAnnotation = field.getAnnotation(PositionalArguments.class);
        minPositionalArguments = positionalArgumentsAnnotation.minElements();
        maxPositionalArguments = positionalArgumentsAnnotation.maxElements();
        if (minPositionalArguments > maxPositionalArguments) {
            throw new CommandLineParserDefinitionException("In @PositionalArguments, minElements cannot be > maxElements");
        }
        try {
            if (field.get(callerOptions) == null) {
                createCollection(field, callerOptions, "@PositionalParameters");
            }
        } catch (IllegalAccessException e) {
            throw new CommandLineParserDefinitionException(field.getName() +
                    " must have public visibility to have @PositionalParameters annotation");

        }
    }

    private boolean isCollectionField(final Field field) {
        try {
            field.getType().asSubclass(Collection.class);
            return true;
        } catch (ClassCastException e) {
            return false;
        }
    }

    private void createCollection(final Field field, final Object callerOptions, String annotationType) throws IllegalAccessException {
        try {
            field.set(callerOptions, field.getType().newInstance());
        } catch (Exception ex) {
            try {
                field.set(callerOptions, new ArrayList());
            } catch (IllegalArgumentException e) {
                throw new CommandLineParserDefinitionException("In collection " + annotationType + " member " + field.getName() +
                        " cannot be constructed or auto-initialized with ArrayList, so collection must be initialized explicitly.");
            }

        }

    }

    /**
     * Returns the type that each instance of the argument needs to be converted to. In
     * the case of primitive fields it will return the wrapper type so that String
     * constructors can be found.
     */
    private Class getUnderlyingType(final Field field) {
        if (isCollectionField(field)) {
            final ParameterizedType clazz = (ParameterizedType)(field.getGenericType());
            final Type[] genericTypes = clazz.getActualTypeArguments();
            if (genericTypes.length != 1) {
                throw new CommandLineParserDefinitionException("Strange collection type for field " + field.getName());
            }
            return (Class)genericTypes[0];

        }
        else {
            final Class type = field.getType();
            if (type == Byte.TYPE)    return Byte.class;
            if (type == Short.TYPE)   return Short.class;
            if (type == Integer.TYPE) return Integer.class;
            if (type == Long.TYPE)    return Long.class;
            if (type == Float.TYPE)   return Float.class;
            if (type == Double.TYPE)  return Double.class;
            if (type == Boolean.TYPE) return Boolean.class;

            return type;
        }
    }

    // True if clazz is an enum, or if it has a ctor that takes a single String argument.
    private boolean canBeMadeFromString(final Class clazz) {
        if (clazz.isEnum()) {
            return true;
        }
        try {
            clazz.getConstructor(String.class);
            return true;
        } catch (NoSuchMethodException e) {
            return false;
        }
    }

    private Object constructFromString(final Class clazz, final String s) {
        try {
            if (clazz.isEnum()) {
                try {
                    return Enum.valueOf(clazz, s);
                } catch (IllegalArgumentException e) {
                    throw new CommandLineParseException("'" + s + "' is not a valid value for " +
                            clazz.getSimpleName() + ".", e);
                }
            }
            final Constructor ctor = clazz.getConstructor(String.class);
            return ctor.newInstance(s);
        } catch (NoSuchMethodException e) {
            // Shouldn't happen because we've checked for presence of ctor
            throw new CommandLineParseException("Cannot find string ctor for " + clazz.getName(), e);
        } catch (InstantiationException e) {
            throw new CommandLineParseException("Abstract class '" + clazz.getSimpleName() +
                    "'cannot be used for an option value type.", e);
        } catch (IllegalAccessException e) {
            throw new CommandLineParseException("String constructor for option value type '" + clazz.getSimpleName() +
                    "' must be public.", e);
        } catch (InvocationTargetException e) {
            throw new CommandLineParseException("Problem constructing " + clazz.getSimpleName() + " from the string '" + s + "'.",
                    e.getCause());
        }
    }

    public String[] getArgv() {
        return argv;
    }

    private static class OptionDefinition {
        final Field field;
        final String name;
        final String shortName;
        final String doc;
        final boolean optional;
        final boolean isCollection;
        final int minElements;
        final int maxElements;
        final String defaultValue;
        boolean hasBeenSet = false;
        boolean hasBeenSetFromOptionsFile = false;
        Set<String> mutuallyExclusive;

        private OptionDefinition(final Field field, final String name, final String shortName, final String doc, final boolean optional, final boolean collection,
                                 final int minElements, final int maxElements, final Object defaultValue, final String[] mutuallyExclusive) {
            this.field = field;
            this.name = name.toUpperCase();
            this.shortName = shortName.toUpperCase();
            this.doc = doc;
            this.optional = optional;
            isCollection = collection;
            this.minElements = minElements;
            this.maxElements = maxElements;
            if (defaultValue != null) {
                if( isCollection && ((Collection) defaultValue).isEmpty()) {
                    //treat empty collections the same as uninitialized primitive types
                    this.defaultValue = "null";
                } else {
                    //this is an intialized primitive type or a non-empty collection
                    this.defaultValue = defaultValue.toString();
                }
            } else {
                this.defaultValue = "null";
            }
            this.mutuallyExclusive = new HashSet<String>(Arrays.asList(mutuallyExclusive));
        }
    }

    public String getProgramVersion() { return programVersion; }

    /**
     * The commandline used to run this program, including any default args that
     * weren't necessarily specified. This is used for logging and debugging.
     *
     * NOTE: {@link #parseOptions(PrintStream, String[])} must be called before
     * calling this method.
     *
     * @return The commandline, or null if {@link #parseOptions(PrintStream, String[])}
     * hasn't yet been called, or didn't complete successfully.
     */
    public String getCommandLine() { return commandLine; }
}
