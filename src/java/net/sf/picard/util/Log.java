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
package net.sf.picard.util;

import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

/**
 * <p>A <em>wafer thin</em> wrapper around System.out that uses var-args to make it
 * much more efficient to call the logging methods in without having to
 * surround every call site with calls to Log.isXXXEnabled().  All the methods on this
 * class take a variable length list of arguments and, only if logging is enabled for
 * the level and channel being logged to, will those arguments be toString()'d and
 * appended together.</p>
 *
 * @author Tim Fennell
 */
public final class Log {
    /** Enumeration for setting log levels. */
    public static enum LogLevel { ERROR, WARNING, INFO, DEBUG }

    private static LogLevel globalLogLevel = LogLevel.DEBUG;

    private final Class<?> clazz;
    private final String className;
    private final PrintStream out = System.out;

    /**
     * Private constructor
     */
    private Log(final Class<?> clazz) {
        this.clazz = clazz;
        this.className = clazz.getSimpleName();
    }

    /**
     * Get a Log instance to perform logging within the Class specified.  Returns an instance
     * of this class which wraps an instance of the commons logging Log class.
     * @param clazz the Class which is going to be doing the logging
     * @return a Log instance with which to log
     */
    public static Log getInstance(final Class<?> clazz) {
        return new Log(clazz);
    }

    public static void setGlobalLogLevel(final LogLevel logLevel) {
        globalLogLevel = logLevel;
    }

    /** Returns true if the specified log level is enabled otherwise false. */
    public final boolean isEnabled(final LogLevel level) {
        return level.ordinal() <= globalLogLevel.ordinal();
    }

    /**
     * Private method that does the actual printing of messages to a PrintWriter. Outputs the log level,
     * class name and parts followed by the stack trace if a throwable is provided.
     *
     * @param level the Log level being logged at
     * @param throwable a Throwable if one is available otherwise null
     * @param parts the parts of the message to be concatenated
     */
    private void emit(final LogLevel level, final Throwable throwable, final Object... parts) {
        if (isEnabled(level)) {
            this.out.print(level.name());
            this.out.print('\t');
            this.out.print(getTimestamp());
            this.out.print('\t');
            this.out.print(this.className);
            this.out.print('\t');

            for (final Object part : parts) {
                if (part != null && part.getClass().isArray()) {
                    final Class<?> component = part.getClass().getComponentType();
                    if (component.equals(Boolean.TYPE))        this.out.print(Arrays.toString( (boolean[]) part));
                    else if (component.equals(Byte.TYPE))      this.out.print(Arrays.toString( (byte[]) part));
                    else if (component.equals(Character.TYPE)) this.out.print(Arrays.toString( (char[]) part));
                    else if (component.equals(Double.TYPE))    this.out.print(Arrays.toString( (double[]) part));
                    else if (component.equals(Float.TYPE))     this.out.print(Arrays.toString( (float[]) part));
                    else if (component.equals(Integer.TYPE))   this.out.print(Arrays.toString( (int[]) part));
                    else if (component.equals(Long.TYPE))      this.out.print(Arrays.toString( (long[]) part));
                    else if (component.equals(Short.TYPE))     this.out.print(Arrays.toString( (short[]) part));
                    else this.out.print(Arrays.toString( (Object[]) part));
                }
                else {
                    this.out.print(part);
                }
            }

            this.out.println();

            // Print out the exception if there is one
            if (throwable != null) {
                throwable.printStackTrace(this.out);
            }
        }
    }

    /**
     * Creates a date string for insertion into the log.  Given that logs are sometimes held statically
     * and SimpleDateFormat is not thread safe, currently creates an instance each time :/
     */
    protected String getTimestamp() {
        final DateFormat fmt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        return fmt.format(new Date());
    }

    /**
     * Logs a Throwable and optional message parts at level error.
     * @param throwable an instance of Throwable that should be logged with stack trace
     * @param messageParts zero or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void error(final Throwable throwable, final Object... messageParts) {
        emit(LogLevel.ERROR, throwable, messageParts);
    }

    /**
     * Logs a Throwable and optional message parts at level warn.
     * @param throwable an instance of Throwable that should be logged with stack trace
     * @param messageParts zero or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void warn(final Throwable throwable, final Object... messageParts) {
        emit(LogLevel.WARNING, throwable, messageParts);
    }

    /**
     * Logs a Throwable and optional message parts at level info.
     * @param throwable an instance of Throwable that should be logged with stack trace
     * @param messageParts zero or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void info(final Throwable throwable, final Object... messageParts) {
        emit(LogLevel.INFO, throwable, messageParts);
    }

    /**
     * Logs a Throwable and optional message parts at level debug.
     * @param throwable an instance of Throwable that should be logged with stack trace
     * @param messageParts zero or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void debug(final Throwable throwable, final Object... messageParts) {
        emit(LogLevel.DEBUG, throwable, messageParts);
    }

    // Similar methods, but without Throwables, follow

    /**
     * Logs one or more message parts at level error.
     * @param messageParts one or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void error(final Object... messageParts) {
        emit(LogLevel.ERROR, null, messageParts);
    }

    /**
     * Logs one or more message parts at level warn.
     * @param messageParts one or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void warn(final Object... messageParts) {
        emit(LogLevel.WARNING, null, messageParts);
    }

    /**
     * Logs one or more message parts at level info.
     * @param messageParts one or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void info(final Object... messageParts) {
        emit(LogLevel.INFO, null, messageParts);
    }

    /**
     * Logs one or more message parts at level debug.
     * @param messageParts one or more objects which should be combined, by calling toString()
     *        to form the log message.
     */
    public final void debug(final Object... messageParts) {
        emit(LogLevel.DEBUG, null, messageParts);
    }
}
