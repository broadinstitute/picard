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

import net.sf.picard.PicardException;

import java.io.*;
import java.util.concurrent.*;

/**
 * Utility class that will execute sub processes via Runtime.getRuntime().exec(...) and read
 * off the output from stderr and stdout of the sub process. This implementation uses a different
 * thread to read each stream: the current thread for stdout and another, internal thread for 
 * stderr. This utility is able to handle concurrent executions, spawning as many threads as
 * are required to handle the concurrent load.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class ProcessExecutor {
    private static final Log log = Log.getInstance(ProcessExecutor.class);
    private static final ExecutorService executorService = Executors.newCachedThreadPool(new ThreadFactory() {
        @Override
        public Thread newThread(final Runnable r) {
            return new Thread(r, "ProcessExecutor Thread");
        }
    });
    
    /**
     * Executes the command via Runtime.getRuntime().exec() then writes stderr to log.error
     * and stdout to log.info and blocks until the command is complete.
     * 
     * @see Runtime#exec(String)
     * 
     * @param command command string
     * @return return code of command
     */
    public static int execute(final String command) {
        try {
            final Process process = Runtime.getRuntime().exec(command);
            return readStreamsAndWaitFor(process);
        } catch (Throwable t) {
            throw new PicardException("Unexpected exception executing [" + net.sf.samtools.util.StringUtil.join(" ", command) + "]", t);
        }
    }

    /**
     * Executes the command via Runtime.getRuntime().exec() then writes stderr to log.error
     * and stdout to log.info and blocks until the command is complete.
     * 
     * @see Runtime#exec(String[])
     * 
     * @param commandParts command string
     * @return return code of command
     */
    public static int execute(final String[] commandParts) {
        return execute(commandParts, null);
    }

    /**
     * Executes the command via Runtime.getRuntime().exec(), writes <code>outputStreamString</code>
     * to the process output stream if it is not null, then writes stderr to log.error
     * and stdout to log.info and blocks until the command is complete.
     *
     * @see Runtime#exec(String[])
     *
     * @param commandParts command string
     * @return return code of command
     */
    public static int execute(final String[] commandParts, String outputStreamString) {
        try {
            final Process process = Runtime.getRuntime().exec(commandParts);
            if (outputStreamString != null) {
                BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
                writer.write(outputStreamString);
                writer.newLine();
                writer.close();
            }
            return readStreamsAndWaitFor(process);
        } catch (Throwable t) {
            throw new PicardException("Unexpected exception executing [" + net.sf.samtools.util.StringUtil.join(" ", commandParts) + "]", t);
        }
    }

    public static String executeAndReturnResult(final String command) {
        try {
            final Process process = Runtime.getRuntime().exec(command);
            final StringBuilderProcessOutputReader err = new StringBuilderProcessOutputReader(process.getErrorStream());
            final Future<?> stderrReader = executorService.submit(err);
            final StringBuilderProcessOutputReader stdout = new StringBuilderProcessOutputReader(process.getInputStream());
            stdout.run();
            // wait for stderr reader to be done
            stderrReader.get();
            final int result = process.waitFor();
            return result == 0 ? stdout.getOutput() : err.getOutput();
        } catch (Throwable t) {
            throw new PicardException("Unexpected exception executing [" + command + "]", t);
        }

    }

    private static int readStreamsAndWaitFor(final Process process)
            throws InterruptedException, ExecutionException {
        final Future<?> stderrReader = executorService.submit(new LogErrorProcessOutputReader(process.getErrorStream()));
        new LogInfoProcessOutputReader(process.getInputStream()).run();
        // wait for stderr reader to be done
        stderrReader.get();
        return process.waitFor();
    }


    /**
     * Runnable that reads off the given stream and logs it somewhere.
     */
    private static abstract class ProcessOutputReader implements Runnable {
        private final BufferedReader reader;
        public ProcessOutputReader(final InputStream stream) {
            reader = new BufferedReader(new InputStreamReader(stream));
        }

        @Override
        public void run() {
            try {
                String line;
                while ((line = reader.readLine()) != null) {
                    write(line);
                }
            } catch (IOException e) {
                throw new PicardException("Unexpected exception reading from process stream", e);
            }
        }
        
        protected abstract void write(String message);
    }


    private static class LogErrorProcessOutputReader extends ProcessOutputReader {
        public LogErrorProcessOutputReader(final InputStream stream) { super(stream); }
        @Override protected void write(final String message) { log.error(message); }
    }

    private static class LogInfoProcessOutputReader extends ProcessOutputReader {
        public LogInfoProcessOutputReader(final InputStream stream) { super(stream); }
        @Override protected void write(final String message) { log.info(message); }
    }

    private static class StringBuilderProcessOutputReader extends ProcessOutputReader {
        private final StringBuilder sb = new StringBuilder();
        public StringBuilderProcessOutputReader(final InputStream stream) { super(stream); }
        @Override protected void write(final String message) { sb.append(message).append("\n"); }
        public String getOutput() { return sb.toString(); }
    }

}
