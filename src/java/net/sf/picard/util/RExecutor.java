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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import net.sf.picard.util.ProcessExecutor;
import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;

/**
 * Util class for executing R scripts.
 * 
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class RExecutor {
    private static final String R_EXE = "Rscript";
    
    /**
     * Executes the given R script that is stored in a file on the classpath. The script file
     * is read from the classpath and written to a temp file then executed by a call to Rscript.
     * Blocks until the R script is complete.
     * 
     * @param rScriptName the fully qualified name of the classpath resource of the script
     * @param arguments any arguments required by the script
     * @return the return code of the R process
     */
    public static int executeFromClasspath(String rScriptName, String... arguments) {
        File scriptFile = writeScriptFile(rScriptName);
        int returnCode = executeFromFile(scriptFile, arguments);
        net.sf.samtools.util.IOUtil.deleteFiles(scriptFile);
        return returnCode;
    }

    /**
     * Executes the given R script that is stored in a file by a call to Rscript.
     * Blocks until the R script is complete.
     * 
     * @param scriptFile the file object for the script
     * @param arguments any arguments required by the script
     * @return the return code of the R process
     */
    public static int executeFromFile(File scriptFile, String... arguments) {
        String[] command = new String[arguments.length + 2];
        command[0] = R_EXE;
        command[1] = scriptFile.getAbsolutePath();
        System.arraycopy(arguments, 0, command, 2, arguments.length);
        return ProcessExecutor.execute(command);
    }

    /**
     * Writes the classpath resource named by rScriptName to the temp dir.
     */
    private static File writeScriptFile(String rScriptName) {
        InputStream scriptStream = null;
        OutputStream scriptFileStream = null;
        try {
            scriptStream = RExecutor.class.getClassLoader().getResourceAsStream(rScriptName);
            if (scriptStream == null) {
                throw new IllegalArgumentException("Script [" + rScriptName + "] not found in classpath");
            }
            File scriptFile = File.createTempFile("script", ".R");
            scriptFileStream = IoUtil.openFileForWriting(scriptFile);
            IoUtil.copyStream(scriptStream, scriptFileStream);
            return scriptFile;
        } catch (IOException e) {
            throw new PicardException("Unexpected exception creating R script file", e);
        } finally {
            if (scriptStream != null) {
                try {
                    scriptStream.close();
                } catch (IOException e) {
                }
            }
            if (scriptFileStream != null) {
                try {
                    scriptFileStream.close();
                } catch (IOException e) {
                }
            }
        }
    }
}
