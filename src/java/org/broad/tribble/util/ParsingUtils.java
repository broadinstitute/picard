/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.util;

import net.sf.picard.util.Log;

import java.awt.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Constructor;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class ParsingUtils {

    private static final Log log = Log.getInstance(ParsingUtils.class);

    public static Map<Object, Color> colorCache = new WeakHashMap<Object, Color>(100);

    // HTML 4.1 color table,  + orange and magenta
    static Map<String, String> colorSymbols = new HashMap();
    public static Class urlHelperClass = HTTPHelper.class;

    static {
        colorSymbols.put("white", "FFFFFF");
        colorSymbols.put("silver", "C0C0C0");
        colorSymbols.put("gray", "808080");
        colorSymbols.put("black", "000000");
        colorSymbols.put("red", "FF0000");
        colorSymbols.put("maroon", "800000");
        colorSymbols.put("yellow", "FFFF00");
        colorSymbols.put("olive", "808000");
        colorSymbols.put("lime", "00FF00");
        colorSymbols.put("green", "008000");
        colorSymbols.put("aqua", "00FFFF");
        colorSymbols.put("teal", "008080");
        colorSymbols.put("blue", "0000FF");
        colorSymbols.put("navy", "000080");
        colorSymbols.put("fuchsia", "FF00FF");
        colorSymbols.put("purple", "800080");
        colorSymbols.put("orange", "FFA500");
        colorSymbols.put("magenta", "FF00FF");
    }


    public static InputStream openInputStream(String path)
            throws IOException {

        InputStream inputStream;

        if (path.startsWith("http:") || path.startsWith("https:") || path.startsWith("ftp:")) {
            inputStream = getURLHelper(new URL(path)).openInputStream();
        } else {
            File file = new File(path);
            inputStream = new FileInputStream(file);
        }

        return inputStream;
    }

    //public static String join(String separator, Collection<String> strings) {
    //    return join( separator, strings.toArray(new String[0]) );
    //}

    public static <T> String join(String separator, Collection<T> objects) {
        if (objects.isEmpty()) {
            return "";
        }
        Iterator<T> iter = objects.iterator();
        final StringBuilder ret = new StringBuilder(iter.next().toString());
        while (iter.hasNext()) {
            ret.append(separator);
            ret.append(iter.next().toString());
        }

        return ret.toString();
    }

    /**
     * a small utility function for sorting a list
     *
     * @param list
     * @param <T>
     * @return
     */
    public static <T extends Comparable> List<T> sortList(Collection<T> list) {
        ArrayList<T> ret = new ArrayList<T>();
        ret.addAll(list);
        Collections.sort(ret);
        return ret;
    }

    public static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {
        List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        List<String> pairs = new ArrayList<String>();
        for (T k : t) {
            pairs.add(k + "=" + c.get(k));
        }

        return "{" + ParsingUtils.join(", ", pairs.toArray(new String[pairs.size()])) + "}";
    }

    /**
     * join an array of strings given a seperator
     *
     * @param separator the string to insert between each array element
     * @param strings   the array of strings
     * @return a string, which is the joining of all array values with the separator
     */
    public static String join(String separator, String[] strings) {
        return join(separator, strings, 0, strings.length);
    }

    /**
     * join a set of strings, using the separator provided, from index start to index stop
     *
     * @param separator the separator to use
     * @param strings   the list of strings
     * @param start     the start position (index in the list)0
     * @param end       the end position (index in the list)
     * @return a joined string, or "" if end - start == 0
     */
    public static String join(String separator, String[] strings, int start, int end) {
        if ((end - start) == 0) {
            return "";
        }
        StringBuilder ret = new StringBuilder(strings[start]);
        for (int i = start + 1; i < end; ++i) {
            ret.append(separator);
            ret.append(strings[i]);
        }
        return ret.toString();
    }


    /**
     * Split the string into tokesn separated by the given delimiter.  Profiling has
     * revealed that the standard string.split() method typically takes > 1/2
     * the total time when used for parsing ascii files.
     *
     * @param aString the string to split
     * @param tokens  an array to hold the parsed tokens
     * @param delim   character that delimits tokens
     * @return the number of tokens parsed
     */
    public static int split(String aString, String[] tokens, char delim) {
        return split(aString, tokens, delim, false);
    }

    /**
     * Split the string into tokens separated by the given delimiter.  Profiling has
     * revealed that the standard string.split() method typically takes > 1/2
     * the total time when used for parsing ascii files.
     *
     * @param aString                the string to split
     * @param tokens                 an array to hold the parsed tokens
     * @param delim                  character that delimits tokens
     * @param condenseTrailingTokens if true and there are more tokens than will fit in the tokens array,
     *                               condense all trailing tokens into the last token
     * @return the number of tokens parsed
     */
    public static int split(String aString, String[] tokens, char delim, boolean condenseTrailingTokens) {

        int maxTokens = tokens.length;
        int nTokens = 0;
        int start = 0;
        int end = aString.indexOf(delim);

        if (end == 0) {
            if (aString.length() > 1) {
                start = 1;
                end = aString.indexOf(delim, start);
            } else {
                return 0;
            }
        }

        if (end < 0) {
            tokens[nTokens++] = aString;
            return nTokens;
        }

        while ((end > 0) && (nTokens < maxTokens)) {
            tokens[nTokens++] = aString.substring(start, end);
            start = end + 1;
            end = aString.indexOf(delim, start);
        }

        // condense if appropriate
        if (condenseTrailingTokens && nTokens == maxTokens) {
            tokens[nTokens - 1] = tokens[nTokens - 1] + delim + aString.substring(start);
        }
        // Add the trailing string
        else if (nTokens < maxTokens) {
            String trailingString = aString.substring(start);
            tokens[nTokens++] = trailingString;
        }

        return nTokens;
    }


    // trim a string for the given character (i.e. not just whitespace)

    public static String trim(String str, char ch) {
        char[] array = str.toCharArray();
        int start = 0;
        while (start < array.length && array[start] == ch)
            start++;

        int end = array.length - 1;
        while (end > start && array[end] == ch)
            end--;

        return str.substring(start, end + 1);
    }


    /**
     * Split the string into tokens separated by tab or space(s).  This method
     * was added so support wig and bed files, which apparently accept space delimiters.
     * <p/>
     * Note:  TODO REGEX expressions are not used for speed.  This should be re-evaluated with JDK 1.5 or later
     *
     * @param aString the string to split
     * @param tokens  an array to hold the parsed tokens
     * @return the number of tokens parsed
     */
    public static int splitWhitespace(String aString, String[] tokens) {

        int maxTokens = tokens.length;
        int nTokens = 0;
        int start = 0;
        int tabEnd = aString.indexOf('\t');
        int spaceEnd = aString.indexOf(' ');
        int end = tabEnd < 0 ? spaceEnd : spaceEnd < 0 ? tabEnd : Math.min(spaceEnd, tabEnd);
        while ((end > 0) && (nTokens < maxTokens)) {
            //tokens[nTokens++] = new String(aString.toCharArray(), start, end-start); //  aString.substring(start, end);
            tokens[nTokens++] = aString.substring(start, end);

            start = end + 1;
            // Gobble up any whitespace before next token -- don't gobble tabs, consecutive tabs => empty cell
            while (start < aString.length() && aString.charAt(start) == ' ') {
                start++;
            }

            tabEnd = aString.indexOf('\t', start);
            spaceEnd = aString.indexOf(' ', start);
            end = tabEnd < 0 ? spaceEnd : spaceEnd < 0 ? tabEnd : Math.min(spaceEnd, tabEnd);

        }

        // Add the trailing string
        if (nTokens < maxTokens) {
            String trailingString = aString.substring(start);
            tokens[nTokens++] = trailingString;
        }
        return nTokens;
    }

    public static <T extends Comparable<? super T>> boolean isSorted(Iterable<T> iterable) {
        Iterator<T> iter = iterable.iterator();
        if (!iter.hasNext())
            return true;

        T t = iter.next();
        while (iter.hasNext()) {
            T t2 = iter.next();
            if (t.compareTo(t2) > 0)
                return false;

            t = t2;
        }

        return true;
    }

    /**
     * Convert an rgb string, hex, or symbol to a color.
     *
     * @param string
     * @return
     */
    public static Color parseColor(String string) {
        try {
            Color c = colorCache.get(string);
            if (c == null) {
                if (string.contains(",")) {
                    String[] rgb = string.split(",");
                    int red = Integer.parseInt(rgb[0]);
                    int green = Integer.parseInt(rgb[1]);
                    int blue = Integer.parseInt(rgb[2]);
                    c = new Color(red, green, blue);
                } else if (string.startsWith("#")) {
                    c = hexToColor(string.substring(1));
                } else {
                    String hexString = colorSymbols.get(string.toLowerCase());
                    if (hexString != null) {
                        c = hexToColor(hexString);
                    }
                }

                if (c == null) {
                    c = Color.black;
                }
                colorCache.put(string, c);
            }
            return c;

        } catch (NumberFormatException numberFormatException) {
            //TODO Throw this exception?
            return Color.black;
        }
    }


    private static Color hexToColor(String string) {
        if (string.length() == 6) {
            int red = Integer.parseInt(string.substring(0, 2), 16);
            int green = Integer.parseInt(string.substring(2, 4), 16);
            int blue = Integer.parseInt(string.substring(4, 6), 16);
            return new Color(red, green, blue);
        } else {
            return null;
        }

    }

    public static boolean resourceExists(String resource) throws IOException{

        boolean remoteFile = resource.startsWith("http://") || resource.startsWith("https://") || resource.startsWith("ftp://");
        if (remoteFile) {
            URL url = null;
            try {
                url = new URL(resource);
            } catch (MalformedURLException e) {
                // Malformed URLs by definition don't exist
                return false;
            }
            //TODO Make FTP helper, use it where necessary
            URLHelper helper = getURLHelper(url);
            return helper.exists();
        } else {
            return (new File(resource)).exists();
        }
    }

    public static URLHelper getURLHelper(URL url) {
        try {
            Constructor constr = urlHelperClass.getConstructor(URL.class);
            URLHelper helper = (URLHelper) constr.newInstance(url);
            return helper;
        } catch (Exception e) {
            log.error("Error instantiating url helper for class: " + urlHelperClass, e);
            return new HTTPHelper(url);
        }
    }

    public static void registerHelperClass(Class helperClass) {
        if (!URLHelper.class.isAssignableFrom(helperClass)) {
            throw new IllegalArgumentException("helperClass must implement URLHelper");
        }
        urlHelperClass = helperClass;

    }
}
