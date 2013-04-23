package net.sf.samtools.util;

import java.io.IOException;
import java.net.URLConnection;
import java.net.URL;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;

/**
 * User: jrobinso
 * Date: Sep 23, 2009
 */
public class HttpUtils {

    public static String getETag(final URL url) {
        return getHeaderField(url, "ETag");
    }

    private static URLConnection openConnection(final URL url) throws IOException{
        URLConnection conn = url.openConnection();
        conn.setReadTimeout(3000);
        conn.setDefaultUseCaches(false);
        conn.setUseCaches(false);
        return conn;
    }

    public static String getHeaderField(final URL url, final String name) {
        URLConnection conn = null;
        try {
            // Create a URLConnection object for a URL
            conn = openConnection(url);
            return conn.getHeaderField(name);

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (conn != null && conn instanceof HttpURLConnection) {
                ((HttpURLConnection) conn).disconnect();
            }
        }
    }

    public static void printHeaderFields(final URL url) {

        URLConnection conn = null;
        try {
            // Create a URLConnection object for a URL
            conn = openConnection(url);

            for (final String name : conn.getHeaderFields().keySet()) {
                System.out.println(name + "\t" + conn.getHeaderField(name));

            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        finally {
            if (conn != null && conn instanceof HttpURLConnection) {
                ((HttpURLConnection) conn).disconnect();
            }
        }
    }

    public static boolean resourceAvailable(final URL url) {
        return getETag(url) != null;
    }

    public static void main(final String[] args) throws MalformedURLException {
        //printHeaderFields(new URL(
        //        "http://www.broadinstitute.org/igvdata/1KG/DCC_merged/freeze5/NA12891.pilot2.SLX.bam"));
        System.out.println(getETag(new URL(
                 "http://www.broadinstitute.org/igvdata/test/sam/303KY.8.paired1.bam.tdf")));
        System.out.println(resourceAvailable(new URL(
                "http://www.broadinstitute.org/igvdata/test/sam/303KY.8.paired1.bam.tdf")));


    }
}
