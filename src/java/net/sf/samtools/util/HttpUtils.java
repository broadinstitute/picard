package net.sf.samtools.util;

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
        URLConnection conn = null;
        try {
            // Create a URLConnection object for a URL
            conn = url.openConnection();
            conn.setReadTimeout(3000);
            return conn.getHeaderField("ETag");
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (conn != null && conn instanceof HttpURLConnection) {
                ((HttpURLConnection) conn).disconnect();
            }
        }
    }

    public static String getHeaderField(final URL url, final String name) {
        URLConnection conn = null;
        try {
            // Create a URLConnection object for a URL
            conn = url.openConnection();
            conn.setReadTimeout(3000);
            return conn.getHeaderField(name);

        } catch (Exception e) {
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
            conn = url.openConnection();
            conn.setReadTimeout(3000);

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
        URLConnection conn = null;
        try {
            // Create a URLConnection object for a URL
            conn = url.openConnection();
            conn.setReadTimeout(3000);
            return conn.getHeaderField("ETag") != null;
        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
        finally {
            if (conn != null && conn instanceof HttpURLConnection) {
                ((HttpURLConnection) conn).disconnect();
            }
        }
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
