// DateParser.java
// $Id: DateParser.java,v 1.3 2001/01/04 13:26:19 bmahe Exp $
// (c) COPYRIGHT MIT, INRIA and Keio, 2000.

/*
W3C IPR SOFTWARE NOTICE

Copyright 1995-1998 World Wide Web Consortium, (Massachusetts Institute of
Technology, Institut National de Recherche en Informatique et en
Automatique, Keio University). All Rights Reserved.
http://www.w3.org/Consortium/Legal/

This W3C work (including software, documents, or other related items) is
being provided by the copyright holders under the following license. By
obtaining, using and/or copying this work, you (the licensee) agree that you
have read, understood, and will comply with the following terms and
conditions:

Permission to use, copy, and modify this software and its documentation,
with or without modification,  for any purpose and without fee or royalty is
hereby granted, provided that you include the following on ALL copies of the
software and documentation or portions thereof, including modifications,
that you make:

  1. The full text of this NOTICE in a location viewable to users of the
     redistributed or derivative work.
  2. Any pre-existing intellectual property disclaimers, notices, or terms
     and conditions. If none exist, a short notice of the following form
     (hypertext is preferred, text is permitted) should be used within the
     body of any redistributed or derivative code: "Copyright World Wide
     Web Consortium, (Massachusetts Institute of Technology, Institut
     National de Recherche en Informatique et en Automatique, Keio
     University). All Rights Reserved. http://www.w3.org/Consortium/Legal/"
  3. Notice of any changes or modifications to the W3C files, including the
     date changes were made. (We recommend you provide URIs to the location
     from which the code is derived).

In addition, creators of derivitive works must include the full text of this
NOTICE in a location viewable to users of the derivitive work.

THIS SOFTWARE AND DOCUMENTATION IS PROVIDED "AS IS," AND COPYRIGHT HOLDERS
MAKE NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT
LIMITED TO, WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
PURPOSE OR THAT THE USE OF THE SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.

COPYRIGHT HOLDERS WILL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THE SOFTWARE OR
DOCUMENTATION.

The name and trademarks of copyright holders may NOT be used in advertising
or publicity pertaining to the software without specific, written prior
permission. Title to copyright in this software and any associated
documentation will at all times remain with copyright holders.

____________________________________

This formulation of W3C's notice and license became active on August 14
1998. See the older formulation for the policy prior to this date. Please
see our Copyright FAQ for common questions about using materials from our
site, including specific terms and conditions for packages like libwww,
Amaya, and Jigsaw. Other questions about this notice can be directed to
site-policy@w3.org .




webmaster
(last updated 14-Aug-1998)

 */

package net.sf.samtools.util;

import net.sf.samtools.SAMException;

import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;
import java.util.TimeZone;

/**
 * NOTE: This code has been taken from w3.org, and modified slightly to handle timezones of the form [-+]DDDD,
 * and also to fix a bug in the application of time zone to the parsed date.
 *
 * Date parser for ISO 8601 format
 * http://www.w3.org/TR/1998/NOTE-datetime-19980827
 * @version $Revision: 1.3 $
 * @author  bmahe@w3.org
 */
public class DateParser {

    private static boolean check(StringTokenizer st, String token)
            throws InvalidDateException
    {
        if (!st.hasMoreElements()) return false;
        if (st.nextToken().equals(token)) {
            return true;
        } else {
            throw new InvalidDateException("Missing ["+token+"]");
        }
    }

    private static Calendar getCalendar(String isodate)
            throws InvalidDateException
    {
        // YYYY-MM-DDThh:mm:ss.sTZD
        StringTokenizer st = new StringTokenizer(isodate, "-T:.+Z", true);

        Calendar calendar = new GregorianCalendar(TimeZone.getTimeZone("UTC"));
        calendar.clear();
        try {
            // Year
            if (st.hasMoreTokens()) {
                int year = Integer.parseInt(st.nextToken());
                calendar.set(Calendar.YEAR, year);
            } else {
                return calendar;
            }
            // Month
            if (check(st, "-") && (st.hasMoreTokens())) {
                int month = Integer.parseInt(st.nextToken()) -1;
                calendar.set(Calendar.MONTH, month);
            } else {
                return calendar;
            }
            // Day
            if (check(st, "-") && (st.hasMoreTokens())) {
                int day = Integer.parseInt(st.nextToken());
                calendar.set(Calendar.DAY_OF_MONTH, day);
            } else {
                return calendar;
            }
            // Hour
            if (check(st, "T") && (st.hasMoreTokens())) {
                int hour = Integer.parseInt(st.nextToken());
                calendar.set(Calendar.HOUR_OF_DAY, hour);
            } else {
                calendar.set(Calendar.HOUR_OF_DAY, 0);
                calendar.set(Calendar.MINUTE, 0);
                calendar.set(Calendar.SECOND, 0);
                calendar.set(Calendar.MILLISECOND, 0);
                return calendar;
            }
            // Minutes
            if (check(st, ":") && (st.hasMoreTokens())) {
                int minutes = Integer.parseInt(st.nextToken());
                calendar.set(Calendar.MINUTE, minutes);
            } else {
                calendar.set(Calendar.MINUTE, 0);
                calendar.set(Calendar.SECOND, 0);
                calendar.set(Calendar.MILLISECOND, 0);
                return calendar;
            }

            //
            // Not mandatory now
            //

            // Secondes
            if (! st.hasMoreTokens()) {
                return calendar;
            }
            String tok = st.nextToken();
            if (tok.equals(":")) { // secondes
                if (st.hasMoreTokens()) {
                    int secondes = Integer.parseInt(st.nextToken());
                    calendar.set(Calendar.SECOND, secondes);
                    if (! st.hasMoreTokens()) {
                        return calendar;
                    }
                    // frac sec
                    tok = st.nextToken();
                    if (tok.equals(".")) {
                        // bug fixed, thx to Martin Bottcher
                        String nt = st.nextToken();
                        while(nt.length() < 3) {
                            nt += "0";
                        }
                        nt = nt.substring( 0, 3 ); //Cut trailing chars..
                        int millisec = Integer.parseInt(nt);
                        //int millisec = Integer.parseInt(st.nextToken()) * 10;
                        calendar.set(Calendar.MILLISECOND, millisec);
                        if (! st.hasMoreTokens()) {
                            return calendar;
                        }
                        tok = st.nextToken();
                    } else {
                        calendar.set(Calendar.MILLISECOND, 0);
                    }
                } else {
                    throw new InvalidDateException("No secondes specified");
                }
            } else {
                calendar.set(Calendar.SECOND, 0);
                calendar.set(Calendar.MILLISECOND, 0);
            }
            // Timezone
            if (! tok.equals("Z")) { // UTC
                if (! (tok.equals("+") || tok.equals("-"))) {
                    throw new InvalidDateException("only Z, + or - allowed");
                }
                boolean plus = tok.equals("+");
                if (! st.hasMoreTokens()) {
                    throw new InvalidDateException("Missing hour field");
                }
                int tzhour = Integer.parseInt(st.nextToken());
                int tzmin  = 0;
                if (check(st, ":") && (st.hasMoreTokens())) {
                    tzmin = Integer.parseInt(st.nextToken());
                } else {
                    // Modified by AW -- minute field is not required, and minutes may be represented
                    // without colon between hours and minutes
                    // throw new InvalidDateException("Missing minute field");
                    if (tzhour >= 100) {
                        tzmin = tzhour % 100;
                        tzhour /= 100;
                    }
                }
                if (!plus) { // Modified by AW -- !plus instead of plus
                    calendar.add(Calendar.HOUR, tzhour);
                    calendar.add(Calendar.MINUTE, tzmin);
                } else {
                    calendar.add(Calendar.HOUR, -tzhour);
                    calendar.add(Calendar.MINUTE, -tzmin);
                }
            }
        } catch (NumberFormatException ex) {
            throw new InvalidDateException("["+ex.getMessage()+
                    "] is not an integer");
        }
        return calendar;
    }

    /**
     * Parse the given string in ISO 8601 format and build a Date object.
     * @param isodate the date in ISO 8601 format
     * @return a Date instance
     * @exception InvalidDateException if the date is not valid
     */
    public static Date parse(String isodate)
            throws InvalidDateException
    {
        Calendar calendar = getCalendar(isodate);
        return calendar.getTime();
    }

    private static String twoDigit(int i) {
        if (i >=0 && i < 10) {
            return "0"+String.valueOf(i);
        }
        return String.valueOf(i);
    }

    /**
     * Generate a ISO 8601 date
     * @param date a Date instance
     * @return a string representing the date in the ISO 8601 format
     */
    public static String getIsoDate(Date date) {
        Calendar calendar = new GregorianCalendar(TimeZone.getTimeZone("UTC"));
        calendar.setTime(date);
        StringBuffer buffer = new StringBuffer();
        buffer.append(calendar.get(Calendar.YEAR));
        buffer.append("-");
        buffer.append(twoDigit(calendar.get(Calendar.MONTH) + 1));
        buffer.append("-");
        buffer.append(twoDigit(calendar.get(Calendar.DAY_OF_MONTH)));
        buffer.append("T");
        buffer.append(twoDigit(calendar.get(Calendar.HOUR_OF_DAY)));
        buffer.append(":");
        buffer.append(twoDigit(calendar.get(Calendar.MINUTE)));
        buffer.append(":");
        buffer.append(twoDigit(calendar.get(Calendar.SECOND)));
        buffer.append(".");
        buffer.append(twoDigit(calendar.get(Calendar.MILLISECOND) / 10));
        buffer.append("Z");
        return buffer.toString();
    }

    public static void test(String isodate) {
        System.out.println("----------------------------------");
        try {
            Date date = parse(isodate);
            System.out.println(">> "+isodate);
            System.out.println(">> "+date.toString()+" ["+date.getTime()+"]");
            System.out.println(">> "+getIsoDate(date));
        } catch (InvalidDateException ex) {
            System.err.println(isodate+" is invalid");
            System.err.println(ex.getMessage());
        }
        System.out.println("----------------------------------");
    }

    public static void test(Date date) {
        String isodate = null;
        System.out.println("----------------------------------");
        try {
            System.out.println(">> "+date.toString()+" ["+date.getTime()+"]");
            isodate = getIsoDate(date);
            System.out.println(">> "+isodate);
            date = parse(isodate);
            System.out.println(">> "+date.toString()+" ["+date.getTime()+"]");
        } catch (InvalidDateException ex) {
            System.err.println(isodate+" is invalid");
            System.err.println(ex.getMessage());
        }
        System.out.println("----------------------------------");
    }

    public static void main(String args[]) {
        test("1997-07-16T19:20:30.45-02:00");
        test("1997-07-16T19:20:30+01:00");
        test("1997-07-16T19:20:30+01:00");
        test("1997-07-16T19:20");
        test("1997-07-16");
        test("1997-07");
        test("1997");
        test(new Date());
    }

    public static class InvalidDateException extends SAMException {
        public InvalidDateException() {
        }

        public InvalidDateException(final String s) {
            super(s);
        }

        public InvalidDateException(final String s, final Throwable throwable) {
            super(s, throwable);
        }

        public InvalidDateException(final Throwable throwable) {
            super(throwable);
        }
    }
}
