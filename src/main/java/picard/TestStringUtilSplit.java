package picard;

import htsjdk.samtools.util.StopWatch;
import htsjdk.samtools.util.StringUtil;

import java.util.StringTokenizer;

/**
 * A class that can be used to check which is faster, StringUtil.split or StringTokenizer.
 * Currently, StringUtil is faster.
 */
public class TestStringUtilSplit {
    private static final String TEXT = "C0A69ACXX111213:6:1101:10000:144257\t83\t5\t128984606\t60\t76M\t=\t128984542\t-140\tAGTGTTAGAACTTCCTCCCCAAAGCATATACTTCAGTGGCAAGCTGTCCTGGATGAAGGTATGACCAACCAGATCA\t@FFFEECC>EFHBJIGIFGIEIJJJIHED<IEHIGIIJIIIGJIGJIIIIIJGGCJIIGIHHHBHGHFFDFFFC@@\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:37\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:76";

    public static void main(String[] args) {
        new TestStringUtilSplit().run();
    }

    public double run() {
        final int ITERATIONS = 1000000;
        final String[] fields = new String[10000];
        final StopWatch watch = new StopWatch();
        final long stringUtilTime;
        watch.start();
        for (int i=0; i<ITERATIONS; ++i) {
            if (StringUtil.split(TEXT, fields, '\t') > 100) {
                System.out.println("Mama Mia that's a lot of tokens!!");
            }
        }
        watch.stop();
        System.out.println("StringUtil.split() took " + watch.getElapsedTime());
        stringUtilTime = watch.getElapsedTime();
        watch.reset();

        final long tokenizerTime;

        watch.start();
        for (int i=0; i<ITERATIONS; ++i) {
            if (split(TEXT, fields, "\t") > 100) {
                System.out.println("Mama Mia that's a lot of tokens!!");
            }
        }
        watch.stop();
        System.out.println("StringTokenizer took " + watch.getElapsedTime());
        tokenizerTime = watch.getElapsedTime();
        return (tokenizerTime + 0D) / stringUtilTime;
    }

    public int split(final String s, final String[] tokens, final String token) {
        final StringTokenizer tokenizer = new StringTokenizer(s, token, false);
        int i=0;
        while (tokenizer.hasMoreTokens()) {
            tokens[i++] = tokenizer.nextToken();
        }

        return i;
    }
}
