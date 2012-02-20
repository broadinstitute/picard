package net.sf.picard;

import net.sf.samtools.util.StopWatch;
import net.sf.samtools.util.StringUtil;

import java.util.StringTokenizer;
import java.util.regex.Pattern;

/**
 *
 */
public class Test {
    private final String text = "C0A69ACXX111213:6:1101:10000:144257\t83\t5\t128984606\t60\t76M\t=\t128984542\t-140\tAGTGTTAGAACTTCCTCCCCAAAGCATATACTTCAGTGGCAAGCTGTCCTGGATGAAGGTATGACCAACCAGATCA\t@FFFEECC>EFHBJIGIFGIEIJJJIHED<IEHIGIIJIIIGJIGJIIIIIJGGCJIIGIHHHBHGHFFDFFFC@@\tXT:A:U\tNM:i:0\tSM:i:37\tAM:i:37\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:76";

    public static void main(String[] args) {
        new Test().run();
    }

    public void run() {
        final int ITERATIONS = 1000000;
        final String[] fields = new String[10000];
        final StopWatch watch = new StopWatch();

        watch.start();
        for (int i=0; i<ITERATIONS; ++i) {
            if (StringUtil.split(text, fields, '\t') > 100) {
                System.out.println("Mama Mia that's a lot of tokens!!");
            }
        }
        watch.stop();
        System.out.println("StringUtil.split() took " + watch.getElapsedTime());
        watch.reset();
        
        watch.start();
        for (int i=0; i<ITERATIONS; ++i) {
            if (split(text, fields, "\t") > 100) {
                System.out.println("Mama Mia that's a lot of tokens!!");
            }
        }
        watch.stop();
        System.out.println("StringTokenizer took " + watch.getElapsedTime());
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
