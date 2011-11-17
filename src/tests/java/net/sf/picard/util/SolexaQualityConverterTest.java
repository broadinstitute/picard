package net.sf.picard.util;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.util.Arrays;

public class SolexaQualityConverterTest {
    //declared as a staic variable because we reuse it in IlluminaUtilTest
    public static Object[][] SOLEXA_QUALS_TO_PHRED_SCORE = new Object[][] {
                new Object[]{new byte[]{}, new byte[]{}},
                new Object[]{iToB(new int[]{120}), iToB(new int[]{56})},
                new Object[]{iToB(new int[]{64, 65, 120, 121}), iToB(new int[]{3, 4, 56, 57})},
                new Object[]{iToB(new int[]{0, 1, 63}), iToB(new int[]{0, 0, 0})}
    };

    public static Object[][] SOLEXA_QUALS_TO_PHRED_SCORE_1_3 = new Object[][] {
                new Object[]{new byte[]{}, new byte[]{}},
                new Object[]{iToB(new int[]{120}), iToB(new int[]{56})},
                new Object[]{iToB(new int[]{64, 65, 120, 121, 156, 157}), iToB(new int[]{0, 1, 56, 57, 92, 93})},
                new Object[]{iToB(new int[]{0, 1, 63}), iToB(new int[]{-64, -63, -1})}
    };

    public static Object[][] INVALID_SOLEXA_QUALS = new Object[][] {
                new Object[]{iToB(new int[]{-1}), iToB(new int[]{-65})},
                new Object[]{iToB(new int[]{-1, -2}), iToB(new int[]{-65, -66})}
    };

    private static final byte [] iToB(int [] intVals) {
        byte [] byteVals = new byte[intVals.length];
        for(int i = 0; i < byteVals.length; i++) {
            byteVals[i] = (byte) intVals[i];
        }
        return byteVals;
    }

    @DataProvider(name="solexaQualsToPhredScore")
    public Object[][] solexaQualsToPhredScore() {return SOLEXA_QUALS_TO_PHRED_SCORE;}

    @DataProvider(name="solexaQualsToPhredScore_1_3")
    public Object[][] solexaQualsToPhredScore_1_3() {return SOLEXA_QUALS_TO_PHRED_SCORE_1_3;}

    @DataProvider(name="invalidSolexaQuals")
    public Object[][] invalidSolexaQuals() {return INVALID_SOLEXA_QUALS;}

    @Test(dataProvider="solexaQualsToPhredScore")
    public void solexaQualsToPhredScoreTestArray(byte [] solexaQuals, byte [] phredScores) {
       byte [] qualsToPhred = Arrays.copyOf(solexaQuals, solexaQuals.length);

       SolexaQualityConverter.getSingleton().convertSolexaQualityCharsToPhredBinary(qualsToPhred);
       Assert.assertEquals(phredScores, qualsToPhred);
    }

    @Test(dataProvider="solexaQualsToPhredScore_1_3")
    public void solexaQualsToPhredScoreTestArray_1_3(byte [] solexaQuals, byte [] phredScores) {
       byte [] qualsToPhred = Arrays.copyOf(solexaQuals, solexaQuals.length);

       SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(qualsToPhred);
       Assert.assertEquals(phredScores, qualsToPhred);
    }

    @Test(dataProvider="invalidSolexaQuals")
    public void invalidSolexaQualsTestArray(byte [] solexaQuals,  byte [] phredScores) {
       byte [] qualsToPhred = Arrays.copyOf(solexaQuals, solexaQuals.length);

       SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(qualsToPhred);
       Assert.assertEquals(phredScores, qualsToPhred);
    }

    @Test(dataProvider="solexaQualsToPhredScore")
    public void solexaCharToPhredBinaryTest(final byte [] solexaQuals, byte [] phredScores) {
       final SolexaQualityConverter solexaConverter = SolexaQualityConverter.getSingleton();
       for(int i = 0; i < solexaQuals.length; i++) {
           Assert.assertEquals(phredScores[i], solexaConverter.solexaCharToPhredBinary(solexaQuals[i]));
       }
    }

    @Test(dataProvider="invalidSolexaQuals", expectedExceptions = IndexOutOfBoundsException.class)
    public void invalidSolexaCharToPhredBinaryTest(final byte [] solexaQuals,  byte [] phredScores) {
       final SolexaQualityConverter solexaConverter = SolexaQualityConverter.getSingleton();
       for(int i = 0; i < solexaQuals.length; i++) {
           solexaConverter.solexaCharToPhredBinary(solexaQuals[i]);
       }
    }
}
