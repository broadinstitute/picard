package picard.arrays.illumina;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * A series of tests for the static methods in InfiniumDataFile
 */
public class InfiniumDataFileTest {

    @DataProvider(name = "byteArrayToIntDataProvider")
    public Object[][] byteArrayToIntDataProvider() {
        return new Object[][]{
                { new byte[]{45, 112}, 28717 },
                { new byte[]{-117, 1}, 395 },
                { new byte[]{111, 3}, 879 },
                { new byte[]{-1, 127}, Short.MAX_VALUE },
                { new byte[]{0, -128}, Short.MIN_VALUE },
                { new byte[]{0, 0}, 0 }
        };
    }

    @Test(dataProvider = "byteArrayToIntDataProvider")
    public void testByteArrayToInt(byte[] bytes, int expectedValue) {
        short value = (short) InfiniumDataFile.byteArrayToInt(bytes);
        Assert.assertEquals(value, expectedValue);
    }

    @Test(dataProvider = "byteArrayToIntDataProvider")
    public void testShortToByteArray(byte[] expectedBytes, int value) {
        byte[] themBytes = InfiniumDataFile.shortToByteArray((short) value);
        Assert.assertEquals(themBytes, expectedBytes);
    }

    @DataProvider(name = "byteArrayToFloatDataProvider")
    public Object[][] byteArrayToFloatDataProvider() {
        return new Object[][]{
                { new byte[]{0, 0, -72, 65}, 23.0f },
                { new byte[]{-48, 15, 73, 64}, 3.14159f },
                { new byte[]{-1, -1, 127, 127}, Float.MAX_VALUE },
                { new byte[]{1, 0, 0, 0}, Float.MIN_VALUE },
                { new byte[]{0, 0, 0, 0}, 0.0f }
        };
    }

    @Test(dataProvider = "byteArrayToFloatDataProvider")
    public void testByteArrayToFloat(byte[] bytes, float expectedValue) {
        float value = InfiniumDataFile.byteArrayToFloat(bytes);
        Assert.assertEquals(value, expectedValue);
    }

    @Test(dataProvider = "byteArrayToFloatDataProvider")
    public void testFloatToByteArray(byte[] expectedBytes, float value) {
        byte[] themBytes = InfiniumDataFile.floatToByteArray(value);
        Assert.assertEquals(themBytes, expectedBytes);
    }

    @Test
    public void testByteArrayToCharArray() {
        byte[] themBytes = new byte[]{77, 69, 71, 95};
        char[] foo = InfiniumDataFile.byteArrayToCharArray(themBytes);
        Assert.assertEquals(String.valueOf(foo), "MEG_");
    }
}
