package picard.util;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;

public class IntervalListToolsTest {
    @Test
    public void testSecondInputValidation(){
        IntervalListTools intervalListTools = new IntervalListTools();
        String[] errors = intervalListTools.customCommandLineValidation();
        Assert.assertNull(errors);

        intervalListTools.SECOND_INPUT = new ArrayList<>();
        errors = intervalListTools.customCommandLineValidation();
        Assert.assertNull(errors);

        intervalListTools.SECOND_INPUT.add(new File("fakefile"));
        errors = intervalListTools.customCommandLineValidation();
        Assert.assertTrue(errors.length == 1);
    }
}
