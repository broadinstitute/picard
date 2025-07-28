package picard.cmdline.argumentcollections;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class RequesterPaysArgumentCollectionTest {

    public static final String P1 = "project1";
    public static final String P2 = "project2";

    @DataProvider
    public static Object[][] settings() {

        return new Object[][]{
                {null, null, null},
                {"", null, null},
                {null, "", null},
                {"", "", null},
                {P1, null, P1},
                {null, P2, P2},
                {P1, P2, P1},
                {"", P2, P2},
                {P1, "", P1}
        };
    }

    @Test(dataProvider = "settings")
    public void testCorrectValues(String arg, String sys, String expected){
        runWithSystemProperty(
                () -> {
                    final RequesterPaysArgumentCollection rpc = new RequesterPaysArgumentCollection();
                    rpc.requesterPaysProject = arg;
                    Assert.assertEquals(rpc.getProjectForRequesterPays(), expected);
                }, RequesterPaysArgumentCollection.PROPERTY_NAME, sys
        );
    }

    @Test(dataProvider = "settings")
    public void testCorrectDescription(String arg, String sys, String expected){
        runWithSystemProperty(
                () -> {
                    final RequesterPaysArgumentCollection rpc = new RequesterPaysArgumentCollection();
                    rpc.requesterPaysProject = arg;
                    final String description = rpc.getDescription();
                    final String value = rpc.getProjectForRequesterPays();
                    if(expected == null) {Assert.assertEquals(description, "Requester Pays Project not set."); }
                    else if(expected.equals(P1)) { Assert.assertEquals(description, "Requester Pays Project set by argument: " + value); }
                    else if(expected.equals(P2)) { Assert.assertEquals(description, "Requester Pays Project set by system property: " + value); }
                    else { Assert.fail("should have been one of the of the previous"); }
                }, RequesterPaysArgumentCollection.PROPERTY_NAME, sys);
    }
    private static void runWithSystemProperty(Runnable toRun, String name, String value){
        String previousValue = null;
        try {
            if(value != null) {
                previousValue = System.setProperty(name, value);
            }

            toRun.run();

        } finally {
            if(previousValue == null){
                System.clearProperty(name);
            } else {
                System.setProperty(name, previousValue);
            }
        }
    }
}