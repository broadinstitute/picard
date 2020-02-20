package picard.util;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

public class MiscTest {

    @Test
    public void canWriteToDevNull() throws IOException {
        File f = new File("/dev/null");
        Assert.assertTrue(f.canRead());

        final OutputStream stream = new FileOutputStream(f);
        final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stream));

        writer.write("Just a test");
        writer.close();
    }
}
