package picard.util.help;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.ServiceLoader;
import java.util.Set;
import java.util.spi.ToolProvider;

/**
 * Smoke test to run doc gen on a subset of classes to make sure it doesn't regress.
 */
public class DocumentationGenerationIntegrationTest {
    /**
     * Entry point for manually running the picardDoc process on a subset of packages from within picard.
     */
    private static String[] docTestPackages = {
            // note that the fingerprint package can't be in this list, since it requires the google
            // cloud nio library to be present to compile
            "picard.analysis",
            "picard.metrics"
    };

    @Test
    public static void documentationSmokeTest() throws IOException {
        final Path docTestTarget =  Files.createTempDirectory("docgentest");
        final String[] argArray = new String[]{
                //"javadoc",
                "-doclet", PicardHelpDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/picard/helpTemplates",
                "-d", docTestTarget.toFile().getAbsolutePath(), // directory must exist
                "-output-file-extension", "html",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-absolute-version", "1.1-111",
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.add("-cp");
        docArgList.add(System.getProperty("java.class.path"));
        docArgList.addAll(Arrays.asList(docTestPackages));

        // Run javadoc in the current JVM with the custom doclet, and make sure it succeeds (this is a smoke test;
        // we just want to make sure it doesn't blow up).
        final ToolProvider jdProvider = ToolProvider.findFirst("javadoc")
                .orElseThrow(() -> new IllegalStateException("Can't find javadoc tool"));

        final String[] args = docArgList.toArray(new String[] {});
        final int retCode = jdProvider.run(System.out, System.err, args);

        // make sure the task succeeded, and generated at least one index file, plus some other files
        Assert. assertEquals(retCode, 0);

        final File outputDir = docTestTarget.toFile();
        final Set<String> pathNames = new HashSet<>(Arrays.stream(outputDir.list()).toList());
        Assert.assertTrue(pathNames.size() > 1);
        Assert.assertTrue(pathNames.contains("index.html"));
    }
}
