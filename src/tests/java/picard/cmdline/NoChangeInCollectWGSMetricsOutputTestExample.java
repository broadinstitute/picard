package picard.cmdline;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.analysis.CollectWgsMetrics;

import java.io.File;

public class NoChangeInCollectWGSMetricsOutputTestExample {
    private final static String picardJarDir = System.getenv().get("PICARD_JAR_DIR"); //for example "/usr/bin/Picard/dist"
    private final static String referenceFasta = System.getenv().get("REFERENCE"); //for example "/data/hg19.fa"

    static {
        if (picardJarDir == null)
            throw new PicardException("PICARD_JAR_DIR environment variable must be set to the latest released Picard jars directory.");
    }

    @Test(dataProvider = "args")
    public void testRunWithDifferentArgs(final String inputBamPath, final String validationStringency) {

        // for the given input parameters, rerun the CLP (and if necessary, the released CLP also)
        // then compare output files.
        new NoChangeInCommandLineProgramOutputTester(false) {

            @Override
            public InputFileArg[] getInputFileArgs() {
                return new InputFileArg[] {
                        new InputFileArg("INPUT", inputBamPath),
                };
            }


            @Override
            public Arg[] getOtherArgs() {
                return new Arg[] {
                        new Arg("R", referenceFasta),
                        new Arg("VALIDATION_STRINGENCY", validationStringency),
                        new Arg("STOP_AFTER", "200000")
                };
            }

            @Override
            public OutputFileArg[] getOutputFileArgs() {
                String[] linesToIgnore = new String[] { "#.*" };
                return new OutputFileArg[] {
                        new OutputFileArg("OUTPUT", "metrics.txt", new TextOutputFileComparator(linesToIgnore))
                };
            }

            @Override
            public File getReleasedCommandLineProgramJarPath() {
                return new File(picardJarDir, "CollectWgsMetrics.jar");
            }

            @Override
            public CommandLineProgram getProgram() {
                return new CollectWgsMetrics();
            }

        }.runTest();
    }

    @DataProvider(name = "args")
    public Object[][] argsDataProvider() {
        return new Object[][] {
                //{ "example.bam", "SILENT" },
                //{ "example.bam", "LENIENT"},

                //{ "example2.bam", "LENIENT"},
                //{ "example2.bam", "SILENT"}
        };
    }

}