package banchmark;

import org.openjdk.jmh.annotations.*;
import picard.cmdline.PicardCommandLine;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgsAppend = {"-Xms6G", "-Xmx6G"})
@Warmup(iterations = 0)
@Measurement(iterations = 1)
public class SamBenchMark {


    private static final String INPUT = "/Users/artem.balabaev/quantori/test_data/SortSam/input.bam";
    private static final String OUTPUT = "/Users/artem.balabaev/quantori/test_data/SortSam/sorted.bam";
    private static final String module = "SortSam";

    @State(Scope.Benchmark)
    public static class ExecutionPlan {

        @Param({ "50000", "100000", "400000", "500000", "1000000", "1500000", "2000000", "2500000"})
        public String bufferSize;

        @Setup(Level.Invocation)
        public void setUp() {
        }
    }

    @Benchmark
    public void sortingBufferParameter(ExecutionPlan plan) {

        final List<String> listOfArgs = new ArrayList<String>() {{
            add(module);
            add(prepareParam("INPUT", INPUT));
            add(prepareParam("OUTPUT", OUTPUT));
            add(prepareParam("SORT_ORDER", "coordinate"));
            add(prepareParam("VALIDATION_STRINGENCY", "SILENT"));
            add(prepareParam("MAX_RECORDS_IN_RAM", plan.bufferSize));
            add(prepareParam("USE_JDK_DEFLATER", "true"));
            add(prepareParam("USE_JDK_INFLATER", "true"));
        }};


        PicardCommandLine.startProcess(listOfArgs.toArray(new String[0]));
    }

    private static String prepareParam(String key, String value) {
        return key + "=" + value;
    }
}
