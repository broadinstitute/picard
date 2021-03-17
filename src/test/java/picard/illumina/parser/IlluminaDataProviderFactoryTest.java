package picard.illumina.parser;


import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Set;

public class IlluminaDataProviderFactoryTest {

    class TestFactory extends IlluminaDataProviderFactory{
        public TestFactory(final File basecallDirectory, final int lane, final ReadStructure readStructure, final Set<IlluminaDataType> dataTypes) {
            super(basecallDirectory, lane, readStructure, new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY), dataTypes);
        }

    }
}
