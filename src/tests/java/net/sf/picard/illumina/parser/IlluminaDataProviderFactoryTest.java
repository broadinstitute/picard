package net.sf.picard.illumina.parser;


import java.io.File;
import java.util.List;

public class IlluminaDataProviderFactoryTest {

    class TestFactory extends IlluminaDataProviderFactory{
        public TestFactory(final File basecallDirectory, final int lane, final ReadStructure readStructure, final IlluminaDataType... dataTypes) {
            super(basecallDirectory, lane, readStructure, dataTypes);
        }

    }
}
