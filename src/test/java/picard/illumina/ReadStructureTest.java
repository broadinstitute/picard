package picard.illumina;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;

import java.util.ArrayList;
import java.util.List;

import static htsjdk.samtools.util.CollectionUtil.makeList;
import static picard.illumina.parser.ReadType.*;

public class ReadStructureTest {

    //to make construction of lists more intelligible
    public ReadDescriptor rd(final int length, final ReadType rt) {
        return new ReadDescriptor(length, rt);
    }

    //Many of these readStructures would be non-sensical but check different test classes/combinations
    @DataProvider(name="validReadStructures")
    public Object[][] validReadStructures() {
        return new Object[][] {
            {"2T",                    makeList(rd(2, T)),                 1, 0, 0, 0},
            {"1234B",                 makeList(rd(1234, B)),              0, 1, 0, 0},
            {Integer.MAX_VALUE + "S", makeList(rd(Integer.MAX_VALUE, S)), 0, 0, 1, 0},
            {Integer.MAX_VALUE + "M", makeList(rd(Integer.MAX_VALUE, M)), 0, 0, 0, 1},


            {"76T76T", makeList(rd(76, T), rd(76,  T)), 2, 0, 0, 0},
            {"76T1B",  makeList(rd(76, T), rd(1,   B)), 1, 1, 0, 0},
            {"76B1T",  makeList(rd(76, B), rd(1,   T)), 1, 1, 0, 0},
            {"1S1B",   makeList(rd(1,  S), rd(1,   B)), 0, 1, 1, 0},
            {"1S1B1M", makeList(rd(1,  S), rd(1,   B), rd(1, M)), 0, 1, 1, 1},
            {"1T999S", makeList(rd(1,  T), rd(999, S)), 1, 0, 1, 0},

            {"100T20T100T",  makeList(rd(100, T),  rd(20,  T), rd(100,  T)),    3, 0, 0, 0},
            {"2S50S10S",     makeList(rd(2,   S),  rd(50,  S), rd(10,   S)),    0, 0, 3, 0},
            {"10T1B11T",     makeList(rd(10,  T),  rd(1,   B), rd(11,   T)),    2, 1, 0, 0},
            {"201T13T111B",  makeList(rd(201, T),  rd(13,  T), rd(111,  B)),    2, 1, 0, 0},
            {"15B1T1T",      makeList(rd(15,  B),  rd(1,   T), rd(1,    T)),    2, 1, 0, 0},
            {"99B7T6B",      makeList(rd(99,  B),  rd(7,   T), rd(6,    B)),    1, 2, 0, 0},
            {"631B776S638T", makeList(rd(631, B),  rd(776, S), rd(638,  T)),    1, 1, 1, 0},
            {"631M776S638T", makeList(rd(631, M),  rd(776, S), rd(638,  T)),    1, 0, 1, 1},


            {"3T7B60S2T",     makeList(rd(3,  T),    rd(7,  B),   rd(60,   S),    rd(2,  T)),   2, 1, 1, 0},
            {"20B9S100T1T",   makeList(rd(20, B),    rd(9,  S),   rd(100,  T),    rd(1,  T)),   2, 1, 1, 0},
            {"33T42B9T81B",   makeList(rd(33, T),    rd(42, B),   rd(9,    T),    rd(81, B)),   2, 2, 0, 0},
            {"28B56B13T123S", makeList(rd(28, B),    rd(56, B),   rd(13,   T),    rd(123,S)),   1, 2, 1, 0},
            {"92S8B8B32B",    makeList(rd(92, S),    rd(8,  B),   rd(8,    B),    rd(32, B)),   0, 3, 1, 0},
            {"92S8M8M32M",    makeList(rd(92, S),    rd(8,  M),   rd(8,    M),    rd(32, M)),   0, 0, 1, 3},

            {"2S88B7T8S9T9T84B100S2S4B3B", makeList(rd(2,S), rd(88,B), rd(7,T), rd(8,S), rd(9,T), rd(9,T), rd(84,B), rd(100,S), rd(2,S),  rd(4,B), rd(3,B)), 3, 4, 4, 0},
            {"2S88B7T8S9T9T84B3M100S2S4B3M3B", makeList(rd(2,S), rd(88,B), rd(7,T), rd(8,S), rd(9,T), rd(9,T), rd(84,B), rd(3,M), rd(100,S), rd(2,S),  rd(4,B), rd(3, M), rd(3,B)), 3, 4, 4, 2}
        };
    }

    @DataProvider(name="invalidReadStructures")
    public Object[][] invalidReadStructures() {
        return new Object[][]{
                {"",         new ArrayList<ReadDescriptor>()},
                {"0T",       makeList(rd(0,  T))},
                {"-1T",      makeList(rd(-1, T))},
                {"0S" ,      makeList(rd(0,  S))},
                {"0M" ,      makeList(rd(0,  M))},
                {"-1B",      makeList(rd(-1, B))},
                {"-1M",      makeList(rd(-1, M))},
                {"8C",       null},
                {"B5",       null},
                {"SS",       null},
                {"SM",       null},
                {"75TS",     null},
                {"8*T",      null},
                {"-66S1B",   makeList(rd(-66, S), rd(1, B))},
                {"-0T5B8C",  null},
                {"77T82B0S", makeList(rd(77, T), rd(82, B), rd(0, S))}
        };
    }

    @DataProvider(name="invalidReadStructuresFromList")
    public Object[][] invalidReadStructuresFromList() {
        int numTests = 0;
        for(final Object [] args : invalidReadStructures()) {
            if(args[1] != null) ++numTests;
        }

        final Object [][] outObjs = new Object[numTests][2];

        numTests = 0;
        for(final Object [] args : invalidReadStructures()) {
            if(args[1] != null) {
                outObjs[numTests++] = args;
            }
        }

        return outObjs;
    }

    @Test(dataProvider = "validReadStructures")
    public void testValidStructuresFromString(final String rsString, final List<ReadDescriptor> descriptors, final int numTemplates, final int numBarcodes, final int numSkips, final int numMolecularIndexes) {
        final ReadStructure readStructure = new ReadStructure(rsString);
        testReadStructure(readStructure, rsString, descriptors, numTemplates, numBarcodes, numSkips, numMolecularIndexes);
    }
    
    @Test(dataProvider = "validReadStructures")
    public void testValidStructuresFromList(final String rsString, final List<ReadDescriptor> descriptors, final int numTemplates, final int numBarcodes, final int numSkips, final int numMolecularIndexes) {
        final ReadStructure readStructure = new ReadStructure(descriptors);
        testReadStructure(readStructure, rsString, descriptors, numTemplates, numBarcodes, numSkips, numMolecularIndexes);
    }

    private void testReadStructure(final ReadStructure readStructure, final String structureString, final List<ReadDescriptor> descriptors, final int numTemplates, final int numBarcodes, final int numSkips, final int numMolecularIndexes) {
        Assert.assertEquals(readStructure.toString(), structureString);

        int totalCycles = 0;

        int tIndex = 0;
        int bIndex = 0;
        int sIndex = 0;
        int mIndex = 0;

        for(int i = 0; i < descriptors.size(); i++) {
            Assert.assertEquals(readStructure.descriptors.get(i), descriptors.get(i));
            switch(readStructure.descriptors.get(i).type) {
                case T:
                    Assert.assertEquals(i, readStructure.templates.getIndices()[tIndex++]);
                    break;
                case B:
                    Assert.assertEquals(i, readStructure.sampleBarcodes.getIndices()[bIndex++]);
                    break;
                case S:
                    Assert.assertEquals(i, readStructure.skips.getIndices()[sIndex++]);
                    break;
                case M:
                    Assert.assertEquals(i, readStructure.molecularBarcode.getIndices()[mIndex++]);
                    break;

                default:
                    Assert.fail("Unrecognized read type: " + readStructure.descriptors.get(i).type);
            }
            totalCycles += readStructure.descriptors.get(i).length;
        }

        Assert.assertEquals(readStructure.totalCycles,  totalCycles);
        Assert.assertEquals(readStructure.sampleBarcodes.length(),  numBarcodes);
        Assert.assertEquals(readStructure.templates.length(), numTemplates);
        Assert.assertEquals(readStructure.molecularBarcode.length(), numMolecularIndexes);
        Assert.assertEquals(readStructure.skips.length(), numSkips);

    }

    @Test(dataProvider = "invalidReadStructures", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidReadStructureFromString(final String rsString, final List<ReadDescriptor> descriptors) {
        final ReadStructure readStructure = new ReadStructure(rsString);
    }

    @Test(dataProvider = "invalidReadStructuresFromList", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidReadStructureFromList(final String rsString, final List<ReadDescriptor> descriptors) {
        final ReadStructure readStructure = new ReadStructure(descriptors);
    }

    @DataProvider(name="substructuresToReadStructureData")
    public Object [][] substructureToReadStructureData() {
        return new Object[][] {
            {new ReadStructure("10T10T"). templates,     "10T10T"    },
            {new ReadStructure("10T4M10T").templates,    "10T10T"    },
            {new ReadStructure("10T8B10T").    nonSkips, "10T8B10T"  },
            {new ReadStructure("10T8B5M10T").  nonSkips, "10T8B5M10T"},
            {new ReadStructure("8S10T8B8S10T").nonSkips, "10T8B10T"  },
            {new ReadStructure("10T8S8S10T").  skips,    "8S8S"      },
            {new ReadStructure("10T8S8S3M10T").skips,    "8S8S"      },
            {new ReadStructure("8B").sampleBarcodes,     "8B"        },
            {new ReadStructure("10T8S8M10T").molecularBarcode, "8M"  }
        };
    }

    @Test(dataProvider = "substructuresToReadStructureData")
    public void testSubstructureToReadStructure(final ReadStructure.Substructure substructure, final String outputRs) {
        Assert.assertEquals(substructure.toReadStructure().toString(), outputRs);
    }

    @DataProvider(name="substructureToReadStructureNegativeData")
    public Object[][] substructureToReadStructureNegativeData() {
        return new Object[][] {
            {new ReadStructure("10T").sampleBarcodes  },
            {new ReadStructure("10M").sampleBarcodes  },
            {new ReadStructure("10T").molecularBarcode},
            {new ReadStructure("10S").nonSkips        },
            {new ReadStructure("10M").skips           },
            {new ReadStructure("10S8B").templates     },
            {new ReadStructure("10S8B4M").templates   },

        };
    }

    @Test(dataProvider = "substructureToReadStructureNegativeData", expectedExceptions = IllegalArgumentException.class)
    public void testSubstructureToReadStructure(final ReadStructure.Substructure substructure) {
        substructure.toReadStructure().toString();
    }
}
