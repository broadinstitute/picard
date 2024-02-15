/*
 * The MIT License
 *
 * Copyright (c) 2022 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesTester (see getTester).
 */
public class MarkDuplicatesForFlowHelperTest {
    protected static String TEST_BASE_NAME = null;
    protected static File TEST_DATA_DIR = null;
    final static String   FLOW_ORDER = "TGCA";


    @BeforeClass
    public void setUp() {
        TEST_BASE_NAME = "MarkDuplicatesForFlowHelper";
        TEST_DATA_DIR = new File("testdata/picard/sam/MarkDuplicates");
    }

    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() {
        return new MarkDuplicatesTester();
    }

    private interface TesterModifier {
        void modify(AbstractMarkDuplicatesCommandLineProgramTester tester);
    }

    private static class TestRecordInfo {
        final int  length;
        final int  alignmentStart;
        final String  cigar;
        final boolean isDuplicate;
        final String startMod;
        final String endMod;

        TestRecordInfo(final int length, final int alignmentStart, final String cigar, final boolean isDuplicate,
                       final String startMod, final String endMod) {
            this.length = length;
            this.alignmentStart = alignmentStart;
            this.cigar = cigar;
            this.isDuplicate = isDuplicate;
            this.startMod = startMod;
            this.endMod = endMod;
        }
    }

    @DataProvider(name ="forFlowDataProvider")
    public Object[][] forFlowDataProvider() {
        return new Object[][] {
                // testUSE_END_IN_UNPAIRED_READS: End location is not significant
                {
                    null,
                    new TestRecordInfo[] {
                            new TestRecordInfo(76, 12, null, false, null, null),
                            new TestRecordInfo(74, 12, null, true, null, null)
                    },
                    new String[] { "USE_END_IN_UNPAIRED_READS=false" }, null
                },

                // testUSE_END_IN_UNPAIRED_READS: End location is significant
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, null, false, null, null),
                                new TestRecordInfo(74, 12, null, false, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Do not use clipped locations (meaning, use unclipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M", false, null, null),
                                new TestRecordInfo(74, 12, "74M", false, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=false" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Use clipped locations (meaning, use clipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M", false, null, null),
                                new TestRecordInfo(74, 12, "74M", true, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=true" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Use clipped locations (meaning, use clipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M1S", true, null, null),
                                new TestRecordInfo(78, 11, "78M", false, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=false", "USE_END_IN_UNPAIRED_READS=true" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Use clipped locations (meaning, use clipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M1S", false, null, null),
                                new TestRecordInfo(78, 11, "78M", false, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=true", "USE_END_IN_UNPAIRED_READS=true" }, null
                },

                // testFLOW_SKIP_FIRST_N_FLOWS: Do not use clipped locations (meaning, use unclipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", false, "ACGTT", "TTGCA"),
                                new TestRecordInfo(76, 12, "76M", true, "ACGGT", "TGGCA")
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "FLOW_SKIP_FIRST_N_FLOWS=0" }, null
                },

                // testFLOW_SKIP_FIRST_N_FLOWS: Do not use clipped locations (meaning, use unclipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", false, "ACGTT", "TTGCA"),
                                new TestRecordInfo(76, 12, "76M", false, "CCGGT", "TGGCA")
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "FLOW_SKIP_FIRST_N_FLOWS=3" }, null
                },

                // testFLOW_QUALITY_SUM_STRATEGY: normal sum
                {
                        DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", true, "AAAC", null),
                                new TestRecordInfo(76, 12, "76M", false, "AACC", null)
                        },
                        new String[] { "FLOW_DUP_STRATEGY=FLOW_QUALITY_SUM_STRATEGY" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("tp", new int[76]);
                                records[1].setAttribute("tp", new int[76]);
                                records[0].getBaseQualities()[0] = 25; // dip inside AAA
                                records[0].getBaseQualities()[2] = 25; // dip inside AAA

                            }
                        }
                },

                // testFLOW_QUALITY_SUM_STRATEGY: flow (homopolymer based) sum
                {
                        DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", false, "AAAC", null),
                                new TestRecordInfo(76, 12, "76M", true, "AACC", null)
                        },
                        new String[] { "FLOW_DUP_STRATEGY=FLOW_QUALITY_SUM_STRATEGY" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("tp", new int[76]);
                                records[1].setAttribute("tp", new int[76]);
                                records[0].getBaseQualities()[1] = 25; // dip inside AAA
                            }
                        }
                },
                // testFLOW_END_QUALITY_STRATEGY: flow (homopolymer based) minumum
                {
                        DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", true, "AAAC", null),
                                new TestRecordInfo(76, 12, "76M", false, "AACC", null)
                        },
                        new String[] { "FLOW_DUP_STRATEGY=FLOW_END_QUALITY_STRATEGY" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("tp", new int[76]);
                                records[1].setAttribute("tp", new int[76]);
                                records[0].getBaseQualities()[1] = 25; // dip inside AAA
                                records[1].getBaseQualities()[30] = 10;
                            }
                        }
                },

                // testUNPAIRED_END_UNCERTAINTY: End location is significant and uncertain, end sorted
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(74, 12, null, true, null, null),
                                new TestRecordInfo(84, 12, null, true, null, null),
                                new TestRecordInfo(94, 12, null, false, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "UNPAIRED_END_UNCERTAINTY=10" }, null
                },

                // testUNPAIRED_START_UNCERTAINTY: End location is significant and uncertain, end sorted
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(74, 12, null, false, null, null),
                                new TestRecordInfo(64, 22, null, true, null, null),
                                new TestRecordInfo(54, 32, null, true, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "UNPAIRED_START_UNCERTAINTY=10" }, null
                },

                // testUNPAIRED_END_UNCERTAINTY: End location is significant and uncertain, end not sorted
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(174, 12, null, true, null, null),
                                new TestRecordInfo(194, 12, null, false, null, null),
                                new TestRecordInfo(184, 12, null, true, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "UNPAIRED_END_UNCERTAINTY=10" }, null
                },

                // testUNPAIRED_END_UNCERTAINTY: End location is significant and uncertain, end not sorted, multiple non-dup
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(174, 12, null, true, null, null),
                                new TestRecordInfo(294, 12, null, false, null, null),
                                new TestRecordInfo(184, 12, null, false, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "UNPAIRED_END_UNCERTAINTY=10" }, null
                },
                // Barcode
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, null, false, null, null),
                                new TestRecordInfo(74, 12, null, true, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=false", "BARCODE_TAG=BC" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("BC", "A");
                                records[1].setAttribute("BC", "A");
                            }
                        }
                },
                // Barcode
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, null, false, null, null),
                                new TestRecordInfo(74, 12, null, false, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=false", "BARCODE_TAG=BC" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("BC", "A");
                                records[1].setAttribute("BC", "T");
                            }
                        }
                },
        };
    }

    @Test(dataProvider = "forFlowDataProvider")
    public void testForFlowMDCall(final DuplicateScoringStrategy.ScoringStrategy scoringStrategy, final TestRecordInfo[] recInfos, final String[] params, TesterModifier modifier) {

        // get tester, build records
        final AbstractMarkDuplicatesCommandLineProgramTester tester =
                scoringStrategy == null ? getTester() : new MarkDuplicatesTester(scoringStrategy);
        for ( final TestRecordInfo info : recInfos ) {
            tester.getSamRecordSetBuilder().setReadLength(info.length);
            if ( info.cigar != null ) {
                tester.addMappedFragment(0, info.alignmentStart, info.isDuplicate, info.cigar, 50);
            } else {
                tester.addMappedFragment(0, info.alignmentStart, info.isDuplicate, 50);
            }
        }

        // modify records
        final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
        for ( int i = 0 ; i < records.length ; i++ ) {
            final SAMRecord       rec = records[i];
            final TestRecordInfo  info = recInfos[i];

            if ( info.startMod != null ) {
                System.arraycopy(info.startMod.getBytes(), 0, rec.getReadBases(), 0, info.startMod.length());
            }
            if ( info.endMod != null ) {
                System.arraycopy(info.endMod.getBytes(), 0, rec.getReadBases(), rec.getReadBases().length - info.endMod.length(), info.endMod.length());
            }
        }

        // add parames, set flow order
        tester.addArg("FLOW_MODE=true");
        for ( final String param : params ) {
            tester.addArg(param);
        }
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);

        // further modify tester
        if ( modifier != null ) {
            modifier.modify(tester);
        }

        // run test
        tester.runTest();
    }

    @DataProvider(name ="getFlowSumOfBaseQualitiesDataProvider")
    public Object[][] getFlowSumOfBaseQualitiesDataProvider() {
        return new Object[][] {
                { "AAAA", new byte[] {50,50,50,50}, 30, 200 },
                { "AAAA", new byte[] {50,10,10,50}, 30, 200 },
                { "AAAA", new byte[] {20,10,10,20}, 30, 0 },
                { "AABBCC", new byte[] {50,50,10,10,40,40}, 30, 180 },
        };
    }

    @Test(dataProvider = "getFlowSumOfBaseQualitiesDataProvider")
    public void testGetFlowSumOfBaseQualities(final String bases, final byte[] quals, final int threshold, final int expectedScore) {

        // build record
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(bases.length());
        tester.addMappedFragment(0, 12, false, 50);

        // install bases and quals
        final SAMRecord rec = tester.getSamRecordSetBuilder().getRecords().iterator().next();
        System.arraycopy(bases.getBytes(), 0, rec.getReadBases(), 0,bases.length());
        System.arraycopy(quals, 0, rec.getBaseQualities(), 0, quals.length);

        // calculate score
        final int score = MarkDuplicatesForFlowHelper.getFlowSumOfBaseQualities(rec, threshold);
        Assert.assertEquals(score, expectedScore);
    }

    @DataProvider(name ="getFlowEndBaseQualitiesDataProvider")
    public Object[][] getFlowEndBaseQualitiesDataProvider() {
        return new Object[][] {
                { "AAAA", new byte[] {50,50,50,50}, 2, 50 },
                { "AAAA", new byte[] {50,10,10,50}, 4, 10 },
                { "ACCA", new byte[] {20,10,10,20}, 1, 20 },
                { "AABBCC", new byte[] {50,50,10,10,40,40}, 30, 10 },
        };
    }

    @Test(dataProvider = "getFlowEndBaseQualitiesDataProvider")
    public void testGetFlowEndBaseQualities(final String bases, final byte[] quals, final int threshold, final int expectedScore) {

        // build record
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(bases.length());
        tester.addMappedFragment(0, 12, false, 50);

        // install bases and quals
        final SAMRecord rec = tester.getSamRecordSetBuilder().getRecords().iterator().next();
        System.arraycopy(bases.getBytes(), 0, rec.getReadBases(), 0,bases.length());
        System.arraycopy(quals, 0, rec.getBaseQualities(), 0, quals.length);

        // calculate score
        final int score = MarkDuplicatesForFlowHelper.getFlowSumOfBaseQualitiesNearEnds(rec, threshold);
        Assert.assertEquals(score, expectedScore);
    }


}
