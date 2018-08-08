/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package picard.util;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.IntervalList.IntervalListScatterMode;
import picard.util.IntervalList.IntervalListScatterer;

import java.io.File;
import java.util.*;

/**
 * Basic test for scatter functionality in IntervalListTools
 */
public class IntervalListScattererTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/util");

    private static final File INTERVAL_FILE = new File(TEST_DATA_DIR, "scatterable.interval_list");
    private static final IntervalList LIST_TO_SCATTER = IntervalList.fromFile(INTERVAL_FILE);

    private static final File INTERVAL_WITH_OVERFLOW_FILE = new File(TEST_DATA_DIR, "scatterable_with_overflow.interval_list");
    private static final IntervalList LIST_TO_SCATTER_WITH_OVERFLOW = IntervalList.fromFile(INTERVAL_WITH_OVERFLOW_FILE);

    private static final File SCATTER_MANY_INTERVAL_FILE = new File(TEST_DATA_DIR, "scatterable_many_intervals.interval_list");
    private static final IntervalList LIST_TO_SCATTER_MANY = IntervalList.fromFile(SCATTER_MANY_INTERVAL_FILE);

    static class Testcase {
        final File file;
        final IntervalList source;
        final List<IntervalList> expectedScatter;
        final int scatterWidth;
        final IntervalListScatterMode mode;

        @Override
        public String toString() {
            return "Testcase{" +
                    "scatterWidth=" + scatterWidth +
                    ", mode=" + mode +
                    '}';
        }

        private Testcase(final File file, final int scatterWidth, final IntervalListScatterMode mode, final List<IntervalList> expectedScatter) {
            this.source = IntervalList.fromFile(file);
            this.file = file;
            this.expectedScatter = expectedScatter;
            this.scatterWidth = scatterWidth;
            this.mode = mode;
        }
    }

    @DataProvider
    public static Iterator<Object[]> testScatterTestcases() {
        final List<Testcase> testCases = new ArrayList<>();
        Assert.assertEquals(LIST_TO_SCATTER.getUniqueBaseCount(), 200, "Wrong unique base count");
        Assert.assertEquals(LIST_TO_SCATTER_MANY.getUniqueBaseCount(), 32 * 2, "Wrong unique base count");
        Assert.assertEquals(LIST_TO_SCATTER_MANY.getIntervals().size(), 32, "Wrong unique interval count");

        testCases.add(new Testcase(
                INTERVAL_FILE, 2, IntervalListScatterMode.INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30100
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30101, 30150,
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 4, IntervalListScatterMode.INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30049
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30050, 30098,
                                30100, 30100
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30101, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 5, IntervalListScatterMode.INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30039
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30040, 30079
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30080, 30098,
                                30100, 30120
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30121, 30150,
                                30200, 30209
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30210, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 6, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30100, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 2, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249)
                )
        ));
        testCases.add(new Testcase(
                INTERVAL_FILE, 6, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30100, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 2, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30100, 30150,
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 1, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                Collections.singletonList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30150,
                                30200, 30249
                        )
                )
        ));
        testCases.add(new Testcase(
                INTERVAL_FILE, 6, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30100, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 2, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30150
                        ),
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_FILE, 1, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW,
                Collections.singletonList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30150,
                                30200, 30249
                        )
                )
        ));

        testCases.add(new Testcase(
                INTERVAL_WITH_OVERFLOW_FILE, 7, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30000, 30000,
                                30100, 30109,
                                30200, 30209

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30300, 30309

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30400, 30409

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30500, 30509

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30600, 30609

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30700, 30709

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30800, 30808

                        )
                )));

        testCases.add(new Testcase(
                INTERVAL_WITH_OVERFLOW_FILE, 7, IntervalListScatterMode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30000, 30000,
                                30100, 30109
                        ),
                        composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30200, 30209
                        ),
                        composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30300, 30309
                        ),
                        composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30400, 30409

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30500, 30509

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30600, 30609

                        ), composeIntervalList(LIST_TO_SCATTER_WITH_OVERFLOW, "1",
                                30700, 30709,
                                30800, 30808
                        )
                )));

        final IntervalList full = new IntervalList(LIST_TO_SCATTER_MANY.getHeader());
        full.add(new Interval("1", 30000, 30000 + 32 * 2 - 1));

        testCases.add(new Testcase(
                SCATTER_MANY_INTERVAL_FILE, 1, IntervalListScatterMode.SCATTER_BY_INTERVAL_COUNT,
                Collections.singletonList(IntervalList.overlaps(LIST_TO_SCATTER_MANY, full))));

        final IntervalList half = new IntervalList(LIST_TO_SCATTER_MANY.getHeader());
        half.add(new Interval("1", 30000, 30000 + 16 * 2 - 1));
        testCases.add(new Testcase(
                SCATTER_MANY_INTERVAL_FILE, 2, IntervalListScatterMode.SCATTER_BY_INTERVAL_COUNT,
                Arrays.asList(IntervalList.overlaps(LIST_TO_SCATTER_MANY, half),
                        IntervalList.overlaps(LIST_TO_SCATTER_MANY, IntervalList.invert(half)))));

        final IntervalList third = new IntervalList(LIST_TO_SCATTER_MANY.getHeader());
        third.add(new Interval("1", 30000, 30000 + 10 * 2 - 1));
        final IntervalList secondThird = new IntervalList(LIST_TO_SCATTER_MANY.getHeader());
        secondThird.add(new Interval("1", 30000 + 10 * 2, 30000 + 20 * 2 - 1));
        testCases.add(new Testcase(
                SCATTER_MANY_INTERVAL_FILE, 3, IntervalListScatterMode.SCATTER_BY_INTERVAL_COUNT,
                Arrays.asList(IntervalList.overlaps(LIST_TO_SCATTER_MANY, third),

                        IntervalList.overlaps(LIST_TO_SCATTER_MANY, secondThird),

                        IntervalList.overlaps(LIST_TO_SCATTER_MANY, IntervalList.invert(IntervalList.concatenate(Arrays.asList(
                                IntervalList.overlaps(LIST_TO_SCATTER_MANY, third),
                                IntervalList.overlaps(LIST_TO_SCATTER_MANY, secondThird))
                        ))))));

        return testCases.stream().map(tc -> new Object[]{tc}).iterator();
    }

    @Test(dataProvider = "testScatterTestcases")
    public void testScatter(final Testcase tc) {
        final IntervalListScatterer scatterer = tc.mode.make();
        final List<IntervalList> scatter = scatterer.scatter(tc.source, tc.scatterWidth);
        Assert.assertEquals(scatter.size(), tc.expectedScatter.size());

        for (int i = 0; i < scatter.size(); i++) {
            Assert.assertEquals(scatter.get(i).getIntervals(), tc.expectedScatter.get(i).getIntervals(), "Problem with the " + i + " scatter");
        }
    }

    private static IntervalList composeIntervalList(final IntervalList source, final String chromosome, final int... segmentsByPair) {
        final IntervalList intervals = new IntervalList(source.getHeader());
        for (int i = 0; i < segmentsByPair.length; i += 2) {
            final Interval parentInterval = lookupIntervalContainingLocus(source, chromosome, segmentsByPair[i]);
            intervals.add(new Interval(chromosome, segmentsByPair[i], segmentsByPair[i + 1], parentInterval.isNegativeStrand(), parentInterval.getName()));
        }
        return intervals;
    }

    private static Interval lookupIntervalContainingLocus(final IntervalList source, final String chromosome, final int startIndex) {
        for (final Interval interval : source) {
            if (interval.getContig().equals(chromosome) && startIndex >= interval.getStart() && startIndex <= interval.getEnd()) {
                return interval;
            }
        }
        throw new NoSuchElementException();
    }
}
