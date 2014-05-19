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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Very basic test for scatter functionality in IntervalListTools
 */
public class IntervalListScattererTest {
    private static final IntervalList LIST_TO_SCATTER;

    static {
        LIST_TO_SCATTER = IntervalList.fromFile(new File("testdata/picard/util/scatterable.interval_list"));
        Assert.assertEquals(LIST_TO_SCATTER.getUniqueBaseCount(), 200, "Wrong unique base count");
    }

    private static class Testcase {
        final IntervalList source;
        final List<IntervalList> expectedScatter;
        final int scatterWidth;
        final IntervalListScatterer.Mode mode;

        @Override
        public String toString() {
            return "Testcase{" +
                    "scatterWidth=" + scatterWidth +
                    ", mode=" + mode +
                    '}';
        }

        private Testcase(final IntervalList source, final int scatterWidth, final IntervalListScatterer.Mode mode, final List<IntervalList> expectedScatter) {
            this.source = source;
            this.expectedScatter = expectedScatter;
            this.scatterWidth = scatterWidth;
            this.mode = mode;
        }
    }

    private static final List<Testcase> testcases = new ArrayList<Testcase>();

    static {
        testcases.add(new Testcase(
                LIST_TO_SCATTER, 2, IntervalListScatterer.Mode.INTERVAL_SUBDIVISION,
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
        
        testcases.add(new Testcase(
                LIST_TO_SCATTER, 4, IntervalListScatterer.Mode.INTERVAL_SUBDIVISION,
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

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 5, IntervalListScatterer.Mode.INTERVAL_SUBDIVISION,
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

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 6, IntervalListScatterer.Mode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
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

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 2, IntervalListScatterer.Mode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
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

        testcases.add(new Testcase(
                LIST_TO_SCATTER, 1, IntervalListScatterer.Mode.BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
                Arrays.asList(
                        composeIntervalList(LIST_TO_SCATTER, "1",
                                30000, 30098,
                                30100, 30150,
                                30200, 30249
                        )
                )
        ));
    }

    @DataProvider
    public Object[][] testScatterTestcases() {
        final Object[][] objects = new Object[testcases.size()][];
        for (int i = 0; i < objects.length; i++) {
            objects[i] = new Object[]{testcases.get(i)};
        }
        return objects;
    }

    @Test(dataProvider = "testScatterTestcases")
    public void testScatter(final Testcase tc) {
        final IntervalListScatterer scatterer = new IntervalListScatterer(tc.mode);
        final List<IntervalList> scatter = scatterer.scatter(tc.source, tc.scatterWidth);
        Assert.assertEquals(scatter, tc.expectedScatter);
    }
    
    private static IntervalList composeIntervalList(final IntervalList source, final String chromosome, final int... segmentsByPair) {
        final IntervalList intervals = new IntervalList(source.getHeader());
        for (int i = 0; i < segmentsByPair.length; i += 2) {
            final Interval parentInterval = lookupIntervalContainingLocus(source, chromosome, segmentsByPair[i]);
            intervals.add(new Interval("1", segmentsByPair[i], segmentsByPair[i + 1], parentInterval.isNegativeStrand(), parentInterval.getName()));
        }
        return intervals;
    }
    
    private static Interval lookupIntervalContainingLocus(final IntervalList source, final String chromosome, final int startIndex) {
        for (final Interval interval : source) {
            if (interval.getSequence().equals(chromosome) && startIndex >= interval.getStart() && startIndex <= interval.getEnd()) {
                return interval;
            }
        }
        throw new NoSuchElementException();
    }
}
