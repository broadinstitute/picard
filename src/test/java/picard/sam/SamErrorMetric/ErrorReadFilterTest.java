package picard.sam.SamErrorMetric;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ErrorReadFilterTest {

    @DataProvider
    Object[][] provideForTestBooleanFilter() {
        return new Object[][] {
                {
                        new boolean[] {false},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.Equal},
                        new boolean[] {false},
                        new boolean[] {true},
                        true
                },
                {
                        new boolean[] {false},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.Equal},
                        new boolean[] {true},
                        new boolean[] {false},
                        false
                },
                {
                        new boolean[] {false},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.NotEqual},
                        new boolean[] {true},
                        new boolean[] {true},
                        true
                },
                {
                        new boolean[] {true, true},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.Equal, SamErrorReadFilter.Comparator.NotEqual},
                        new boolean[] {true, false},
                        new boolean[] {true, true},
                        true
                },
                {
                        new boolean[] {true, true},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.Equal, SamErrorReadFilter.Comparator.Equal},
                        new boolean[] {true, false},
                        new boolean[] {true, false},
                        false
                },
        };
    }

    @Test(dataProvider = "provideForTestBooleanFilter")
    void testBooleanFilter(final boolean[] values, final SamErrorReadFilter.Comparator[] comparators, final boolean[] thresholds, final boolean[] satisfied, final boolean overallSatisfied) {
        Map<String, ArrayList<SamErrorReadFilterCriterion>> criteria = new HashMap<>();

        for(int i = 0; i < values.length; i++) {
            criteria.put(String.valueOf(i), new ArrayList<>(Collections.singletonList(new BooleanSamErrorReadFilterCriterion(comparators[i], thresholds[i]))));
        }

        SamErrorReadFilter filter = new SamErrorReadFilter("testFilter", criteria);

        for(int valueIndex = 0; valueIndex < values.length; ) {
            filter.processValue(String.valueOf(valueIndex), values[valueIndex]);
            for(int suffixCriteriaIndex = 0; valueIndex < values.length && suffixCriteriaIndex < filter.criteria.get(String.valueOf(valueIndex)).size(); suffixCriteriaIndex++) {
                Assert.assertEquals(filter.criteria.get(String.valueOf(valueIndex)).get(suffixCriteriaIndex).isSatisifed(), satisfied[valueIndex]);
                valueIndex++;
            }
        }

        Assert.assertEquals(filter.isSatisfied(), overallSatisfied);
    }

    @DataProvider
    Object[][] provideForTestNumericFilter() {
        return new Object[][] {
                {
                        new int[] {0},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.Equal},
                        new int[] {0},
                        new boolean[] {true},
                        true
                },
                {
                        new int[] {0},
                        new SamErrorReadFilter.Comparator[] {SamErrorReadFilter.Comparator.Equal},
                        new int[] {1},
                        new boolean[] {false},
                        false
                },
        };
    }

    @Test(dataProvider = "provideForTestNumericFilter")
    void testNumericFilter(final int[] values, final SamErrorReadFilter.Comparator[] comparators, final int[] thresholds, final boolean[] satisfied, final boolean overallSatisfied) {
        Map<String, ArrayList<SamErrorReadFilterCriterion>> criteria = new HashMap<>();

        for(int i = 0; i < values.length; i++) {
            criteria.put(String.valueOf(i), new ArrayList<>(Collections.singletonList(new NumericSamErrorReadFilterCriterion(comparators[i], thresholds[i]))));
        }

        SamErrorReadFilter filter = new SamErrorReadFilter("testFilter", criteria);

        for(int valueIndex = 0; valueIndex < values.length; ) {
            filter.processValue(String.valueOf(valueIndex), values[valueIndex]);
            for(int suffixCriteriaIndex = 0; valueIndex < values.length && suffixCriteriaIndex < filter.criteria.get(String.valueOf(valueIndex)).size(); suffixCriteriaIndex++) {
                Assert.assertEquals(filter.criteria.get(String.valueOf(valueIndex)).get(suffixCriteriaIndex).isSatisifed(), satisfied[valueIndex]);
                valueIndex++;
            }
        }

        Assert.assertEquals(filter.isSatisfied(), overallSatisfied);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    void testIllegalComparator() {
        Map<String, ArrayList<SamErrorReadFilterCriterion>> criteria = new HashMap<>();
        criteria.put("test", new ArrayList<>(Collections.singletonList(new BooleanSamErrorReadFilterCriterion(SamErrorReadFilter.Comparator.Greater, true))));

        SamErrorReadFilter filter = new SamErrorReadFilter("testFilter", criteria);

        // Expect exception here
        filter.processValue("test", true);

        Assert.fail();
    }

    private Map<String, ArrayList<SamErrorReadFilterCriterion>> createCriteriaMap(String[] suffixes, SamErrorReadFilterCriterion[] criteria) {
        HashMap<String, ArrayList<SamErrorReadFilterCriterion>> criteriaMap = new HashMap<>();
        for(int i = 0; i < suffixes.length; i++) {
            if (criteriaMap.containsKey(suffixes[i])) {
                criteriaMap.get(suffixes[i]).add(criteria[i]);
            }
            else {
                criteriaMap.put(suffixes[i], new ArrayList<>(Collections.singletonList(criteria[i])));
            }
        }
        return criteriaMap;
    }

    @DataProvider
    Object[][] provideForTestFilterFileParsing() {
        return new Object[][] {
                {
                        new String[] {
                                "testname",
                                "criterionname\tint\t>=\t5"
                        },
                        new SamErrorReadFilter(
                                "testname",
                                createCriteriaMap(
                                        new String[] {"criterionname"},
                                        new SamErrorReadFilterCriterion[] {new NumericSamErrorReadFilterCriterion(SamErrorReadFilter.Comparator.GreaterOrEqual, 5)}
                                ))
                },
                {
                        new String[] {
                                "#",
                                "",
                                "  testname  ",
                                "",
                                "#",
                                " criterionname \t int \t >= \t 5 ",
                                "",
                                "#",
                                ""
                        },
                        new SamErrorReadFilter(
                                "testname",
                                createCriteriaMap(
                                        new String[] {"criterionname"},
                                        new SamErrorReadFilterCriterion[] {new NumericSamErrorReadFilterCriterion(SamErrorReadFilter.Comparator.GreaterOrEqual, 5)}
                                ))
                },
                {
                        new String[] {
                                "testname",
                                "criterionname\tint\t>>\t5"
                        },
                        new SamErrorReadFilter(
                                "testname",
                                createCriteriaMap(new String[] {}, new SamErrorReadFilterCriterion[] {})
                        )
                },
                {
                        new String[] {
                                "testname",
                                "criterionname\tint\t>\ta"
                        },
                        new SamErrorReadFilter(
                                "testname",
                                createCriteriaMap(new String[] {}, new SamErrorReadFilterCriterion[] {})
                        )
                },
        };
    }

    @Test(dataProvider = "provideForTestFilterFileParsing")
    void testFilterFileParsing(String[] lines, SamErrorReadFilter expected) throws IOException {
        final File file = File.createTempFile("test", ".tmp");
        file.deleteOnExit();
        FileWriter fileWriter = new FileWriter(file);
        for(String line : lines) {
            fileWriter.write(line + "\n");
        }
        fileWriter.close();

        SamErrorReadFilter filter = SamErrorReadFilter.fromFile(file);

        if(expected == null)
        {
            Assert.assertNull(filter);
        }
        else {
            Assert.assertEquals(filter.getName(), expected.getName());
            Assert.assertEquals(filter.criteria.size(), expected.criteria.size());

            for (String criterionKey : expected.criteria.keySet()) {
                for(int suffixCriterionIndex = 0; suffixCriterionIndex < expected.criteria.get(criterionKey).size(); suffixCriterionIndex++) {
                    SamErrorReadFilterCriterion criterion = filter.criteria.get(criterionKey).get(suffixCriterionIndex);
                    SamErrorReadFilterCriterion expectedCriterion = expected.criteria.get(criterionKey).get(suffixCriterionIndex);
                    Assert.assertEquals(criterion.value, expectedCriterion.value);
                    Assert.assertEquals(criterion.comparator, expectedCriterion.comparator);
                }
            }
        }
    }
}
