package picard.sam.SamErrorMetric;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SamErrorReadFilter {
    public enum Comparator {
        Equal,
        NotEqual,
        Smaller,
        SmallerOrEqual,
        Greater,
        GreaterOrEqual;

        public static Comparator fromString(final String str) {
            switch (str) {
                case "=":
                case "==": return Comparator.Equal;
                case "!=": return Comparator.NotEqual;
                case "<": return Comparator.Smaller;
                case "<=": return Comparator.SmallerOrEqual;
                case ">": return Comparator.Greater;
                case ">=": return Comparator.GreaterOrEqual;
            }
            throw new IllegalArgumentException("Cannot convert given string to Comparator.");
        }
    }

    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    private final String name;

    @VisibleForTesting
    protected final Map<String, SamErrorReadFilterCriterion> criteria;
    private Map<String, Integer> filteredReads = new HashMap<>();

    /**
     * Creates a new SamErrorReadFilter.
     * @param name The name of the filter, which will determine the output file name
     * @param criteria A Map of criteria with their corresponding suffixes as keys
     */
    public SamErrorReadFilter(final String name, final Map<String, SamErrorReadFilterCriterion> criteria) {
        this.name = name;
        this.criteria = criteria;
    }

    public final String getName() {
        return name;
    }

    /**
     * Creates a new SamErrorReadFilter with an empty list of criteria.
     * @param name The name of the filter, which will determine the output file name
     */
    private SamErrorReadFilter(final String name) {
        this.name = name;
        this.criteria = new HashMap<>();
    }

    /**
     * Parses a SamErrorReadFilter from a file.
     *
     * @param file the file containing the filter
     * @return a SamErrorReadFilter object
     */
    public static SamErrorReadFilter fromFile(final File file) {
        SamErrorReadFilter filter = null;
        try (final BufferedReader reader = IOUtil.openFileForBufferedReading(file.toPath())) {
            String line = reader.readLine();
            // Read the filter name first
            while (line != null) {
                if (line.isEmpty() || line.startsWith("#")) {
                    line = reader.readLine();
                    continue;
                }
                filter = new SamErrorReadFilter(line.trim());
                break;
            }

            if(line == null) {
                log.warn("Empty filter input file: " + file.toPath().toUri());
                return null;
            }

            // Then read the criteria line by line
            line = reader.readLine();
            while (line != null) {
                if (line.isEmpty() || line.startsWith("#")) {
                    line = reader.readLine();
                    continue;
                }

                final char separator = '\t';
                final String[] splitUnits = new String[ReadBaseStratification.Stratifier.values().length + 1];
                final int numberOfTerms = ParsingUtils.split(line, splitUnits, separator, false);
                if(numberOfTerms == 4) {
                    try {
                        SamErrorReadFilterCriterion criterion = null;
                        Comparator comparator = Comparator.fromString(splitUnits[2].trim());

                        switch (splitUnits[1].trim()) {
                            case "boolean":
                                if(!splitUnits[3].trim().toLowerCase().equals(String.valueOf(true)) && !splitUnits[3].trim().toLowerCase().equals(String.valueOf(false)))
                                    throw new NumberFormatException();
                                criterion = new BooleanSamErrorReadFilterCriterion(comparator, Boolean.valueOf(splitUnits[3].trim())); break;
                            case "int": criterion = new NumericSamErrorReadFilterCriterion(comparator, Integer.valueOf(splitUnits[3].trim())); break;
                            default: log.warn("Invalid criterion type in line \""+ line + "\" in filter input file: " + file.toPath().toUri()); continue;
                        }
                        filter.criteria.put(splitUnits[0].trim(), criterion);
                    } catch (IllegalArgumentException ex) {
                        // Also catches NumberFormatException, which is a subclass of IllegalArgumentException
                        log.warn("Invalid value in line \""+ line + "\" in filter input file: " + file.toPath().toUri() + ". Skipping this criterion.");
                    }
                }
                else {
                    log.warn("Invalid line in filter input file: " + file.toPath().toUri());
                }
                line = reader.readLine();
            }
        } catch (final IOException e) {
            throw new SAMException(String.format("Failed to close file %s after reading", file.toPath().toUri().toString()));
        }
        return filter;
    }

    /**
     * Add a read to an internal container to be able to keep track of the number of occurrences.
     * @param readId The ID of the read to be added.
     */
    public void addReadById(String readId) {
        if (filteredReads.containsKey(readId)) {
            filteredReads.put(readId, filteredReads.get(readId) + 1);
        }
        else {
            filteredReads.put(readId, 1);
        }
    }

    /**
     * Returns all filters that have been added by the addReadById method with its corresponding number of occurrences.
     */
    public final Map<String, Integer> getFilteredReads() {
        return filteredReads;
    }

    /**
     * Processes a value for a specific suffix. If the filter does contain criteria for that suffix and the value
     * satisfies that criterion, the isSatisfied flag for that criterion is automatically set.
     * @param suffix The suffix for the specific criterion
     * @param stratus The value to be checked
     */
    public void processValue(final String suffix, final Object stratus) {
        if(criteria.containsKey(suffix)) {
            SamErrorReadFilterCriterion criterion = criteria.get(suffix);

            criterion.checkCriterion(stratus);
        }
    }

    /**
     * Returns whether or not all criteria have been satisfied after the last call the filter has been reset.
     */
    public boolean isSatisfied() {
        return criteria.values().stream().allMatch(SamErrorReadFilterCriterion::isSatisifed);
    }

    /**
     * Resets the filter. This must be called each time a new RecordAndOffset is considered.
     */
    public void reset() {
        criteria.forEach((suffix, criterion) -> criterion.reset());
    }
}
