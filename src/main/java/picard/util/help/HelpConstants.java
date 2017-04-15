package picard.util.help;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public final class HelpConstants {

    private HelpConstants() {}

    /**
     * Definition of the group names / descriptions for documentation/help purposes.
     */
    public final static String DOC_CAT_FASTA = "Fasta File Tools";
    public final static String DOC_CAT_FASTA_SUMMARY = "Tools for analysis and manipulation of files in fasta format";

    public final static String DOC_CAT_INTERVALS = "Interval Tools";
    public final static String DOC_CAT_INTERVALS_SUMMARY = "Tools for processing intervals and associated overlapping records";

    public final static String DOC_CAT_READS = "SAM/BAM/CRAM Tools";
    public final static String DOC_CAT_READS_SUMMARY = "Tools for manipulating read-level data (SAM/BAM/CRAM)";

    public final static String DOC_CAT_VARIANT = "VCF Tools";
    public final static String DOC_CAT_VARIANT_SUMMARY = "Tools for manipulating variants and associated metadata";

    public final static String DOC_CAT_TEST = "Test Tools";
    public final static String DOC_CAT_TEST_SUMMARY = "Tools for internal test purposes";

    /**
     * List of "supercategory" values used for doc purposes. Every doc group name can/should be put into
     * one of the following supercategories.
     */
    public final static String DOC_SUPERCAT_TOOLS = "tools";
    public final static String DOC_SUPERCAT_UTILITIES = "utilities";
    public final static String DOC_SUPERCAT_EXCLUDE = "exclude";

    private static Map<String, String> groupToSuperCategory;

    public static Map<String, String> getSuperCategoryMap() {
        if (groupToSuperCategory == null) {

            // do this only on demand since we only need it during docgen
            groupToSuperCategory = new HashMap<>();

            // supercat Tools
            groupToSuperCategory.put(DOC_CAT_FASTA, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_READS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT, DOC_SUPERCAT_TOOLS);

            // supercat Exclude
            groupToSuperCategory.put(DOC_CAT_TEST, DOC_SUPERCAT_EXCLUDE);
        }
        return groupToSuperCategory;
    }

    /**
     * Given a group name, return a supercategory string for use by the online doc system to determine which
     * supercateogry the group is in. The strings returned by this method should match those used in the
     * corresponding help template.
     *
     * @param groupName
     * @return supercategory string corresponding to {@code groupName} for use in for use in determining
     * which supercategory the group is in for online doc purposes.
     */
    public static String getSuperCategoryProperty(final String groupName) {
        return getSuperCategoryMap().getOrDefault(groupName, "other");
    }

}