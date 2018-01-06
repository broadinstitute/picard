package picard.util.help;

import java.util.HashMap;
import java.util.Map;

public final class HelpConstants {

    private HelpConstants() {}

    /**
     * Definition of the group names / descriptions for documentation/help purposes.
     */
    public final static String DOC_CAT_BASE_CALLING = "Base Calling";
    public final static String DOC_CAT_BASE_CALLING_SUMMARY = "Tools that process sequencing machine data, e.g. Illumina base calls, and detect sequencing level attributes, e.g. adapters";

    public final static String DOC_CAT_DIAGNOSTICS_AND_QC = "Diagnostics and Quality Control";
    public final static String DOC_CAT_DIAGNOSTICS_AND_QC_SUMMARY = "Tools that collect sequencing quality related and comparative metrics";

    public final static String DOC_CAT_INTERVALS_MANIPULATION = "Intervals Manipulation";
    public final static String DOC_CAT_INTERVALS_MANIPULATION_SUMMARY = "Tools that process genomic intervals in various formats";

    public final static String DOC_CAT_OTHER = "Other";
    public final static String DOC_CAT_OTHER_SUMMARY = "Miscellaneous tools, e.g. those that aid in data streaming";

    public final static String DOC_CAT_READ_DATA_MANIPULATION = "Read Data Manipulation";
    public final static String DOC_CAT_READ_DATA_MANIPULATION_SUMMARY = "Tools that manipulate read data in SAM, BAM or CRAM format";

    public final static String DOC_CAT_REFERENCE = "Reference";
    public final static String DOC_CAT_REFERENCE_SUMMARY = "Tools that analyze and manipulate FASTA format references";

    public final static String DOC_CAT_TEST = "Test Tools";
    public final static String DOC_CAT_TEST_SUMMARY = "Tools for internal test purposes";

    public final static String DOC_CAT_VARIANT_FILTERING = "Variant Filtering";
    public final static String DOC_CAT_VARIANT_FILTERING_SUMMARY = "Tools that filter variants by annotating the FILTER column";

    public final static String DOC_CAT_VARIANT_EVALUATION = "Variant Evaluation and Refinement";
    public final static String DOC_CAT_VARIANT_EVALUATION_SUMMARY = "Tools that evaluate and refine variant calls, e.g. with annotations not offered by the engine";

    public final static String DOC_CAT_VARIANT_MANIPULATION = "Variant Manipulation";
    public final static String DOC_CAT_VARIANT_MANIPULATION_SUMMARY = "Tools that manipulate variant call format (VCF) data";

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
            groupToSuperCategory.put(DOC_CAT_BASE_CALLING, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_DIAGNOSTICS_AND_QC, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_INTERVALS_MANIPULATION, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_OTHER, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_READ_DATA_MANIPULATION, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_REFERENCE, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_READ_DATA_MANIPULATION, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT_FILTERING, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT_EVALUATION, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT_MANIPULATION, DOC_SUPERCAT_TOOLS);

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