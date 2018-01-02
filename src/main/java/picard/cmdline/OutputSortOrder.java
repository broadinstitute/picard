package picard.cmdline;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import org.broadinstitute.barclay.argparser.CommandLineParser;

/**
 * A OutputSortOrder class intended to expose the various options available as inputs to CommandLinePrograms
 * that need to control the output OutputSortOrder of the SAM files they generate.
 *
 * In particular this enables to add a description and also to not expose "duplicate", "unsorted" and "unknown"
 * as they are not appropriate values to sort a file into.
 */
public enum OutputSortOrder implements CommandLineParser.ClpEnum {
    queryname("Sorts according to the readname. This will place read-pairs and other derived reads (secondary and " +
            "supplementary) adjacent to each other. Note that the readnames are compared lexicographically, even though " +
            "they may include numbers. In paired reads, Read1 sorts before Read2."),
    coordinate("Sorts primarily according to the SEQ and POS fields of the record. The sequence will sorted according to " +
            "the order in the sequence dictionary, taken from from the header of the file. Within each reference sequence, the " +
            "reads are sorted by the position. Unmapped reads whose mates are mapped will be placed near their mates. " +
            "Unmapped read-pairs are placed after all the mapped reads and their mates.");
    private String description;

    OutputSortOrder(String description) {
        this.description = description;
    }

    public SortOrder getSortOrder() {
        return SortOrder.valueOf(this.name());
    }

    @Override
    public String getHelpDoc() {
        return description;
    }
}