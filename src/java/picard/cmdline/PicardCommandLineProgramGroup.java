package picard.cmdline;

/**
 * A set of program groups to which Picard command line programs can belong.
 * @author Nils Homer
 */
public enum PicardCommandLineProgramGroup implements CommandLineProgramGroup {

    SamOrBam {
        public String getName() { return "SAM/BAM"; }
        public String getDescription() { return "Tools for manipulating SAM, BAM, or related data."; }
    },
    VcfOrBcf {
        public String getName() { return "VCF/BCF"; }
        public String getDescription() { return "Tools for manipulating VCF, BCF, or related data."; }
    },
    Fasta {
        public String getName() { return "Fasta"; }
        public String getDescription() { return "Tools for manipulating FASTA, or related data."; }
    },
    Metrics {
        public String getName() { return "Metrics"; }
        public String getDescription() { return "Tools for reporting metrics on various data types."; }
    },
    Intervals {
        public String getName() { return "Interval Tools"; }
        public String getDescription() { return "Tools for manipulating Picard interval lists."; }
    },
    Illumina {
        public String getName() { return "Illumina Tools"; }
        public String getDescription() { return "Tools for manipulating data specific to Illumina sequencers."; }
    },
    None {
        public String getName() { return "Miscellaneous Tools"; }
        public String getDescription() { return "A miscellaneous set of uncategorized tools."; }
    };
}
