package picard.vcf;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Class that holds detail metrics about Genotype Concordance
 *
 * @author George Grant
 */
public class GenotypeConcordanceDetailMetrics extends MetricBase {
    /** The type of event SNP/Indel */
    public String EVENT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE;

    /** The state of the 'truth' sample (i.e. HET, HOMREF...) */
    public String TRUTH_STATE;

    /** The state of the 'call' sample (i.e. HET, HOMREF...) */
    public String CALL_STATE;

    /** Whether the alternate alleles agree for between truth and call genotypes */
    public boolean ALT_ALLELES_AGREE;

    /** The number of events of type TRUTH_STATE and CALL_STATE for the EVENT_TYPE and SAMPLEs */
    public long COUNT;
}
