package picard.analysis.replicates;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.MergeSamFiles;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.testng.Assert.assertEquals;

/**
 * Created by farjoun on 6/24/15.
 */
public class CollectIndependentReplicatesMetricTest {
    private final static File testdir = new File("testdata/picard/independent_replicates");
    private final static File bamOutDir = IOUtil.createTempDir("convertSamToBam", "dir");

    private final static Map<String, String> sams = new ImmutableMap.Builder<String, String>()
            .put("twoPairs", "twopairs.sam")
            .put("twoPairsWithUMIs", "twopairsWithUMIs.sam")
            .put("twoPairsWithBadUMIs", "twopairsWithBadUMIs.sam")
            .put("aTriple", "aTriple.sam")
            .put("aTripleWithUMIs", "aTripleWithUMIs.sam")
            .put("multipleContigs", "multipleContigs.sam")
            .put("twopairsWithUMIsMultipleOrientations","twopairsWithUMIsMultipleOrientations.sam").build();

    private final static Map<String, File> bams = new HashMap<>(sams.size());

    @BeforeTest
    public void prepareBams() throws IOException {

        sams.keySet().forEach(key -> {
            try {
                bams.put(key, convertSamToBam(sams.get(key)));
                bams.get(key).deleteOnExit();
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

        bamOutDir.deleteOnExit();
    }

    @AfterTest
    public void tearDown() {
        IOUtil.recursiveDelete(bamOutDir.toPath());
    }

    @DataProvider(name = "simpleTests")
    public Iterator<Object[]> simpleTestsData() {
        final List<Object[]> tests = new ArrayList<>(3);
        {
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nSites", 3);
            map.put("nDuplicateSets", 3);
            map.put("nMismatchingUMIsInCoOrientedBiDups",2);
            map.put("nMismatchingUMIsInContraOrientedBiDups",1);

            tests.add(new Object[]{"multipleContigs.vcf", "twopairsWithUMIsMultipleOrientations", map});
        }
        {
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nSites", 1);
            map.put("nDuplicateSets", 1);
            map.put("nDifferentAllelesTriDups", 0);
            map.put("nDifferentAllelesBiDups", 1);
            map.put("nReferenceAllelesBiDups", 0);
            map.put("nAlternateAllelesBiDups", 0);
            map.put("biSiteHeterogeneityRate", 1.0);
            map.put("biSiteHomogeneityRate", 0.0);
            map.put("nAlternateReads", 1);
            map.put("nReferenceReads", 1);

            tests.add(new Object[]{"hets.vcf", "twoPairs", map});
        }
        {// this tests the GQ cutoff
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nDifferentAllelesBiDups", 1);
            map.put("biSiteHeterogeneityRate", 1.0);
            map.put("biSiteHomogeneityRate", 0.0);
            map.put("nDifferentAllelesTriDups", 0);
            map.put("nSites", 1);

            tests.add(new Object[]{"twoSamplesHet.vcf", "twoPairs", map});
        }
        {
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nExactlyTriple", 1);
            map.put("nExactlyDouble", 1);
            map.put("nDuplicateSets", 2);
            map.put("nAlternateReads", 2);
            map.put("nReferenceReads", 3);

            tests.add(new Object[]{"hets.vcf", "aTriple", map});
        }
        {
            final Map<String, Object> map = new HashMap<>();

            map.put("nSites", 4);
            map.put("nTotalReads", 20);
            map.put("nDuplicateSets", 8);

            tests.add(new Object[]{"multipleContigs.vcf", "multipleContigs", map});
        }
        {
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nAlternateAllelesTriDups", 1);
            map.put("nMismatchingAllelesBiDups", 0); //we remove sites that have mismatching alleles, so this should be zero.

            tests.add(new Object[]{"hets_pos20.vcf", "aTriple", map});
        }
        tests.add(new Object[]{"hets_pos21_HOMREF_G.vcf", "aTriple", Collections.singletonMap("nReferenceAllelesTriDups", 1)});
        tests.add(new Object[]{"hets_pos20.vcf", "twoPairs", Collections.singletonMap("nAlternateAllelesBiDups", 1)});
        tests.add(new Object[]{"hets_pos21_HOMREF_G.vcf", "twoPairs", Collections.singletonMap("nReferenceAllelesBiDups", 1)});
        {
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nSites", 1);
            map.put("nThreeAllelesSites", 1);
            map.put("nAlternateAllelesTriDups", 0);
            map.put("nMismatchingAllelesBiDups", 0); //we remove sites that have mismatching alleles, so this should be zero.

            tests.add(new Object[]{"hets_pos22_IncorrectAlleles.vcf", "twoPairs", map});
        }
        //This tests the BQ cutoff
        tests.add(new Object[]{"hets_pos22_IncorrectAlleles.vcf", "aTriple", Collections.singletonMap("nMismatchingAllelesTriDups", 0)});

        {// tests for UMIs
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nSites", 1);
            map.put("nMismatchingUMIsInDiffBiDups", 1);
            map.put("nMatchingUMIsInDiffBiDups", 0);
            map.put("nMismatchingUMIsInSameBiDups", 0);
            map.put("nMatchingUMIsInSameBiDups", 0);
            map.put("nGoodBarcodes",1);
            map.put("nBadBarcodes",0);


            tests.add(new Object[]{"hets.vcf", "twoPairsWithUMIs", map});
        }
        {// tests for UMIs
            final Map<String, Object> map = new LinkedHashMap<>();

            map.put("nSites", 1);
            map.put("nMismatchingUMIsInDiffBiDups", 0);
            map.put("nMatchingUMIsInDiffBiDups", 0);
            map.put("nMismatchingUMIsInSameBiDups", 0);
            map.put("nMatchingUMIsInSameBiDups", 0);
            map.put("nGoodBarcodes",0);
            map.put("nBadBarcodes",1);

            tests.add(new Object[]{"hets.vcf", "twoPairsWithBadUMIs", map});
        }
        return tests.iterator();
    }

    @Test(dataProvider = "simpleTests")
    public void simpleTest(final String vcf, final String bam, final Map<String, Object> fieldValueMap) throws IOException, NoSuchFieldException, IllegalAccessException {

        final CollectIndependentReplicateMetrics est = new CollectIndependentReplicateMetrics();
        est.INPUT = bams.get(bam);
        est.VCF = new File(testdir, vcf);
        est.OUTPUT = IOUtil.newTempFile("singleHet", ".duplication_metric", new File[]{bamOutDir});
        est.MATRIX_OUTPUT = IOUtil.newTempFile("singleHet", ".duplication_matrix", new File[]{bamOutDir});
        est.SAMPLE = "SAMPLE1";

        est.OUTPUT.deleteOnExit();
        est.MATRIX_OUTPUT.deleteOnExit();

        est.doWork();

        final MetricsFile<IndependentReplicateMetric, Integer> retval = new MetricsFile<>();
        retval.read(new FileReader(est.OUTPUT));

        for (final Map.Entry<String, Object> fieldValue : fieldValueMap.entrySet()) {
            final String field = fieldValue.getKey();
            final Object expectedValue = fieldValue.getValue();
            final Field o = IndependentReplicateMetric.class.getField(field);
            assertEquals(o.get(retval.getMetrics().get(0)), expectedValue, field);
        }
    }

    /**
     * Converts a sam-file to a bam-file changing the extension from .sam to .bam
     *
     */
    private static File convertSamToBam(final String sam) throws IOException {
        final MergeSamFiles msf = new MergeSamFiles();
        final File bam = new File(bamOutDir, sam.replaceAll("sam$", "bam"));
        final int returnCode = msf.instanceMain(
                new String[]{
                        "INPUT=" + (new File(testdir, sam).getAbsolutePath()),
                        "CREATE_INDEX=true",
                        "OUTPUT=" + bam.getAbsolutePath()});
        Assert.assertEquals(returnCode, 0);

        return bam;
    }
}