package picard.arrays.illumina;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class CreateExtendedIlluminaManifestTest {
    private static final Path TEST_DATA_DIR = Paths.get("testdata/picard/arrays/illumina/");

    @Test
    public void testFlagDuplicates() throws IOException {
        final File extendedManifestFile = TEST_DATA_DIR.resolve("GDA-8v1-0_A5.2.0.extended.csv").toFile();
        Build37ExtendedIlluminaManifest extendedManifest = new Build37ExtendedIlluminaManifest(extendedManifestFile);
        final List<Build37ExtendedIlluminaManifestRecord> records = new ArrayList<>();

        final Iterator<Build37ExtendedIlluminaManifestRecord> iterator = extendedManifest.extendedIterator();
        while (iterator.hasNext()) {
            records.add(iterator.next());
        }
        final Map<String, Float> nameToGentrainScore = new HashMap<>();
        // These two are simple SNPs in the manifest at the same chr/pos
        nameToGentrainScore.put("1:100343253-CT", 0.761f);      // Note - this is the 2nd entry in the manifest
        nameToGentrainScore.put("rs199660743", 0.887f);

        // These three exist in the manifest at the same chr/pos
        // This one is a SNP (so different alleles than the next two) and so is not flagged as a dupe
        nameToGentrainScore.put("rs118203419", 0.808f);
        // These two are indels, and although on different strands, have the same alleles.
        // So one will be flagged as a duplicate.  That will be the one with the lower GenTrain score
        nameToGentrainScore.put("9:132921819_IlmnFwd", 0.823f);     // Note - this is the 5th entry in the manifest
        nameToGentrainScore.put("9:132921819_IlmnRev", 0.848f);

        final CreateExtendedIlluminaManifest clp = new CreateExtendedIlluminaManifest();
        Set<Integer> duplicateIndices = clp.flagDuplicates(records, nameToGentrainScore);
        Assert.assertEquals(duplicateIndices.size(), 2);
        Assert.assertTrue(duplicateIndices.contains(1));        // Finds the 2nd entry (0-based)
        Assert.assertTrue(duplicateIndices.contains(4));        // Finds the 5th entry (0-based)
        int index = 0;
        for (Build37ExtendedIlluminaManifestRecord record : records) {
            if (!record.isFail() && record.isDupe()) {
                Assert.assertTrue(duplicateIndices.contains(index));
            } else {
                Assert.assertFalse(duplicateIndices.contains(index));
            }
            index++;
        }
    }

    // Test that we retrieve the rsIds properly from dbSnp VCF for a SNP and an indel
    @Test
    public void testGenerateLocusToRsidMap() throws IOException {
        final File dbSnpVcf = TEST_DATA_DIR.resolve("dbsnp_138.b37.fortesting.vcf").toFile();
        File indexedDbSnpVcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(dbSnpVcf, "createExtendedIlluminaManifestTest.tmp.");
        final SAMSequenceDictionary sequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(indexedDbSnpVcf.toPath());
        IntervalList manifestSnpIntervals = new IntervalList(sequenceDictionary);
        IntervalList manifestIndelIntervals = new IntervalList(sequenceDictionary);

        final VCFFileReader vcfFileReader = new VCFFileReader(indexedDbSnpVcf, false);

        int intervals = 0;
        for(final VariantContext vc : vcfFileReader){
            String name = vc.getID();
            final int intervalEnd=vc.getCommonInfo().getAttributeAsInt("END",vc.getEnd());
            if (".".equals(name) || name == null)
                name = "interval-" + (++intervals);
            Interval interval = new Interval(vc.getContig(), vc.getStart(), intervalEnd, false, name);
            if (vc.isSNP()) {
                manifestSnpIntervals.add(interval);
            } else {
                manifestIndelIntervals.add(interval);
            }
        }

        vcfFileReader.close();

        final CreateExtendedIlluminaManifest clp = new CreateExtendedIlluminaManifest();
        Map<String, String> snpLocusToRsId = clp.generateLocusToRsidMap(indexedDbSnpVcf, manifestSnpIntervals);
        Assert.assertTrue(snpLocusToRsId.containsKey("1.10109"));
        Assert.assertEquals("rs376007522", snpLocusToRsId.get("1.10109"));
        Map<String, String> indelLocusToRsId = clp.generateLocusToRsidMap(indexedDbSnpVcf, manifestIndelIntervals);
        Assert.assertTrue(indelLocusToRsId.containsKey("1.10110"));     // indel spans two bases.
        Assert.assertEquals("rsDummyIndel", indelLocusToRsId.get("1.10110"));
        Assert.assertEquals("rsDummyIndel", indelLocusToRsId.get("1.10111"));
    }
}
