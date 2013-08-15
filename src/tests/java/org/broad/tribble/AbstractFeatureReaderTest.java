package org.broad.tribble;

import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.readers.LineIterator;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.testng.annotations.Test;

import java.io.IOException;

import static org.testng.Assert.*;

/**
 * @author jacob
 * @date 2013-Apr-10
 */
public class AbstractFeatureReaderTest {

    final static String HTTP_INDEXED_VCF_PATH = "http://www.broadinstitute.org/~picard/testdata/ex2.vcf";
    final static String LOCAL_MIRROR_HTTP_INDEXED_VCF_PATH = VariantBaseTest.variantTestDataRoot + "ex2.vcf";

    /**
     * Asserts readability and correctness of VCF over HTTP.  The VCF is indexed and requires and index.
     */
    @Test
    public void testVcfOverHTTP() throws IOException {
        final VCFCodec codec = new VCFCodec();
        final AbstractFeatureReader<VariantContext, LineIterator> featureReaderHttp =
                AbstractFeatureReader.getFeatureReader(HTTP_INDEXED_VCF_PATH, codec, true); // Require an index to
        final AbstractFeatureReader<VariantContext, LineIterator> featureReaderLocal =
                AbstractFeatureReader.getFeatureReader(LOCAL_MIRROR_HTTP_INDEXED_VCF_PATH, codec, false);
        final CloseableTribbleIterator<VariantContext> localIterator = featureReaderLocal.iterator();
        for (final Feature feat : featureReaderHttp.iterator()) {
            assertEquals(feat.toString(), localIterator.next().toString());
        }
        assertFalse(localIterator.hasNext());
    }

    @Test
    public void testLoadBEDFTP() throws Exception {
        final String path = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands with spaces.hg18.bed";
        final BEDCodec codec = new BEDCodec();
        final AbstractFeatureReader<BEDFeature, LineIterator> bfs = AbstractFeatureReader.getFeatureReader(path, codec, false);
        for (final Feature feat : bfs.iterator()) {
            assertNotNull(feat);
        }
    }
}
