package picard.vcf;

import java.io.File;

/**
 * Created by farjoun on 4/1/17.
 */
public class VcfUtils {

    static String UNCOMPRESSED_VCF_ENDING=".vcf";
    static String COMPRESSED_VCF_ENDING=".vcf.gz";
    static String BCF_ENDING=".bcf";

    /**
     * Checks if the suffix is one of those that are allowed for the various
     * formats that contain variants (currently vcf and bcf)
     */
    static public boolean isVariantFile(final File file){
        final String name = file.getName();

        if (name.endsWith(UNCOMPRESSED_VCF_ENDING)) return true;
        if (name.endsWith(COMPRESSED_VCF_ENDING)) return true;
        if (name.endsWith(BCF_ENDING)) return true;
        return false;
    }
}
