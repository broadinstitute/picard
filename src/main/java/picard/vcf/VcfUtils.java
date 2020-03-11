package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import java.io.File;

/**
 * Created by farjoun on 4/1/17.
 */
public class VcfUtils {

    /**
     * Checks if the suffix is one of those that are allowed for the various
     * formats that contain variants (currently vcf and bcf)
     */
    static public boolean isVariantFile(final File file){
        final String name = file.getName();

        return name.endsWith(IOUtil.VCF_FILE_EXTENSION) ||
                name.endsWith(IOUtil.COMPRESSED_VCF_FILE_EXTENSION) ||
                name.endsWith(IOUtil.BCF_FILE_EXTENSION);
    }
}
