package picard.arrays;


import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;

class ZCallPedFile {
    private final Log log = Log.getInstance(ZCallPedFile.class);
    private final ProgressLogger logger = new ProgressLogger(log, 10000);

    private static final int OFFSET = 6;
    private final Map<String, String> snpToAlleleMap = new HashMap<>();

    private void addAllele(String snp, String allele) {
        logger.record("0", 0);
        snpToAlleleMap.put(snp, allele);
    }

    String getAlleles(String snp) {
        return snpToAlleleMap.get(snp);
    }

    public static ZCallPedFile fromFile(File pedFile, File mapFile) throws FileNotFoundException {
        String[] pedFileFields = IOUtil.slurp(pedFile).split(" ");
        String[] mapFileLines = IOUtil.slurpLines(mapFile).toArray(new String[0]);

        ZCallPedFile zCallPedFile = new ZCallPedFile();

        /* first six fields are ignored
            Family ID
            Individual ID
            Paternal ID
            Maternal ID
            Sex (1=male; 2=female; other=unknown)
            Phenotype
         */
        //two fields for each snp (each allele)
        for (int i = 0; i < mapFileLines.length; i++) {
            int index = (i * 2) + OFFSET;
            String alleles = pedFileFields[index] + pedFileFields[index + 1];
            zCallPedFile.addAllele(mapFileLines[i].split("\t")[1], alleles);
        }
        return zCallPedFile;
    }
}
