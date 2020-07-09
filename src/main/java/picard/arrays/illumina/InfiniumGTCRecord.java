package picard.arrays.illumina;

public class InfiniumGTCRecord {
    public final int rawXIntensity;
    public final int rawYIntensity;
    public final byte genotype;
    public final float genotypeScore;

    // Derived values:
    public final float normalizedXIntensity;
    public final float normalizedYIntensity;
    public final float RIlmn;
    public final float thetaIlmn;
    public final float bAlleleFreq;
    public final float logRRatio;

    InfiniumGTCRecord(final int rawXIntensity, final int rawYIntensity, final byte genotype, float genotypeScore,
                      final float normalizedXIntensity, final float normalizedYIntensity, final float RIlmn, final float thetaIlmn,
                      final float bAlleleFreq, final float logRRatio) {
        this.rawXIntensity = rawXIntensity;
        this.rawYIntensity = rawYIntensity;
        this.genotype = genotype;
        this.genotypeScore = genotypeScore;
        this.normalizedXIntensity = normalizedXIntensity;
        this.normalizedYIntensity = normalizedYIntensity;
        this.RIlmn = RIlmn;
        this.thetaIlmn = thetaIlmn;
        this.bAlleleFreq = bAlleleFreq;
        this.logRRatio = logRRatio;
    }

    public InfiniumGTCRecord(final String csvString) {
        String[] components = csvString.split(",");
        this.rawXIntensity = Integer.parseInt(components[0]);
        this.rawYIntensity = Integer.parseInt(components[1]);
        this.genotype = (byte)Integer.parseInt(components[2]);
        this.genotypeScore = Float.parseFloat(components[3]);
        this.normalizedXIntensity = Float.parseFloat(components[4]);
        this.normalizedYIntensity = Float.parseFloat(components[5]);
        this.RIlmn = Float.parseFloat(components[6]);
        this.thetaIlmn = Float.parseFloat(components[7]);
        this.bAlleleFreq = Float.parseFloat(components[8]);
        this.logRRatio = Float.parseFloat(components[9]);
    }
}
