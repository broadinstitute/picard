package picard.arrays;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import picard.PicardException;

public class ArraysUtils {

    static public ArraysCluster contaminateClusters(ArraysCluster main, ArraysCluster contaminant, final double contamination) {

        final double c = contamination;
        return new ArraysCluster(
                (1 - c) * main.meanX + c * contaminant.meanX,
                (1 - c) * main.meanY + c * contaminant.meanY,
                rmsLinearCompose(main.devX, 1 - c, contaminant.devX, c),
                rmsLinearCompose(main.devY, 1 - c, contaminant.devY, c));
    }

    static private double rmsLinearCompose(double lhs, double lhsScalar, double rhs, double rhsScalar) {
        return Math.sqrt(
                lhsScalar * lhs * lhs +
                        rhsScalar * rhs * rhs);
    }

    static public ArraysAssay getArraysAssay(VariantContext vc) {
        ArraysCluster AA = getCluster("AA", vc);
        ArraysCluster AB = getCluster("AB", vc);
        ArraysCluster BB = getCluster("BB", vc);

        String alleleA = vc.getAttributeAsString("ALLELE_A", "");
        String alleleB = vc.getAttributeAsString("ALLELE_B", "");

        if (alleleA.equals(vc.getReference().toString()) &&
                alleleB.equals(vc.getAlternateAllele(0).toString())) {
            return new ArraysAssay(AA, AB, BB);
        } else if (alleleB.equals(vc.getReference().toString()) &&
                alleleA.equals(vc.getAlternateAllele(0).toString())) {
            return new ArraysAssay(BB, AB, AA);
        } else {
            throw new PicardException("None of the alleles in the Assay matched the reference and the alternate allele");
        }
    }

    static private ArraysCluster getCluster(String clusterSuffix, final VariantContext vc) {

        if (!vc.hasAttribute("meanX_" + clusterSuffix))
            throw new PicardException("Missing attribute " + "meanX_" + clusterSuffix + " in " + vc);
        if (!vc.hasAttribute("meanY_" + clusterSuffix))
            throw new PicardException("Missing attribute " + "meanY_" + clusterSuffix + " in " + vc);
        if (!vc.hasAttribute("devX_" + clusterSuffix))
            throw new PicardException("Missing attribute " + "devX_" + clusterSuffix + " in " + vc);
        if (!vc.hasAttribute("devY_" + clusterSuffix))
            throw new PicardException("Missing attribute " + "devY_" + clusterSuffix + " in " + vc);

        return new ArraysCluster(
                vc.getAttributeAsDouble("meanX_" + clusterSuffix, 0),
                vc.getAttributeAsDouble("meanY_" + clusterSuffix, 0),
                vc.getAttributeAsDouble("devX_" + clusterSuffix, 0),
                vc.getAttributeAsDouble("devY_" + clusterSuffix, 0));
    }
}
