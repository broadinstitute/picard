package picard.arrays;

public class ArraysUtils {

    static public ArraysCluster contaminateClusters(ArraysCluster main, ArraysCluster contaminant, final double contamination) {

        final double c = contamination;
        return new ArraysCluster(
                (1 - c) * main.meanX + c * contaminant.meanX,
                (1 - c) * main.meanY + c * contaminant.meanY,
                rmsLinearCompose(main.devX, 1-c, contaminant.devX, c),
                rmsLinearCompose(main.devY, 1-c, contaminant.devY, c));
    }

    static private double rmsLinearCompose(double lhs, double lhsScalar, double rhs, double rhsScalar) {
        return Math.sqrt(
                lhsScalar * lhs * lhs +
                rhsScalar * rhs * rhs);
    }
}
