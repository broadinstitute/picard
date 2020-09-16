package picard.arrays;

import org.apache.commons.math3.analysis.function.Gaussian;

public class ArraysCluster {

    final double meanX;
    final double meanY;
    final double devX;
    final double devY;

    final Gaussian xDistribution;
    final Gaussian yDistribution;

    public ArraysCluster(double meanX, double meanY, double devX, double devY) {
        this.meanX = meanX;
        this.meanY = meanY;
        this.devX = devX;
        this.devY = devY;

        xDistribution = new Gaussian(meanX, devX);
        yDistribution = new Gaussian(meanY, devY);
    }

    public double likelihood(final double x, final double y) {
        return xDistribution.value(x)*yDistribution.value(y);
    }
}
