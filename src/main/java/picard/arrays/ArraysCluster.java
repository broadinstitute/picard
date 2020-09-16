package picard.arrays;


import org.apache.commons.math3.analysis.function.Gaussian;

public class ArraysCluster {
    final public float meanX;
    final public float meanY;

    final public float devX;
    final public float devY;

//    final Gaussian xDistribution;
//    final Gaussian yDistribution;


    public ArraysCluster(float meanX, float meanY, float devX, float devY) {
        this.meanX = meanX;
        this.meanY = meanY;
        this.devX = devX;
        this.devY = devY;


    }

//    public float likelihood(final float x, final float y) {
//        return Gaussian()

//    }
}
