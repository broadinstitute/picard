package picard.arrays;

import picard.PicardException;

public class ArraysAssay {

    final private ArraysCluster AA;
    final private ArraysCluster AB;
    final private ArraysCluster BB;

    public ArraysAssay(final ArraysCluster AA, final ArraysCluster AB, final ArraysCluster BB) {
        this.AA = AA;
        this.AB = AB;
        this.BB = BB;
    }

    public ArraysCluster getAA() {
        return AA;
    }

    public ArraysCluster getAB() {
        return AB;
    }

    public ArraysCluster getBB() {
        return BB;
    }
    // 0-AA 1-AB 2-BB
    public ArraysCluster getCluster(int genotype){
        switch(genotype){
            case 0:
                return getAA();
            case 1:
                return getAB();
            case 2:
                return getBB();
            default:
                throw new PicardException("genotype must be 0,1, or 2. Got " + genotype);
        }
    }
}
