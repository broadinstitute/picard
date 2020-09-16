package picard.arrays;

import picard.PicardException;

public class ArraysAssay {

    final private ArraysCluster HomRef;
    final private ArraysCluster Het;
    final private ArraysCluster HomVar;

    public ArraysAssay(final ArraysCluster AA, final ArraysCluster AB, final ArraysCluster BB) {
        this.HomRef = AA;
        this.Het = AB;
        this.HomVar = BB;
    }

    public ArraysCluster getHomRef() {
        return HomRef;
    }

    public ArraysCluster getHet() {
        return Het;
    }

    public ArraysCluster getHomVar() {
        return HomVar;
    }
    // 0-HOM_REF 1-HET 2-HOM_VAR
    public ArraysCluster getCluster(int genotype){
        switch(genotype){
            case 0:
                return getHomRef();
            case 1:
                return getHet();
            case 2:
                return getHomVar();
            default:
                throw new PicardException("genotype must be 0,1, or 2. Got " + genotype);
        }
    }
}
