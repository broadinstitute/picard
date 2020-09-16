package picard.arrays;

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
}
