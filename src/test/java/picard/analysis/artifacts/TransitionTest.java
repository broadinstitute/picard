package picard.analysis.artifacts;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.stream.Stream;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class TransitionTest {

    @DataProvider
    public Iterator<Object[]> allTransitions() {
        return Stream.of(Transition.values()).map(t -> new Object[]{t}).iterator();
    }

    @Test(dataProvider = "allTransitions")
    public void testIdentityAfterTwoComplement(final Transition transition) {
        Assert.assertEquals(transition.complement().complement(), transition);
    }

    @Test(dataProvider = "allTransitions")
    public void testTransitionOfSelf(final Transition transition) {
        final Transition ofSelf = Transition.transitionOf(transition.ref(), transition.call());
        Assert.assertEquals(ofSelf, transition);
        Assert.assertSame(ofSelf, transition);
    }

    @DataProvider
    public Object[][] badBases() {
        return new Object[][] {{Character.MIN_VALUE}, {Transition.Base.A.base - 1}, {'Z'}, {Character.MAX_VALUE}};
    }

    @Test(dataProvider = "badBases", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidRef(final char wrongBase) {
        Transition.transitionOf(wrongBase, 'A');
    }

    @Test(dataProvider = "badBases", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidCall(final char wrongBase) {
        Transition.transitionOf('A', wrongBase);
    }
}
