/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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