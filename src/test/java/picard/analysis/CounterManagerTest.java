package picard.analysis;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

public class CounterManagerTest {

    private CounterManager testCounterManager;
    private CounterManager.Counter testCounter;
    private int arrayLength = 10;
    private int readLength = 10;
    private int OFFSET = 7;
    private CounterManager secondTestCounterManager;
    private CounterManager.Counter secondCounter;

    @BeforeTest
    public void setUp() {
        secondTestCounterManager = new CounterManager(arrayLength, readLength);
        secondTestCounterManager.setOffset(OFFSET);
        secondCounter = secondTestCounterManager.newCounter();
        testCounterManager = new CounterManager(arrayLength, readLength);
        testCounterManager.setOffset(OFFSET);
        testCounter = testCounterManager.newCounter();

        for (int i = 0; i < arrayLength; i++) {
            testCounter.increment(i + OFFSET);
            secondCounter.increment(i + OFFSET);
        }
        testCounter.increment(OFFSET);
        secondCounter.increment(OFFSET);
    }

    @Test
    public void testCounterInc() {
        Assert.assertEquals(2, testCounter.get(OFFSET), "Test method increment:");
    }

    @Test
    public void testForClearCounter() {
        testCounterManager.clear();
        Assert.assertEquals(0, testCounter.get(OFFSET), "The value of the array with index 0 must be 0 after clear manager:");
    }

    @Test
    public void testForCorrectOffsetAfterRebase() {
        secondTestCounterManager.checkOutOfBounds(11);
        Assert.assertEquals(11, secondTestCounterManager.getOffset(), "After rebase offset must be new int");
    }

    @Test
    public void testForCorrectCounterAfterRebase() {
        secondTestCounterManager.checkOutOfBounds(11);
        Assert.assertEquals(1, secondCounter.get(11), "The value of the array with index 0 must be 1 after rebase:");
    }

    @Test
    public void testForOutOfBoundCounter() {
        secondTestCounterManager.checkOutOfBounds(44);
        Assert.assertEquals(44, secondTestCounterManager.getOffset(), "New offset after clear must be 44:");
    }

    @Test
    public void testForCleanCounterAfter() {
        testCounterManager.checkOutOfBounds(88);
        Assert.assertEquals(0, testCounter.get(88), "The value of the array with index 0 must be 1 after clean:");
    }

    @Test
    public void testForCheckIncrement(){
        CounterManager testCounterManager = new CounterManager(arrayLength, readLength);
        testCounterManager.setOffset(0);
        CounterManager.Counter counter = testCounterManager.newCounter();
        for (int i=0; i<10; i++){
            counter.increment(1);
        }
        int[] testArray = new int[arrayLength];
        for (int i = 0; i< arrayLength; i++){
            testArray[i] = counter.get(i);
        }
        int[]templateArray = new int[arrayLength];
        templateArray[1] = 10;
        Assert.assertEquals(templateArray, testArray);
    }

    @Test(expectedExceptions = IndexOutOfBoundsException.class)
    public void testForWrongIndexInInc(){
        testCounter.increment(40);
    }

    @Test(expectedExceptions = IndexOutOfBoundsException.class)
    public void testForWrongIndexInGet(){
        testCounter.get(40);
    }
}
