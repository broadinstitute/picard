package picard.annotation;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;

import java.util.Collection;

import org.testng.Assert;
import org.testng.annotations.Test;

public class GeneTest {

	
	// Genes with the same chromosome/start/end had the same hashcode/equals, so were overwritten when added to the overlap detector.
	@Test
	public void hashCodeEqualsTest () {
		OverlapDetector<String> od = new OverlapDetector<String>(0, 0);
		Interval i1= new Interval("1", 1, 10);
		Interval i2= new Interval("1", 1, 10);
		Interval testInterval= new Interval("1", 1, 10);
		
		od.addLhs("Foo", i1);
		od.addLhs("Bar", i2);
		Collection<String> result = od.getOverlaps(testInterval);
		Assert.assertEquals(result.size(), 2);
		
		OverlapDetector<Gene> od2 = new OverlapDetector<Gene>(0, 0);
		
		Gene g1 = new Gene (i1.getSequence(), i1.getStart(), i1.getEnd(), true, "Foo");
		Gene g2 = new Gene (i2.getSequence(), i2.getStart(), i2.getEnd(), true, "Bar");
		
		od2.addLhs(g1, i1);
		od2.addLhs(g2, i2);
		
		Collection<Gene> result2 = od2.getOverlaps(testInterval);
		Assert.assertEquals(result2.size(), 2);
	}
}
