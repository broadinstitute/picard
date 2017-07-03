package picard.analysis.replicates;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.*;
import htsjdk.tribble.annotation.Strand;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;

public class ConverterTest {


    @Test
    public void refFlatToGenePredConvertionTest() {

        RefFlatRecord refFlatRecord = new RefFlatRecord("myGene", "myGene", "sequence", Strand.POSITIVE, 49, 500, 74, 400, 2, new int[]{49, 249}, new int[]{200, 500});
        Stream<RefFlatRecord> refFlatRecordStream = Stream.of(refFlatRecord);

        GenePredRecord genePredRecord = new GenePredRecord("myGene", "sequence", Strand.POSITIVE, 49, 500, 74, 400, 2, new int[]{49, 249}, new int[]{200, 500});

        List<GenePredRecord> genePredRecordList = new ArrayList<>();
        genePredRecordList.add(genePredRecord);

        Assert.assertEquals(Converter.refFlatToGenePred(refFlatRecordStream).collect(toList()), genePredRecordList);
    }

    @Test
    public void gtfToGenePredConvertionTest() {
        String source = "mySource";
        Map<String, String> firstTranscriptAttributes = new HashMap<>();
        firstTranscriptAttributes.put("gene_id", "myGene1");
        firstTranscriptAttributes.put("transcript_id", "myTranscript1");
        Strand strand = Strand.POSITIVE;
        String sequenceName = "sName";

        GtfRecord firstGtf = new GtfRecord(sequenceName, source, Feature.EXON, 60, 210, ".", strand, ".", firstTranscriptAttributes);
        GtfRecord secondGtf = new GtfRecord(sequenceName, source, Feature.CDS, 85, 210, ".", strand, ".", firstTranscriptAttributes);

        Map<String, String> secondTranscriptAttributes = new HashMap<>();
        secondTranscriptAttributes.put("gene_id", "myGene1");
        secondTranscriptAttributes.put("transcript_id", "myTranscript2");

        GtfRecord thirdGtf = new GtfRecord(sequenceName, source, Feature.EXON, 260, 510, ".", strand, ".", secondTranscriptAttributes);
        GtfRecord fourthGtf = new GtfRecord(sequenceName, source, Feature.CDS, 260, 410, ".", strand, ".", secondTranscriptAttributes);

        List<GtfRecord> gtfRecordList = new ArrayList<>();
        gtfRecordList.add(firstGtf);
        gtfRecordList.add(secondGtf);
        gtfRecordList.add(thirdGtf);
        gtfRecordList.add(fourthGtf);

        Stream<GtfRecord> gtfRecordsStream = gtfRecordList.stream();

        int[] firstExonStarts = new int[]{59};
        int[] firstExonEnds = new int[]{210};
        GenePredRecord firstGenePred = new GenePredRecord("myTranscript1","sName", Strand.POSITIVE, 59, 210,84, 210,1, firstExonStarts, firstExonEnds);

        int[] secondExonStarts = new int[]{259};
        int[] secondExonEnds = new int[]{510};
        GenePredRecord secondGenePred = new GenePredRecord("myTranscript2","sName", Strand.POSITIVE, 259, 510,259, 410,1, secondExonStarts, secondExonEnds);

        List<GenePredRecord> genePredRecordList = new ArrayList<>();
        genePredRecordList.add(firstGenePred);
        genePredRecordList.add(secondGenePred);

        Assert.assertEquals(Converter.gtfToGenePred(gtfRecordsStream).collect(toList()), genePredRecordList);
    }
}
