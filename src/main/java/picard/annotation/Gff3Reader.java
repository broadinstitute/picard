package picard.annotation;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Loads gene annotations from a gff3 file into an OverlapDetector<Gene>.  Discards annotations that are not
 * internally consistent, e.g. transcripts on different chromosomes or different strands.
 */

public class Gff3Reader {
    final static Gff3FeatureEvaluator featureEvaluator = new Gff3FeatureEvaluator();

    public static OverlapDetector<Gene> load(final File gffFile) {
        final OverlapDetector<Gene> overlapDetector = new OverlapDetector<Gene>(0, 0);
        try (AbstractFeatureReader<GtfFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(gffFile.getAbsolutePath(), null, new Gff3Codec(), false)) {
            for (final GtfFeature feature : reader.iterator()) {
                if(featureEvaluator.isASubClassOf(feature.getType(), "http://purl.obolibrary.org/obo/SO_0001263")) {
                    continue;
                }
                if (featureEvaluator.isASubClassOf(feature.getType(), Gff3FeatureEvaluator.DEFAULT_GENE_IRI)) {
                    final Gene gene = new Gene(feature.getContig(), feature.getStart(), feature.getEnd(), feature.getStrand() == Strand.NEGATIVE, feature.getAttribute("ID"));
                    for (final GtfFeature child : feature.getChildren()) {
                        if (featureEvaluator.isASubClassOf(child.getType(), Gff3FeatureEvaluator.DEFAULT_CDS_IRI)) {
                            final Gene.Transcript tx = gene.addTranscript(child.getAttribute("ID"), child.getStart(), child.getEnd(), child.getStart(), child.getEnd(), 1);
                            tx.addExon(child.getStart(), child.getEnd());
                        }
                        if (featureEvaluator.isASubClassOf(child.getType(), Gff3FeatureEvaluator.DEFAULT_TRANSCRIPT_IRI)) {
                            loadTranscriptIntoGene(child, gene);
                        }
                    }
                    overlapDetector.addLhs(gene, gene);
                }
            }
        } catch (IOException e) {
            throw new PicardException("Error reading GFF3 file " + gffFile.getAbsolutePath());
        }
        return overlapDetector;
    }

    private static void loadTranscriptIntoGene(final GtfFeature transcript, final Gene gene) {
        if (!featureEvaluator.isASubClassOf(transcript.getType(), Gff3FeatureEvaluator.DEFAULT_TRANSCRIPT_IRI)) {
            throw new PicardException("tried to load non-transcript type " + transcript.getType() + " as transcript");
        }
        final List<GtfFeature> cdsChildren = transcript.getChildren().stream().filter(f -> featureEvaluator.isASubClassOf(f.getType(), Gff3FeatureEvaluator.DEFAULT_CDS_IRI)).collect(Collectors.toList());
        final int cdsStart = cdsChildren.size()>0? cdsChildren.stream().map(GtfFeature::getStart).min(Comparator.comparing(Integer::valueOf)).get() : 0;
        final int cdsEnd = cdsChildren.size()>0? cdsChildren.stream().map(GtfFeature::getEnd).max(Comparator.comparing(Integer::valueOf)).get() : 0;

        final List<GtfFeature> exons = transcript.getChildren().stream().filter(f -> featureEvaluator.isASubClassOf(f.getType(), Gff3FeatureEvaluator.DEFAULT_EXON_IRI)).collect(Collectors.toList());

        Gene.Transcript tx = gene.addTranscript(transcript.getAttribute("ID"), transcript.getStart(), transcript.getEnd(), cdsStart, cdsEnd, exons.size() > 0 ? exons.size() : 1);
        if (exons.size() > 0) {
            for (final GtfFeature exon : exons) {
                tx.addExon(exon.getStart(), exon.getEnd());
            }
        } else {
            tx.addExon(cdsStart, cdsEnd);
        }

    }

    public static OverlapDetector<Interval> loadRibosomalIntervals(final File gffFile) {
        final OverlapDetector<Interval> overlapDetector = new OverlapDetector<>(0, 0);
        try (AbstractFeatureReader<GtfFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(gffFile.getAbsolutePath(), null, new Gff3Codec(), false)) {
            final String rRNAIRI = "http://purl.obolibrary.org/obo/SO_0000252";
            for (final GtfFeature topFeature : reader.iterator()) {
                for(final GtfFeature feature : topFeature.flatten()) {
                    if (featureEvaluator.isASubClassOf(feature.getType(), rRNAIRI)) {
                        final Interval interval = new Interval(feature.getContig(), feature.getStart(), feature.getEnd());
                        overlapDetector.addLhs(interval, interval);
                    }
                }
            }
        } catch (IOException e) {
            throw new PicardException("Error reading GFF3 file " + gffFile.getAbsolutePath());
        }
        return overlapDetector;
    }
}
