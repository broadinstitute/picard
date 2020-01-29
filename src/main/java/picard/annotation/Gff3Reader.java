package picard.annotation;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
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
        try (AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(gffFile.getAbsolutePath(), null, new Gff3Codec(), false)) {
            for (final Gff3Feature feature : reader.iterator()) {
                if (featureEvaluator.isASubClassOf(feature.getType(), "gene")) {
                    final Gene gene = new Gene(feature.getContig(), feature.getStart(), feature.getEnd(), feature.getStrand() == Strand.NEGATIVE, feature.getAttribute("ID"));
                    final Set<String> usedTranscriptNames = new HashSet<>();
                    for (final Gff3Feature child : feature.getChildren()) {
                        if (featureEvaluator.isASubClassOf(child.getType(), "CDS")) {
                            String transcriptName = child.getAttribute("ID");
                            if (usedTranscriptNames.contains(transcriptName)) {
                                final String transcriptNameBase = transcriptName;
                                int index = 1;
                                while (usedTranscriptNames.contains(transcriptName)) {
                                    transcriptName = transcriptNameBase + "_" + index;
                                    index++;
                                }
                            }
                            final Gene.Transcript tx = gene.addTranscript(transcriptName, child.getStart(), child.getEnd(), child.getStart(), child.getEnd(), 1);
                            tx.addExon(child.getStart(), child.getEnd());
                            usedTranscriptNames.add(transcriptName);
                        }
                        if (featureEvaluator.isASubClassOf(child.getType(), "transcript")) {
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

    private static void loadTranscriptIntoGene(final Gff3Feature transcript, final Gene gene) {
        if (!featureEvaluator.isASubClassOf(transcript.getType(), "transcript")) {
            throw new PicardException("tried to load non-transcript type " + transcript.getType() + " as transcript");
        }
        final List<Gff3Feature> cdsChildren = transcript.getChildren().stream().filter(f -> featureEvaluator.isASubClassOf(f.getType(), "CDS")).collect(Collectors.toList());
        final int cdsStart = cdsChildren.size()>0? cdsChildren.stream().map(Gff3Feature::getStart).min(Comparator.comparing(Integer::valueOf)).get() : 0;
        final int cdsEnd = cdsChildren.size()>0? cdsChildren.stream().map(Gff3Feature::getEnd).max(Comparator.comparing(Integer::valueOf)).get() : 0;

        final List<Gff3Feature> exons = transcript.getChildren().stream().filter(f -> featureEvaluator.isASubClassOf(f.getType(), "exon")).collect(Collectors.toList());

        boolean notYetLoaded = true;
        String transcriptName = transcript.getAttribute("ID");
        final String transcriptBaseName = transcriptName;
        int index = 1;
        while (notYetLoaded) {
            try {
                Gene.Transcript tx = gene.addTranscript(transcriptName, transcript.getStart(), transcript.getEnd(), cdsStart, cdsEnd, exons.size() > 0 ? exons.size() : 1);
                if (exons.size() > 0) {
                    for (final Gff3Feature exon : exons) {
                        tx.addExon(exon.getStart(), exon.getEnd());
                    }
                } else {
                    tx.addExon(cdsStart, cdsEnd);
                }
                notYetLoaded = false;
            } catch (final AnnotationException ex) {
                transcriptName = transcriptBaseName + "_"+index;
                index++;
            }
        }

    }

    public static OverlapDetector<Interval> loadRibosomalIntervals(final File gffFile) {
        final OverlapDetector<Interval> overlapDetector = new OverlapDetector<>(0, 0);
        try (AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(gffFile.getAbsolutePath(), null, new Gff3Codec(), false)) {
            final String rRNAIRI = "rRNA";
            for (final Gff3Feature topFeature : reader.iterator()) {
                for(final Gff3Feature feature : topFeature.flatten()) {
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
