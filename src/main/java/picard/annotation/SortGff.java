package picard.annotation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3Writer;
import htsjdk.tribble.gff.SequenceRegion;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * <h3> Summary </h3>
 * <p> This tool sorts a gff3 file by coordinates, so that it can be indexed.
 * It additionally adds flush directives where possible, which can significantly reduce the memory footprint of downstream tools.
 * Sorting of multiple contigs can be specified by a sequence dictionary; if no sequence dictionary is specified, contigs are sorted lexicographically. </p>
 *
 * <h3> Usage Examples </h3>
 * <h4> 1. Sort gff3 file, add flush directives.  Contigs will be sorted lexicographically.</h4>
 * <pre>
 * java -jar picard.jar SortGff
 *      I=input.gff3
 *      O=output.gff3
 * </pre>
 *
 * <h4> 2. Sort gff3 file, add flush directives.  Contigs will be sorted according to order in sequence dictionary</h4>
 * <pre>
 * java -jar picard.jar SortGff
 *      I=input.gff3
 *      O=output.gff3
 *      SD=dictionary.dict
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = SortGff.USAGE_DETAILS,
        oneLineSummary = SortGff.USAGE_SUMMARY,
        programGroup = OtherProgramGroup.class)
public class SortGff extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Sorts a gff3 file, and adds flush directives";
    static final String USAGE_DETAILS = "<h3> Summary </h3>\n" +
            "  <p> This tool sorts a gff3 file by coordinates, so that it can be indexed.\n" +
            " It additionally adds flush directives where possible, which can significantly reduce the memory footprint of downstream tools.\n" +
            " Sorting of multiple contigs can be specified by a sequence dictionary; if no sequence dictionary is specified, contigs are sorted lexicographically. </p>\n" +
            "\n" +
            " <h3> Usage Examples </h3>\n" +
            " <h4> 1. Sort gff3 file, add flush directives.  Contigs will be sorted lexicographically.</h4>\n" +
            " <pre>\n" +
            " java -jar picard.jar SortGff\n" +
            "      I=input.gff3\n" +
            "      O=output.gff3\n" +
            " </pre>\n" +
            "\n" +
            " <h4> 2. Sort gff3 file, add flush directives.  Contigs will be sorted according to order of sequence dictionary</h4>\n" +
            " <pre>\n" +
            " java -jar picard.jar SortGff\n" +
            "      I=input.gff3\n" +
            "      O=output.gff3\n" +
            "      SD=dictionary.dict\n" +
            " </pre>";

    @Argument(doc = "Input Gff3 file to sort.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Sorted Gff3 output file.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "Dictionary to sort contigs by.  If dictionary is not provided, contigs are sorted lexicographically.", shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, optional = true)
    public File SEQUENCE_DICTIONARY;

    private final Log log = Log.getInstance(SortGff.class);

    private final Map<String, Integer> latestStartMap = new HashMap<>();
    private int latestStart = 0;
    String latestChrom;

    private static class FeatureComparator implements Comparator<Gff3Feature> {
        final SAMSequenceDictionary dict;

        FeatureComparator(final SAMSequenceDictionary dict) {
            this.dict = dict;
        }

        public int compare(final Gff3Feature f1, final Gff3Feature f2) {
            int comp = dict == null ? f1.getContig().compareTo(f2.getContig()) : dict.getSequenceIndex(f1.getContig()) - dict.getSequenceIndex(f2.getContig());
            if (comp == 0) {
                comp = f1.getStart() - f2.getStart();
            }

            return comp;
        }
    }

    @Override
    protected int doWork() {
        final Gff3Codec inputCodec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
        if (!(inputCodec.canDecode(INPUT.toString()))) {
            throw new IllegalArgumentException("Input file " + INPUT + " cannot be read by Gff3Codec");
        }

        final SAMSequenceDictionary dict = SEQUENCE_DICTIONARY == null? null : SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY.toPath());
        SortingCollection<Gff3Feature> sorter = SortingCollection.newInstance(
                Gff3Feature.class, new Gff3SortingCollectionCodec(), new FeatureComparator(dict), 500000
        );

        final ProgressLogger progressRead = new ProgressLogger(log, (int) 1e4, "Read");



        try (final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(INPUT.getAbsolutePath(), null, inputCodec, false)) {
            for (final Gff3Feature shallowFeature : reader.iterator()) {
                if (shallowFeature != null) {
                    sorter.add(shallowFeature);
                    progressRead.record(shallowFeature.getContig(), shallowFeature.getStart());

                    final String featureID = shallowFeature.getID();
                    final List<String> parentIDs = shallowFeature.getAttribute("Parent");

                    if (featureID != null) {
                        //update latestStartMap for this feature based on its position
                        latestStartMap.compute(featureID, (k, v) -> v == null ? shallowFeature.getStart() : Math.max(shallowFeature.getStart(), v));
                    }
                    final Integer currentLatestStart = featureID == null ? shallowFeature.getStart() : latestStartMap.get(featureID);

                    //update latestStartMap for each parent based on this feature's latest start
                    for (final String parentID : parentIDs) {
                        latestStartMap.compute(parentID, (k, v) -> v == null ? currentLatestStart : Math.max(currentLatestStart, v));
                    }
                }
            }

        } catch(final IOException e) {
            throw new PicardException("Error reading file " + INPUT, e);
        }

        final ProgressLogger progressWrite = new ProgressLogger(log, (int) 1e4, "Wrote");
        try(final Gff3Writer writer = new Gff3Writer(OUTPUT.toPath())) {
            //add comments and sequence regions
            for (final String comment : inputCodec.getComments()) {
                writer.addComment(comment);
            }

            for (final SequenceRegion sequenceRegion : inputCodec.getSequenceRegions()) {
                writer.addDirective(Gff3Codec.Gff3Directive.SEQUENCE_REGION_DIRECTIVE, sequenceRegion);
            }

            //add features
            final Iterator<Gff3Feature> sortedIterator = sorter.iterator();
            if (sortedIterator.hasNext()) {
                final Gff3Feature firstFeature = sortedIterator.next();
                writer.addFeature(firstFeature);
                updateLatestStart(firstFeature);
                while (sortedIterator.hasNext()) {
                    final Gff3Feature feature = sortedIterator.next();
                    //check if we can add a flush directive before this feature
                    if (feature.getStart() > latestStart || !feature.getContig().equals(latestChrom)) {
                        writer.addDirective(Gff3Codec.Gff3Directive.FLUSH_DIRECTIVE);
                    }
                    writer.addFeature(feature);
                    updateLatestStart(feature);
                    progressWrite.record(feature.getContig(), feature.getStart());
                }
            }
        } catch (final IOException ex) {
            throw new PicardException("Error opening  " + OUTPUT + " to write to", ex);
        }

        return 0;
    }

    private void updateLatestStart(final Gff3Feature feature) {
        final String featureID = feature.getID();
        if (featureID != null) {
            latestStart = Math.max(latestStart, latestStartMap.get(feature.getID()));
        }
        //check for parents which may have been found after this feature
        final List<String> parentIDs = feature.getAttribute("Parent");
        latestStart = Math.max(latestStart, parentIDs.stream().map(latestStartMap::get).max(Integer::compareTo).orElse(feature.getStart()));

        latestChrom = feature.getContig();
    }

    //a sorting collection to sort gff3 features.  Iterator returns SHALLOW decoded features; these features are not linked to any parents, children, or co-features, so any linking required must be handled elsewhere
    class Gff3SortingCollectionCodec implements SortingCollection.Codec<Gff3Feature> {
        LineIterator lineIteratorIn;
        Gff3Codec gff3Codec;
        Gff3Writer gff3Writer;


        Gff3SortingCollectionCodec() {
            gff3Codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
        }

        @Override
        public SortingCollection.Codec<Gff3Feature> clone() {
            return new Gff3SortingCollectionCodec();
        }

        @Override
        public void setOutputStream(final OutputStream os) {
            gff3Writer = new Gff3Writer(new PrintStream(os));
        }

        @Override
        public void setInputStream(final InputStream is) {
            lineIteratorIn = gff3Codec.makeSourceFromStream(is);
        }

        @Override
        public Gff3Feature decode() {
            while (!gff3Codec.isDone(lineIteratorIn)) {
                try {
                    final Gff3Feature feature = gff3Codec.decode(lineIteratorIn);

                    if (feature == null) {
                        continue;
                    }
                    return feature;
                } catch (final IOException ex) {
                    throw new PicardException("Error decoding feature in Gff3SortingCollectionCodec", ex);
                }
            }
            return null;
        }

        @Override
        public void encode(final Gff3Feature val) {
            gff3Writer.addFeature(val);
        }
    }
}
