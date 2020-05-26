package picard.annotation;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3Writer;
import htsjdk.tribble.gff.SequenceRegion;
import htsjdk.tribble.readers.LineIterator;
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

@CommandLineProgramProperties(
        summary = SortGff.USAGE_DETAILS,
        oneLineSummary = SortGff.USAGE_SUMMARY,
        programGroup = OtherProgramGroup.class)
public class SortGff extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Sorts a Gff3 file, and adds flush directives";
    static final String USAGE_DETAILS = "Sorts a Gff3 file, and adds flush directives";

    @Argument(doc = "Input Gff3 file to sort.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Sorted Gff3 output file.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    private final Log log = Log.getInstance(SortGff.class);

    private final Map<String, Integer> latestStartMap = new HashMap<>();
    private final Map<String, Set<String>> parentChildrenMap = new HashMap<>();

    @Override
    protected int doWork() {
        final Gff3Codec codec = new Gff3Codec();
        if (!(codec.canDecode(INPUT.toString()))) {
            throw new IllegalArgumentException("Input file " + INPUT + " cannot be read by Gff3Codec");
        }

        SortingCollection<Gff3Feature> sorter = SortingCollection.newInstance(
                Gff3Feature.class,
                new Gff3SortingCollectionCodec(),
                (Gff3Feature f1, Gff3Feature f2) -> f1.getContig().compareTo(f2.getContig()) == 0 ? f1.getStart() - f2.getStart() : f1.getContig().compareTo(f2.getContig()),
                500000
        );

        final ProgressLogger progressRead = new ProgressLogger(log, (int) 1e4, "Read");

        final Gff3Codec inputCodec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);

        try (final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(INPUT.getAbsolutePath(), null, inputCodec, false)) {
            for (final Gff3Feature feature : reader.iterator()) {
                if (feature != null) {
                    sorter.add(feature);
                    progressRead.record(feature.getContig(), feature.getStart());

                    final String featureID = feature.getID();
                    final List<String> parentIDs = feature.getAttribute("Parent");

                    if (featureID != null) {
                        //update latestStartMap for this feature based on its position
                        latestStartMap.compute(featureID, (k, v) -> v == null ? feature.getStart() : Math.max(feature.getStart(), v));

                        //update latestStartMap for this feature based on its parents' latest starts
                        final Integer parentsLatestStart = parentIDs.stream().map(latestStartMap::get).max(Integer::compareTo).orElse(0);
                        final Integer currentLatestStart = latestStartMap.compute(featureID, (k, v) -> Math.max(v, parentsLatestStart));

                        //check for any child features of this feature which have been seen before it
                        if (parentChildrenMap.containsKey(featureID)) {
                            for (final String childID : parentChildrenMap.get(featureID)) {
                                //update latestStartMap for child
                                latestStartMap.compute(featureID, (k, v) -> Math.max(v, currentLatestStart));
                            }
                        }
                    }
                    final Integer currentLatestStart = featureID == null? feature.getStart() : latestStartMap.get(featureID);

                    //update latestStartMap for each parent based on this feature's latest start, and add this feature as child to each of its parents
                    for (final String parentID : parentIDs) {
                        latestStartMap.compute(parentID, (k, v) -> v == null? currentLatestStart : Math.max(currentLatestStart, v));

                        if (featureID != null) {
                            final Set<String> childrenOfParent = parentChildrenMap.computeIfAbsent(parentID, id -> new HashSet<>());
                            childrenOfParent.add(featureID);
                        }
                    }
                }
            }

        } catch(final IOException e) {
            throw new PicardException("Error reading file " + INPUT, e);
        }

        final ProgressLogger progressWrite = new ProgressLogger(log, (int) 1e4, "Wrote");
        try(final Gff3Writer writer = new Gff3Writer(OUTPUT.toPath())) {
            int latestStart = -1;
            String latestChrom = "";
            //add comments and sequence regions
            for (final String comment : codec.getComments()) {
                writer.addComment(comment);
            }

            for (final SequenceRegion sequenceRegion : codec.getSequenceRegions()) {
                writer.addSequenceRegionDirective(sequenceRegion);
            }

            //add features
            for (final Gff3Feature feature : sorter) {
                if (latestStart>=0 && (feature.getStart()>latestStart || !feature.getContig().equals(latestChrom))) {
                    writer.addFlushDirective();
                    latestStart = 0;
                }
                writer.addFeature(feature);
                final String featureID = feature.getID();
                if (featureID != null) {
                    latestStart = Math.max(latestStart, latestStartMap.get(feature.getID()));
                } else {
                    //if featureID was null we get latestStart for it based on its parents, since they are the only features it could be linked with
                    final List<String> parentIDs = feature.getAttribute("Parent");
                    latestStart = Math.max(latestStart, parentIDs.stream().map(latestStartMap::get).max(Integer::compareTo).orElse(feature.getStart()));
                }
                latestChrom = feature.getContig();
                progressWrite.record(feature.getContig(), feature.getStart());
            }
        } catch (final IOException ex) {
            throw new PicardException("Error opening  " + OUTPUT + " to write to", ex);
        }

        return 0;
    }

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
            //new LineIteratorImpl(new SynchronousLineReader(is));
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
