package picard.annotation;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
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
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

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

    @Override
    protected int doWork() {
        if (!(new Gff3Codec().canDecode(INPUT.toString()))) {
            throw new IllegalArgumentException("Input file " + INPUT + " cannot be read by Gff3Codec");
        }

        SortingCollection<Gff3FeatureWithFlushableLocation> sorter = SortingCollection.newInstance(
                Gff3FeatureWithFlushableLocation.class,
                new Gff3SortingCollectionCodec(),
                (Gff3FeatureWithFlushableLocation f1, Gff3FeatureWithFlushableLocation f2) -> f1.feature.getContig().compareTo(f2.feature.getContig()) == 0 ? f1.feature.getStart() - f2.feature.getStart() : f1.feature.getContig().compareTo(f2.feature.getContig()),
                500000
        );

        final ProgressLogger progressRead = new ProgressLogger(log, (int) 1e4, "Read");

        final Object header;

        try (final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(INPUT.toString(), null, new Gff3Codec(), false)) {

            header = reader.getHeader();
            for (final Gff3Feature feature : reader.iterator()) {
                sorter.add(new Gff3FeatureWithFlushableLocation(feature));
                progressRead.record(feature.getContig(), feature.getStart());
            }
        } catch(final IOException e) {
            throw new PicardException("Error reading file " + INPUT, e);
        }

        final ProgressLogger progressWrite = new ProgressLogger(log, (int) 1e4, "Wrote");
        try(final Gff3Writer writer = new Gff3Writer(OUTPUT, header)) {
            int latestStart = -1;
            for (final Gff3FeatureWithFlushableLocation featureWithFlushableLocation : sorter) {
                if (latestStart>=0 && featureWithFlushableLocation.feature.getStart()>latestStart) {
                    writer.addFlushDirective();
                }
                writer.addFeature(featureWithFlushableLocation.feature);
                latestStart = Math.max(latestStart, featureWithFlushableLocation.latestStart);
                progressWrite.record(featureWithFlushableLocation.feature.getContig(), featureWithFlushableLocation.feature.getStart());
            }
        } catch (final IOException ex) {
            throw new PicardException("Error opening  " + OUTPUT + " to write to", ex);
        }

        return 0;
    }

    class Gff3SortingCollectionCodec implements SortingCollection.Codec<Gff3FeatureWithFlushableLocation> {
        PrintStream outStream;
        LineIterator lineIteratorIn;
        Gff3Codec gff3Codec;
        Gff3Writer gff3Writer;


        Gff3SortingCollectionCodec() {
            gff3Codec = new Gff3Codec();
        }

        @Override
        public SortingCollection.Codec<Gff3FeatureWithFlushableLocation> clone() {
            return new Gff3SortingCollectionCodec();
        }

        @Override
        public void setOutputStream(final OutputStream os) {
            outStream = new PrintStream(os);
            gff3Writer = new Gff3Writer(outStream);
        }

        @Override
        public void setInputStream(final InputStream is) {
            lineIteratorIn = gff3Codec.makeSourceFromStream(is);
            //new LineIteratorImpl(new SynchronousLineReader(is));
        }

        @Override
        public Gff3FeatureWithFlushableLocation decode() {
            if (!lineIteratorIn.hasNext()) {
                return null;
            }
            try {
                final Gff3Feature feature = gff3Codec.decode(lineIteratorIn);
                final int latestStart = Integer.parseInt(lineIteratorIn.next());
                return new Gff3FeatureWithFlushableLocation(feature, latestStart);
            } catch (IOException e) {
                throw new PicardException("Gff3SortingCollection error while reading from disk",e);
            }
        }

        @Override
        public void encode(final Gff3FeatureWithFlushableLocation val) {
            gff3Writer.addFeature(val.feature);
            gff3Writer.addFlushDirective();
            //write start position of latest directly linked feature
            outStream.println(val.latestStart);
        }
    }

    private static class Gff3FeatureWithFlushableLocation {
        private final Gff3Feature feature;
        private final int latestStart;

        Gff3FeatureWithFlushableLocation(final Gff3Feature feature) {
            this.feature = feature;

            int latestStartParents = feature.getParents().stream().mapToInt(Gff3Feature::getStart).max().orElse(0);
            int latestStartChildren = feature.getChildren().stream().mapToInt(Gff3Feature::getStart).max().orElse(0);
            int latestStartCoFeatures = feature.getCoFeatures().stream().mapToInt(Gff3Feature::getStart).max().orElse(0);

            latestStart = Math.max(latestStartChildren, Math.max(latestStartParents, latestStartCoFeatures));
        }

        Gff3FeatureWithFlushableLocation(final Gff3Feature feature, final int latestStart) {
            this.feature = feature;
            this.latestStart = latestStart;
        }

    }

}
