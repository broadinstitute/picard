package picard.annotation;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.util.ParsingUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

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

    private final Map<String, Integer> positionMap = new HashMap<>();
    private final Map<String, Integer> latestStartMap = new HashMap<>();
    private final Map<String, String> unfoundParentsMap = new HashMap<>();

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

        final Object header;

        try {
            final String path = INPUT.getAbsolutePath();
            final PositionalBufferedStream pbs;
            final InputStream inputStream = ParsingUtils.openInputStream(path, null);
            if (IOUtil.hasBlockCompressedExtension(path)) {
                // Gzipped -- we need to buffer the GZIPInputStream methods as this class makes read() calls,
                // and seekableStream does not support single byte reads
                final InputStream is = new GZIPInputStream(new BufferedInputStream(inputStream, 512000));
                pbs = new PositionalBufferedStream(is, 1000);  // Small buffer as this is buffered already.
            } else {
                pbs = new PositionalBufferedStream(inputStream, 512000);
            }

            final LineIterator lineIterator = codec.makeSourceFromStream(pbs);

            header = codec.readHeader(lineIterator).getHeaderValue();
            while (!codec.isDone(lineIterator)) {
                final Gff3Feature feature = codec.decodeShallow(lineIterator);
                if (feature != null) {
                    sorter.add(feature);
                    progressRead.record(feature.getContig(), feature.getStart());
                    positionMap.put(feature.getID(), feature.getStart());
                    final Integer currentLatestStart = latestStartMap.getOrDefault(feature.getID(),0);
                    final String parentIDAttribute = feature.getAttributes().get("Parent");
                    final List<String> parentIDs = parentIDAttribute != null? ParsingUtils.split(parentIDAttribute, ',') : new ArrayList<>();
                    final Integer parentsLatestStart = parentIDs.stream().map(id -> positionMap.get(id)).max(Integer::compareTo).orElse(0);

                    latestStartMap.put(feature.getID(), Math.max(currentLatestStart, parentsLatestStart));
                    for (final String parentID : parentIDs) {
                        final Integer parentCurrentLatestStart = latestStartMap.getOrDefault(parentID,0);
                        latestStartMap.put(parentID, Math.max(parentCurrentLatestStart, feature.getStart()));
                    }
                }
            }

        } catch(final IOException e) {
            throw new PicardException("Error reading file " + INPUT, e);
        }

        final ProgressLogger progressWrite = new ProgressLogger(log, (int) 1e4, "Wrote");
        try(final Gff3Writer writer = new Gff3Writer(OUTPUT, header)) {
            int latestStart = -1;
            String latestChrom = "";
            for (final Gff3Feature feature : sorter) {
                if (latestStart>=0 && (feature.getStart()>latestStart || !feature.getContig().equals(latestChrom))) {
                    writer.addFlushDirective();
                    latestStart = 0;
                }
                writer.addFeature(feature);
                latestStart = Math.max(latestStart, latestStartMap.get(feature.getID()));
                latestChrom = feature.getContig();
                progressWrite.record(feature.getContig(), feature.getStart());
            }
        } catch (final IOException ex) {
            throw new PicardException("Error opening  " + OUTPUT + " to write to", ex);
        }

        return 0;
    }

    class Gff3SortingCollectionCodec implements SortingCollection.Codec<Gff3Feature> {
        PrintStream outStream;
        LineIterator lineIteratorIn;
        Gff3Codec gff3Codec;
        Gff3Writer gff3Writer;


        Gff3SortingCollectionCodec() {
            gff3Codec = new Gff3Codec();
        }

        @Override
        public SortingCollection.Codec<Gff3Feature> clone() {
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
        public Gff3Feature decode() {
            while (!gff3Codec.isDone(lineIteratorIn)) {
                final Gff3Feature feature = gff3Codec.decodeShallow(lineIteratorIn);
                if (feature == null) {
                    continue;
                }
                return feature;
            }
            return null;
        }

        @Override
        public void encode(final Gff3Feature val) {
            gff3Writer.addFeature(val);
        }
    }
}
