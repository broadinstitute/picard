package picard.annotation;


import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import picard.PicardException;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class Gff3Writer implements Closeable {

    final private PrintStream out;
    final static String versionDirective = "##gff-version 3";
    final static String flushDirective = "###";
    final static String sequenceRegionDirective ="##sequence-region";

    Gff3Writer(final File file) throws FileNotFoundException {
        this(file, null);
    }

    Gff3Writer(final File file, final Object header) throws FileNotFoundException {
        this(new PrintStream(file), header);
    }

    Gff3Writer(final OutputStream stream) {
        this(stream, null);
    }

    Gff3Writer(final OutputStream stream, final Object header) {
        out = new PrintStream(stream);
        out.println(versionDirective);
        if (header != null) {
            if ((header instanceof List<?>)) {
                List<?> headerStrings = (ArrayList<?>) header;

                for (final Object headerLine : headerStrings) {
                    if ((headerLine instanceof String)) {
                        String headerLineString = (String) headerLine;
                        if (headerLineString.startsWith(versionDirective)) {
                            continue;
                        }

                        out.println(headerLineString);
                    } else {
                        throw new PicardException("Header passed to Gff3Writer is not a list of header lines");
                    }
                }
            } else {
                throw new PicardException("Header passed to Gff3Writer is not a list of header lines");
            }
        }
    }

    public void addFeature(final Gff3Feature feature) {
        final String lineNoAttributes = String.join("\t",
                feature.getContig(),
                feature.getSource(),
                feature.getType(),
                Integer.toString(feature.getStart()),
                Integer.toString(feature.getEnd()),
                ".",
                feature.getStrand().toString(),
                feature.getPhase()<0 ? "." : Integer.toString(feature.getPhase())
        );
        final List<String> attributesStrings = feature.getAttributes().entrySet().stream().map(e -> String.join("=", new String[] {e.getKey(), e.getValue()})).collect(Collectors.toList());
        final String attributesString = String.join(";", attributesStrings);

        final String lineString = lineNoAttributes + "\t" + attributesString;
        out.println(lineString);
    }

    public void addFlushDirective() {
        out.println(flushDirective);
    }

    public void addSequenceRegionDirective(final String contig, final int start, final int end) {
        out.println(String.join(" ", sequenceRegionDirective , contig, Integer.toString(start), Integer.toString(end)));
    }

    @Override
    public void close() {
        out.close();
    }
}
