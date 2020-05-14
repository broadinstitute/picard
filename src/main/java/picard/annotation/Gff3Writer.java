package picard.annotation;



import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.util.ParsingUtils;
import picard.PicardException;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

//to-do: move this class to htsjdk
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
        try {
            final String lineNoAttributes = String.join("\t",
                    URLEncoder.encode(feature.getContig(), "UTF-8"),
                    URLEncoder.encode(feature.getSource(), "UTF-8"),
                    URLEncoder.encode(feature.getType(), "UTF-8"),
                    Integer.toString(feature.getStart()),
                    Integer.toString(feature.getEnd()),
                    feature.getScore() < 0 ? "." : Double.toString(feature.getScore()),
                    feature.getStrand().toString(),
                    feature.getPhase() < 0 ? "." : Integer.toString(feature.getPhase())
            );
            final List<String> attributesStrings = feature.getAttributes().entrySet().stream().map(e -> {
                try {
                    return String.join("=", new String[]{encodeForNinthColumn(e.getKey()), encodeForNinthColumn(e.getValue())});
                    } catch (final URISyntaxException ex) {
                    throw new PicardException("Exception writing out gff",ex);
                    }
                }
            ).collect(Collectors.toList());
            final String attributesString = attributesStrings.isEmpty() ? "." : String.join(";", attributesStrings);

            final String lineString = lineNoAttributes + "\t" + attributesString;
            out.println(lineString);
        } catch(final UnsupportedEncodingException ex) {
            throw new PicardException("Exception writing out gff",ex);
        }
    }

    private String encodeForNinthColumn(final String decodedString) throws URISyntaxException {
        //in the ninth column of Gff certain characters have special meaning
        final List<String> splitString = ParsingUtils.split(decodedString, ',');
        final List<String> encodedSplitString = new ArrayList<>();
        for (final String string : splitString) {
            final URI uri = new URI(string);
            encodedSplitString.add(uri.toASCIIString());
        }

        return String.join(",", encodedSplitString);
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
