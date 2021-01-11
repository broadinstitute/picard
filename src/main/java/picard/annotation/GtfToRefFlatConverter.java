package picard.annotation;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.lang.StringUtils;
import picard.PicardException;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.annotation.Strand;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class GtfToRefFlatConverter {

    private final File GTF;

    private File refFlat = null;

    private static final String COLUMN_DELIMITER = "\t";
    private static final String COORDINATE_DELIMITER = ",";
    private static final String NEW_LINE_DELIMITER = "\n";
    private static final String ATTRIBUTE_DELIMITER = " ";

    private String gene_id = "";
    private String transcriptId = "";
    private String chromosome = "";
    private Strand strand = null;
    private String type = "";
    private int cdsStart = -1;
    private int cdsEnd = -1;
    private List<Integer> exonStarts = new ArrayList<>();
    private List<Integer> exonEnds = new ArrayList<>();

    private final List<String> rows = new ArrayList<>();

    int minExonStart = Integer.MAX_VALUE;
    int maxExonEnd = Integer.MIN_VALUE;
    boolean hasExon = false;

    /**
     * @param GTF             Gene annotations in GTF form
     */
    public GtfToRefFlatConverter(File GTF) {
        this.GTF = GTF;

        this.doWork();
    }

    private int doWork() {
        if (GTF != null) {

            // convert the Gtf to a Gff3 and create a Gff3Feature
            final File GFF3 = convertToGFF3(GTF);
            final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(GFF3.getAbsolutePath(), null, new Gff3Codec(), false);

            String currentTranscriptId = "";
            String ignore_id = "";
            boolean has_stopCodon = false;

            try {

                for (final Gff3Feature feature : reader.iterator()) {

                    // if the line has no transcript_id, move onto the next line
                    if(feature.getAttribute("transcript_id").isEmpty()) {
                        continue;
                    }

                    currentTranscriptId = feature.getAttribute("transcript_id").get(0);

                    // Since information is grouped by transcipt_id, once the transcript_id is different
                    // create a row with the collected information for the refFlat
                    if (!transcriptId.equals(currentTranscriptId) && !transcriptId.equals("")
                            && !ignore_id.equals(transcriptId)) {
                        this.addRow();
                    } else if (strand != feature.getStrand() && strand != null) {
                        System.out.println("Error: all group members must be on the same strand");
                        ignore_id = currentTranscriptId;
                        resetVariables();
                        continue;
                    }

                    gene_id = feature.getID();
                    chromosome = feature.getContig();
                    strand = feature.getStrand();
                    type = feature.getType().toLowerCase();

                    int newStart = feature.getStart() - 1;
                    int newEnd = feature.getEnd();

                    calculateExonLists(newStart, newEnd);

                    if (strand.equals(Strand.POSITIVE)) {
                        if (type.equals("start_codon")) {
                            cdsStart = newStart + 1;
                        }
                        if (type.equals("stop_codon")) {
                            cdsEnd = newEnd;
                            has_stopCodon = true;
                        }
                    } else {
                        if (type.equals("stop_codon")) {
                            cdsStart = newStart + 1;
                        }
                        if (type.equals("start_codon")) {
                            cdsEnd = newEnd;
                            has_stopCodon = true;
                        }
                    }
                    if (type.equals("cds") && cdsStart == -1) {
                        cdsStart = newStart + 1;
                    }
                    if (type.equals("cds") && !has_stopCodon) {
                        cdsEnd = newEnd;
                    }

                    transcriptId = currentTranscriptId;
                }

                // add last row, format them into refFlat format, and create the refFlat file
                this.addRow();
                String data = String.join(NEW_LINE_DELIMITER, rows);
                refFlat = writeToFile(GTF.getName(), ".refflat", data);

            } catch (Exception e) {
                System.out.println("There was an error while converting the given GFT to a refFlat for CollectRnaSeqMetrics. " +
                        "Make sure the GTF file is tab separated.");
            }
        }
        return 0;
    }

    private void calculateExonLists(int newStart, int newEnd) {
        if (type.equals("exon")) {
            if (!hasExon && (!exonStarts.isEmpty() || !exonEnds.isEmpty())) {
                exonStarts.clear();
                exonEnds.clear();
            }
            exonStarts.add(newStart);
            exonEnds.add(newEnd);
            hasExon = true;
        } else if (hasExon == false) {

            if (maxExonEnd == Integer.MIN_VALUE && minExonStart == Integer.MAX_VALUE) {
                minExonStart = newStart;
                maxExonEnd = newEnd;
            } else {
                // if the intervals overlap, find the min start position and max end position
                if (newStart <= maxExonEnd || newEnd <= maxExonEnd) {
                    maxExonEnd = Math.max(newEnd, maxExonEnd);
                    minExonStart = Math.min(newStart, minExonStart);
                } else {
                    exonStarts.add(minExonStart);
                    exonEnds.add(maxExonEnd);

                    maxExonEnd = newEnd;
                    minExonStart = newStart;
                }
            }
        }
    }

    // calculate start and end variables and the exon count, format the variables into a refFlat line
    private void addRow() {

        if (!hasExon) {
            exonStarts.add(minExonStart);
            exonEnds.add(maxExonEnd);
        }

        Collections.sort(exonStarts);
        Collections.sort(exonEnds);

        cdsStart = cdsStart == -1 ? exonStarts.get(0) : cdsStart - 1;
        cdsEnd = cdsEnd == -1 ? exonEnds.get(exonEnds.size() - 1) : cdsEnd;

        rows.add(String.join(
                COLUMN_DELIMITER,
                gene_id, transcriptId, chromosome, strand.toString(),
                Integer.toString(exonStarts.get(0)), Integer.toString(exonEnds.get(exonEnds.size() - 1)),
                Integer.toString(cdsStart), Integer.toString(cdsEnd),
                Integer.toString(exonStarts.size()),
                exonStarts.stream().map(Object::toString).collect(Collectors.joining(COORDINATE_DELIMITER)
                ), exonEnds.stream().map(Object::toString).collect(Collectors.joining(COORDINATE_DELIMITER))));

        resetVariables();
    }

    // reset all the start end and count variables for the next line in the Gff3Feature
    private void resetVariables() {
        exonStarts = new ArrayList<>();
        exonEnds = new ArrayList<>();
        cdsStart = -1;
        cdsEnd = -1;
        maxExonEnd = Integer.MIN_VALUE;
        minExonStart = Integer.MAX_VALUE;
        transcriptId = "";
    }

    public File getRefFlat() {
        return refFlat;
    }

    private File writeToFile(String fileName, String suffix, String data) {
        File newFile;
        FileWriter fr = null;
        try {
            newFile = File.createTempFile(fileName, suffix);
            fr = new FileWriter(newFile);
            fr.write(data);
        } catch (IOException e) {
            throw new PicardException("Could not write to file " + fileName);
        } finally {
            try {
                fr.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return newFile;
    }

    private File convertToGFF3(File gtf) {
        List<String> rows = new ArrayList<>();

        try {
            Scanner myReader = new Scanner(gtf);
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                boolean isEmpty = data.matches("");
                boolean isComment = data.startsWith("#");
                if (!isEmpty && !isComment) {
                    rows.add(useGff3Syntax(data));
                }
            }
            myReader.close();
        } catch (FileNotFoundException e) {
            System.out.println("An error occurred. Could not find the gtf file.");
            e.printStackTrace();
        }
        String data = String.join(NEW_LINE_DELIMITER, rows);

        return writeToFile(GTF.getName(), ".gff3", data);
    }

    private String useGff3Syntax(String row) {
        String[] values = row.split(COLUMN_DELIMITER);
        String[] attributes = values[values.length - 1].split(ATTRIBUTE_DELIMITER);

        List<String> resultValues = new ArrayList<>(Arrays.asList(values).subList(0, values.length - 1));
        List<String> resultAttributes = new ArrayList<>();

        for (int i = 0; i < attributes.length; i++) {
            switch (attributes[i]) {
                case "gene_id":
                    resultAttributes.add("ID=" + attributes[i + 1].replace("\"", ""));
                    i++;
                    break;
                case "transcript_id":
                    resultAttributes.add("transcript_id=" + attributes[i + 1].replace("\"", ""));
                    i++;
                    break;
                default:
                    break;
            }
        }
        resultValues.add(StringUtils.chop(String.join("", resultAttributes)));

        return String.join(COLUMN_DELIMITER, resultValues);
    }
}
