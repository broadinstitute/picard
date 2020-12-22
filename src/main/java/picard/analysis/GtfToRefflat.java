package picard.analysis;

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

public class GtfToRefflat {

    private final File GTF;


    private static final String COLUMN_DELIMITER = "\t";
    private static final String COORDINATE_DELIMITER = ",";
    private static final String NEW_LINE_DELIMITER = "\n";
    private static final String ATTRIBUTE_DELIMITER = " ";

    private File refflat = null;

    private List<String> rows = new ArrayList<>();
    private String transcriptId = "";
    private List<Integer> exonStarts = new ArrayList<>();
    private List<Integer> exonEnds = new ArrayList<>();
    private int codingStart = -1;
    private int codingEnd = -1;
    private String gene_id = "";
    private String chromosome = "";

    private Strand strand = null;
    private String type = "";

    public GtfToRefflat(File GTF) {
        this.GTF = GTF;

        this.doWork();
    }


    protected int doWork() {
        if (GTF != null) {

            // convert the Gtf to a Gff3 and create a Gff3Feature
            final File GFF3 = convertToGFF3(GTF);
            final AbstractFeatureReader<Gff3Feature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(GFF3.getAbsolutePath(), null, new Gff3Codec(), false);

            String currentTranscriptId = "";

            try {

                boolean hasExon = false;
                String ignore_id = "";
                boolean has_stopCodon = false;

                // iterate through each line in the Gff3Feature
                for (final Gff3Feature feature : reader.iterator()) {

                    // if the line has no transcipt_id, move onto the next line
                    if(!feature.getAttribute("transcript_id").isEmpty()) {
                        currentTranscriptId = feature.getAttribute("transcript_id").get(0);
                    } else {
                        continue;
                    }

                    // Since information is grouped by transcipt_id, once the transcript_id is different
                    // create a line in the refflat file with the collected information
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

                    if (type.equals("exon")) {
                        exonStarts.add(feature.getStart() - 1);
                        exonEnds.add(feature.getEnd());
                        hasExon = true;
                    } else if (hasExon == false) {
                        if (!type.equals("start_codon")) {
                            exonEnds.add(feature.getEnd());
                        }
                        if (!type.equals("stop_codon")) {
                            exonStarts.add(feature.getStart() - 1);
                        }
                    }

                    if (strand.equals(Strand.POSITIVE)) {
                        if (type.equals("start_codon")) {
                            codingStart = feature.getStart();
                        }
                        if (type.equals("stop_codon")) {
                            codingEnd = feature.getEnd();
                        }
                        if (type.equals("cds") && codingStart == -1) {
                            codingStart = feature.getStart();
                        }
                        if (type.equals("cds")) {
                            codingEnd = feature.getEnd();
                        }
                    } else {
                        if (type.equals("stop_codon")) {
                            codingStart = feature.getStart();
                        }
                        if (type.equals("start_codon")) {
                            codingEnd = feature.getEnd();
                            has_stopCodon = true;
                        }
                        if (type.equals("cds") && codingStart == -1) {
                            codingStart = feature.getStart();
                        }
                        if (type.equals("cds") && !has_stopCodon) {
                            codingEnd = feature.getEnd();
                        }
                    }

                    transcriptId = currentTranscriptId;
                }

                // add last row, format them into refflat format, and create the refflat file
                this.addRow();
                String data = String.join(NEW_LINE_DELIMITER, rows);
                refflat = writeToFile("new" + GTF.getName().substring(0, GTF.getName().length() - 4) + ".refflat", data);

            } catch (Exception e) {
                System.out.println("There was an error while converting the given GFT to a refFlat for CollectRnaSeqMetrics. " +
                        "Make sure the GTF file is tab separated.");
            }
        }
        return 0;
    }

    // calculate start and end variables and the exon count, format the variables into a refflat line
    private void addRow() {
        Collections.sort(exonStarts);
        exonStarts = exonStarts.stream().distinct().collect(Collectors.toList());
        Collections.sort(exonEnds);
        exonEnds = exonEnds.stream().distinct().collect(Collectors.toList());

        codingStart = codingStart == -1 ? exonStarts.get(exonStarts.size() - 1) : codingStart - 1;
        codingEnd = codingEnd == -1 ? exonEnds.get(exonEnds.size() - 1) : codingEnd;

        rows.add(String.join(
                COLUMN_DELIMITER,
                gene_id, transcriptId, chromosome, strand.toString(),
                Integer.toString(exonStarts.get(0)), Integer.toString(exonEnds.get(exonEnds.size() - 1)),
                Integer.toString(codingStart), Integer.toString(codingEnd),
                Integer.toString(exonStarts.size()),
                exonStarts.stream().map(Object::toString).collect(Collectors.joining(COORDINATE_DELIMITER)
                ), exonEnds.stream().map(Object::toString).collect(Collectors.joining(COORDINATE_DELIMITER))));

        resetVariables();
    }

    // reset all the start end and count variables for the next line in the Gff3Feature
    public void resetVariables() {
        exonStarts = new ArrayList<>();
        exonEnds = new ArrayList<>();
        codingStart = -1;
        codingEnd = -1;
    }

    // return the Refflat file
    public File getRefflat() {
        return refflat;
    }

    // create a file and write data in it
    public File writeToFile(String fileName, String data) {
        File newFile = new File(fileName);
        FileWriter fr = null;
        try {
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

    // convert the Gtf file to a Gff3 file
    public File convertToGFF3(File gtf) {
        List<String> rows = new ArrayList<>();

        // read in the gtf file and if the current line is not a comment or blank, convert it to gff3 format
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

        return writeToFile("new" + GTF.getName().substring(0, GTF.getName().length() - 4) + ".gff3", data);
    }

    // format the given gtf line to gff3 format. Everything is the same except the attribute column
    public String useGff3Syntax(String row) {
        String[] values = row.split(COLUMN_DELIMITER);
        String[] attributes = values[values.length - 1].split(ATTRIBUTE_DELIMITER);
        List<String> resultValues = new ArrayList<>();
        List<String> resultAttributes = new ArrayList<>();
        for (int i = 0; i < values.length - 1; i++) {
            resultValues.add(values[i]);
        }
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
                case "gene_name":
                    resultAttributes.add("Name=" + attributes[i + 1].replace("\"", ""));
                    i++;
                    break;
                default:
                    break;
            }
        }

        String newAttributes = String.join("", resultAttributes);
        String newRow = String.join(COLUMN_DELIMITER, resultValues);
        newAttributes = StringUtils.chop(newAttributes);

        String newValues = newRow.concat(COLUMN_DELIMITER + newAttributes);

        return newValues;
    }

}
