package net.sf.picard.pedigree;

import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.util.RuntimeIOException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;

/**
 * Represents a .ped file of family information as documented here:
 *    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
 *
 * Stores the information in memory as a map of individualId -> Pedigree information for that individual
 */
public class PedFile extends TreeMap<String,PedTrio> {
    private static final Log log = Log.getInstance(PedFile.class);
    static final Pattern WHITESPACE = Pattern.compile("\\s+");

    /** Adds a trio to the PedFile keyed by the individual id. */
    public void add(final PedTrio trio) {
        put(trio.getIndividualId(), trio);
    }

    /**
     * Writes a set of pedigrees out to disk.
     */
    public void write(final File file) {
        IoUtil.assertFileIsWritable(file);
        final BufferedWriter out = IoUtil.openFileForBufferedWriting(file);

        try {
            for (final PedTrio trio : values()) {
                out.write(trio.getFamilyId());
                out.write("\t");
                out.write(trio.getIndividualId());
                out.write("\t");
                out.write(trio.getPaternalId());
                out.write("\t");
                out.write(trio.getMaternalId());
                out.write("\t");
                out.write(String.valueOf(trio.getSex().toCode()));
                out.write("\t");
                out.write(trio.getPhenotype().toString());
                out.newLine();
            }

            out.close();
        }
        catch (IOException ioe) {
            throw new RuntimeIOException("IOException while writing to file " + file.getAbsolutePath(), ioe);
        }
    }

    /**
     * Attempts to read a pedigree file into memory.
     */
    public static PedFile fromFile(final File file) {
        final PedFile pedfile = new PedFile();

        IoUtil.assertFileIsReadable(file);
        for (final String line : IoUtil.readLines(file)) {
            final String[] fields = WHITESPACE.split(line);
            if (fields.length != 6) {
                log.error("Ped file line contained invalid number of fields, skipping: " + line);
                continue;
            }

            final PedTrio trio = new PedTrio(fields[0],
                                             fields[1],
                                             fields[2],
                                             fields[3],
                                             Sex.fromCode(Integer.parseInt(fields[4])),
                                             fields[5].contains(".") ? Double.parseDouble(fields[5]) : Integer.parseInt(fields[5])
                                            );
            pedfile.add(trio);
        }

        return pedfile;
    }

    /**
     * Scans through the pedigrees and removes all entries that do not have both paternal and maternal ids set.
     */
    public PedFile removeIncompleteTrios() {
        final Iterator<Map.Entry<String,PedTrio>> iterator = entrySet().iterator();

        while (iterator.hasNext()) {
            if (!iterator.next().getValue().hasBothParents()) iterator.remove();
        }

        return this;
    }
}
