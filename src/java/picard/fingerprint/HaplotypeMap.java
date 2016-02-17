/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.fingerprint;

import picard.PicardException;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.StringLineReader;

import java.io.*;
import java.util.*;

/**
 * A collection of metadata about Haplotype Blocks including multiple in memory "indices" of the data
 * to make it easy to query the correct HaplotypeBlock or Snp by snp names, positions etc. Also has the
 * ability to read and write itself to and from files.
 *
 * @author Tim Fennell / Kathleen Tibbetts
 */
public class HaplotypeMap {
    private final List<HaplotypeBlock> haplotypeBlocks = new ArrayList<HaplotypeBlock>();
    private final Map<Snp, HaplotypeBlock> haplotypesBySnp = new HashMap<Snp, HaplotypeBlock>();
    private final Map<String, HaplotypeBlock> haplotypesBySnpName = new HashMap<String, HaplotypeBlock>();
    private final Map<String, HaplotypeBlock> haplotypesBySnpLocus = new HashMap<String, HaplotypeBlock>();
    private final Map<String,Snp> snpsByPosition = new HashMap<String,Snp>();
    private final IntervalList intervals;
    private final SAMFileHeader header;

    /**
     * Constructs a HaplotypeMap from the provided file.
     */
    public HaplotypeMap(final File file) {
        BufferedReader in = null;
        try {
            in = new BufferedReader(new InputStreamReader(IOUtil.openFileForReading(file)));
            // Setup a reader and parse the header
            final StringBuilder builder = new StringBuilder(4096);
            String line = null;

            while ((line = in.readLine()) != null) {
                if (line.startsWith("@")) {
                    builder.append(line).append('\n');
                }
                else {
                    break;
                }
            }

            if (builder.length() == 0) {
                throw new IllegalStateException("Haplotype map file must contain header: " + file.getAbsolutePath());
            }

            this.header = new SAMTextHeaderCodec().decode(new StringLineReader(builder.toString()), "BufferedReader");
            this.intervals = new IntervalList(header);

            // Then read in the file
            final FormatUtil format = new FormatUtil();
            final List<HaplotypeMapFileEntry> entries = new ArrayList<HaplotypeMapFileEntry>();
            final Map<String, HaplotypeBlock> anchorToHaplotype = new HashMap<String, HaplotypeBlock>();

            do {
                if (line.trim().isEmpty()) continue; // skip over blank lines
                if (line.startsWith("#")) continue;      // ignore comments/headers

                // Make sure we have the right number of fields
                final String[] fields = line.split("\\t");
                if (fields.length < 6 || fields.length > 8) {
                    throw new PicardException("Invalid haplotype map record contains " +
                                              fields.length + " fields: " + line);
                }

                // Then parse them out
                final String chrom = fields[0];
                final int pos = format.parseInt(fields[1]);
                final String name = fields[2];
                final byte major =  (byte)fields[3].charAt(0);
                final byte minor =  (byte)fields[4].charAt(0);
                final double maf = format.parseDouble(fields[5]);
                final String anchor = fields.length > 6 ? fields[6] : null;
                final String fpPanels = fields.length > 7 ? fields[7] : null;
                List<String> panels = null;
                if (fpPanels != null) {
                    panels = new ArrayList<String>();
                    for (final String panel : fpPanels.split(",")) {
                        panels.add(panel);
                    }
                }

                // If it's the anchor snp, start the haplotype
                if (anchor == null || anchor.trim().equals("") || name.equals(anchor)) {
                    final HaplotypeBlock type = new HaplotypeBlock(maf);
                    type.addSnp(new Snp(name, chrom, pos, major, minor, maf, panels));
                    anchorToHaplotype.put(name, type);
                }
                else {  // Otherwise save it for later
                    final HaplotypeMapFileEntry entry = new HaplotypeMapFileEntry(
                            chrom, pos, name, major, minor, maf, anchor, panels);
                    entries.add(entry);
                }
            }
            while ((line = in.readLine()) != null);

            // Now, go through and add all the anchored snps
            for (final HaplotypeMapFileEntry entry : entries) {
                final HaplotypeBlock block = anchorToHaplotype.get(entry.anchorSnp);

                if (block == null) {
                    throw new PicardException("No haplotype found for anchor snp " + entry.anchorSnp);
                }

                block.addSnp(new Snp(entry.snpName, entry.chromosome, entry.position,
                                 entry.majorAllele, entry.minorAllele,
                                 entry.minorAlleleFrequency, entry.panels));
            }

            // And add them all
            for (final HaplotypeBlock block : anchorToHaplotype.values()) {
                addHaplotype(block);
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Error parsing haplotype map.", ioe);
        }
        finally {
            if (in != null) {
                try { in.close(); } catch (Exception e) { /* do nothing */ }
            }
        }
    }

    /** Constructs an empty HaplotypeMap using the provided SAMFileHeader's sequence dictionary. */
    public HaplotypeMap(final SAMFileHeader header) {
        this.header = header;
        this.intervals = new IntervalList(header);
    }

    /**
     * Adds a HaplotypeBlock to the map and updates all the relevant caches/indices.
     */
    public void addHaplotype(final HaplotypeBlock haplotypeBlock) {
        this.haplotypeBlocks.add(haplotypeBlock);

        for (final Snp snp : haplotypeBlock.getSnps()) {
            this.haplotypesBySnp.put(snp, haplotypeBlock);
            this.haplotypesBySnpName.put(snp.getName(), haplotypeBlock);
            this.haplotypesBySnpLocus.put(toKey(snp.getChrom(), snp.getPos()), haplotypeBlock);
            this.snpsByPosition.put(toKey(snp.getChrom(), snp.getPos()), snp);
            this.intervals.add(new Interval(snp.getChrom(), snp.getPos(), snp.getPos(), false, snp.getName()));
        }
    }

    /** Queries a HaplotypeBlock by Snp object. Returns NULL if none found. */
    public HaplotypeBlock getHaplotype(final Snp snp) {
        return this.haplotypesBySnp.get(snp);
    }

    /** Queries a HaplotypeBlock by Snp name. Returns NULL if none found. */
    public HaplotypeBlock getHaplotype(final String snpName) {
        return this.haplotypesBySnpName.get(snpName);
    }

    /** Queries a HaplotypeBlock by Snp chromosome and position. Returns NULL if none found. */
    public HaplotypeBlock getHaplotype(final String chrom, final int pos) {
        return this.haplotypesBySnpLocus.get(toKey(chrom, pos));
    }

    /** Returns an unmodifiable collection of all the haplotype blocks in the map. */
    public List<HaplotypeBlock> getHaplotypes() {
        return Collections.unmodifiableList(this.haplotypeBlocks);
    }

    /** Queries a Snp by chromosome and position. Returns NULL if none found. */
    public Snp getSnp(final String chrom, final int pos) {
        return this.snpsByPosition.get(toKey(chrom, pos));
    }

    /** Returns an unmodifiable collection of all SNPs in all Haplotype blocks. */
    public Set<Snp> getAllSnps() {
        return Collections.unmodifiableSet(haplotypesBySnp.keySet());
    }

    /** Returns an IntervalList with an entry for every SNP in every Haplotype in the map. */
    public IntervalList getIntervalList() {
        this.intervals.sort(); // TODO: should probably do this elsewhere
        return this.intervals;
    }

    private String toKey(final String chrom, final int pos) {
        return chrom + ":" + pos;
    }

    /**
     * Returns a copy of this haplotype map that excludes haplotypes on the chromosomes provided.
     * @param chroms a set of zero or more chromosome names
     */
    public HaplotypeMap withoutChromosomes(final Set<String> chroms) {
        final HaplotypeMap out = new HaplotypeMap(getHeader());
        for (final HaplotypeBlock block : this.haplotypeBlocks) {
            if (!chroms.contains(block.getFirstSnp().getChrom())) {
                out.addHaplotype(block);
            }
        }

        return out;
    }

    /** Writes out a HaplotypeMap file with the contents of this map. */
    public void writeToFile(final File file) {
        try {
            final BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(file)));
            final FormatUtil format = new FormatUtil();

            // Write out the header
            if (this.header != null) {
                final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                codec.encode(out, this.header);
            }

            // Write the header for the entries.
            out.write("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS");
            out.newLine();

            final List<HaplotypeMapFileEntry> entries = new ArrayList<HaplotypeMapFileEntry>();
            for (final HaplotypeBlock block : this.getHaplotypes()) {
                String anchor = null;
                final SortedSet<Snp> snps = new TreeSet<Snp>(block.getSnps());

                for (final Snp snp : snps) {
                    entries.add(new HaplotypeMapFileEntry(snp.getChrom(), snp.getPos(), snp.getName(),
                            snp.getAllele1(), snp.getAllele2(), snp.getMaf(), anchor, snp.getFingerprintPanels()));

                    if (anchor == null) {
                        anchor = snp.getName();
                    }
                }
            }

            Collections.sort(entries);
            for (final HaplotypeMapFileEntry entry : entries) {
                out.write(entry.chromosome + "\t");
                out.write(format.format(entry.position) + "\t");
                out.write(entry.snpName + "\t");
                out.write((char)entry.majorAllele + "\t");
                out.write((char)entry.minorAllele + "\t");
                out.write(format.format(entry.minorAlleleFrequency) + "\t");
                if (entry.anchorSnp != null) {
                    out.write(entry.anchorSnp);
                }
                out.write("\t");
                if (entry.getPanels() != null) {
                    out.write(entry.getPanels());
                }
                out.newLine();
            }
            out.flush();
            out.close();
        }
        catch (IOException ioe) {
            throw new PicardException("Error writing out maplotype map to file: " + file.getAbsolutePath(), ioe);
        }
    }

    public SAMFileHeader getHeader() { return header; }

    /** Class used to represent all the information for a row in a haplotype map file, used in reading and writing. */
    private class HaplotypeMapFileEntry implements Comparable {
        private final String chromosome;
        private final int position;
        private final String snpName;
        private final byte majorAllele;
        private final byte minorAllele;
        private final double minorAlleleFrequency;
        private final String anchorSnp;
        private final List<String> panels;

        public HaplotypeMapFileEntry(final String chrom, final int pos, final String name,
                                     final byte major, final byte minor, final double maf,
                                     final String anchorSnp, final List<String> fingerprintPanels) {
            this.chromosome = chrom;
            this.position = pos;
            this.snpName = name;
            this.majorAllele = major;
            this.minorAllele = minor;
            this.minorAlleleFrequency = maf;
            this.anchorSnp = anchorSnp;

            // Always sort the list of panels so they are in a predictable order
            this.panels = new ArrayList<String>();
            if (fingerprintPanels != null) {
                this.panels.addAll(fingerprintPanels);
                Collections.sort(this.panels);
            }            
        }

        public String getPanels() {
            if (panels == null) return "";
            final StringBuilder sb = new StringBuilder();

            for (final String panel : panels) {
                if (sb.length() > 0) sb.append(",");
                sb.append(panel);
            }

            return sb.toString();
        }

        public int compareTo(final Object o) {
            final HaplotypeMapFileEntry that = (HaplotypeMapFileEntry) o;
            int diff = header.getSequenceIndex(this.chromosome) - header.getSequenceIndex(that.chromosome);
            if (diff != 0) return diff;

            diff = this.position - that.position;
            if (diff != 0) return diff;

            diff = this.snpName.compareTo(that.snpName);
            if (diff != 0) return diff;

            diff = this.majorAllele - that.majorAllele;
            if (diff != 0) return diff;

            diff = this.minorAllele - that.minorAllele;
            if (diff != 0) return diff;

            diff = Double.compare(this.minorAlleleFrequency, that.minorAlleleFrequency);
            if (diff != 0) return diff;

            if (this.anchorSnp != null) {
                if (that.anchorSnp != null) {
                    diff = this.anchorSnp.compareTo(that.anchorSnp);
                }
                else {
                    diff = 1;
                }
            }
            else {
                if (that.anchorSnp != null) {
                    diff = -1;
                }
                else {
                    diff = 0;
                }

            }
            if (diff != 0) return diff;

            final String p1 = this.getPanels();
            final String p2 = that.getPanels();
            if (p1 != null) {
                if (p2 != null) {
                    return p1.compareTo(p2);
                }
                return 1;
            }
            else if (p2 != null) {
                return -1;
            }
            else {
                return 0;
            }
        }
    }
}
