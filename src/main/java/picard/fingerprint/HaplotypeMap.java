/*
 * The MIT License
 *
 * Copyright (c) 2010-2017 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import picard.PicardException;
import picard.vcf.VcfUtils;

import java.io.*;
import java.util.*;

/**
 * A collection of metadata about Haplotype Blocks including multiple in memory "indices" of the data
 * to make it easy to query the correct HaplotypeBlock or Snp by snp names, positions etc. Also has the
 * ability to read and write itself to and from files.
 *
 * @author Tim Fennell / Kathleen Tibbetts / Yossi Farjoun
 */
public class HaplotypeMap {
    public static final String HET_GENOTYPE_FOR_PHASING = "HetGenotypeForPhasing";
    public static final String SYNTHETIC_PHASESET_PREFIX = "Synthetic";
    public static final String PHASESET_PREFIX = "PhaseSet";

    private final List<HaplotypeBlock> haplotypeBlocks = new ArrayList<>();
    private final Map<Snp, HaplotypeBlock> haplotypesBySnp = new HashMap<>();
    private final Map<String, HaplotypeBlock> haplotypesBySnpName = new HashMap<>();
    private final Map<String, HaplotypeBlock> haplotypesBySnpLocus = new HashMap<>();
    private final Map<String,Snp> snpsByPosition = new HashMap<>();
    private IntervalList intervals;
    private SAMFileHeader header;



    /**
     * Constructs a HaplotypeMap from the provided file.
     */

    private void fromVcf(final File file) {
        try ( final VCFFileReader reader = new VCFFileReader(file, false)) {

            final SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();

            if (dict == null || dict.getSequences().isEmpty()) {
                throw new IllegalStateException("Haplotype map VCF file must contain header: " + file.getAbsolutePath());
            }

            initialize(new SAMFileHeader(dict));

            final Map<String, HaplotypeBlock> anchorToHaplotype = new HashMap<>();

            for (final VariantContext vc : reader) {

                if (vc.getNSamples() > 1) {
                    throw new IllegalStateException("Haplotype map VCF file must contain at most one sample: " + file.getAbsolutePath());
                }

                final Genotype gc = vc.getGenotype(0); // may be null
                final boolean hasGc = gc != null;

                if (vc.getAlternateAlleles().size() != 1) {
                    throw new IllegalStateException("Haplotype map VCF file must contain exactly one alternate allele per site: " + vc.toString());
                }

                if (!vc.isSNP()) {
                    throw new IllegalStateException("Haplotype map VCF file must contain only SNPs: " + vc.toString());
                }

                if (!vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                    throw new IllegalStateException("Haplotype map VCF Variants must have an '"+ VCFConstants.ALLELE_FREQUENCY_KEY + "' INFO field: " + vc.toString());
                }


                if (hasGc && gc.isPhased() && !gc.hasExtendedAttribute(VCFConstants.PHASE_SET_KEY)) {
                    throw new IllegalStateException("Haplotype map VCF Variants' genotypes that are phased must have a PhaseSet (" + VCFConstants.PHASE_SET_KEY+")" + vc.toString());
                }

                if (hasGc && gc.isPhased() && !gc.isHet()) {
                    throw new IllegalStateException("Haplotype map VCF Variants' genotypes that are phased must be HET" + vc.toString());
                }

                // Then parse them out
                final String chrom = vc.getContig();
                final int pos = vc.getStart();
                final String name = vc.getID();

                final byte ref = vc.getReference().getBases()[0];
                final byte var = vc.getAlternateAllele(0).getBases()[0];

                final double temp_maf = vc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0D);
                final boolean swapped = hasGc && !gc.getAllele(0).equals(vc.getReference());

                final byte major, minor;
                final double maf;

                if (swapped) {
                    major = var;
                    minor = ref;
                    maf = 1 - temp_maf;
                } else {
                    major = ref;
                    minor = var;
                    maf = temp_maf;
                }

                final String anchor = anchorFromVc(vc);

                // If it's the anchor snp, start the haplotype
                if (!anchorToHaplotype.containsKey(anchor)) {
                    final HaplotypeBlock newBlock = new HaplotypeBlock(maf);
                    anchorToHaplotype.put(anchor, newBlock);
                }
                final HaplotypeBlock block = anchorToHaplotype.get(anchor);
                block.addSnp(new Snp(name, chrom, pos, major, minor, maf, null));
            }

            // And add them all now that they are all ready.
            anchorToHaplotype.values().forEach(this::addHaplotype);
        }
    }

    static private String anchorFromVc(final VariantContext vc) {
        final Genotype genotype = vc.getGenotype(0);

        if (genotype == null || !genotype.hasExtendedAttribute(VCFConstants.PHASE_SET_KEY)) {
            return SYNTHETIC_PHASESET_PREFIX + "_" + vc.getContig() + "_" + vc.getStart();
        } else {
            return PHASESET_PREFIX + "_" + vc.getContig() + "_" + genotype.getExtendedAttribute(VCFConstants.PHASE_SET_KEY);
        }
    }

    private void fromHaplotypeDatabase(final File file) {

        BufferedReader in = null;
        try {
            in = new BufferedReader(new InputStreamReader(IOUtil.openFileForReading(file)));
            // Setup a reader and parse the header
            final StringBuilder builder = new StringBuilder(4096);
            String line = null;

            while ((line = in.readLine()) != null) {
                if (line.startsWith("@")) {
                    builder.append(line).append('\n');
                } else {
                    break;
                }
            }

            if (builder.length() == 0) {
                throw new IllegalStateException("Haplotype map file must contain header: " + file.getAbsolutePath());
            }

            final SAMFileHeader header = new SAMTextHeaderCodec().decode(new StringLineReader(builder.toString()), "BufferedReader");

            initialize(header);
            // Then read in the file
            final FormatUtil format = new FormatUtil();
            final List<HaplotypeMapFileEntry> entries = new ArrayList<>();
            final Map<String, HaplotypeBlock> anchorToHaplotype = new HashMap<>();

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
                final byte major = (byte) fields[3].charAt(0);
                final byte minor = (byte) fields[4].charAt(0);
                final double maf = format.parseDouble(fields[5]);
                final String anchor = fields.length > 6 ? fields[6] : null;
                final String fpPanels = fields.length > 7 ? fields[7] : null;
                List<String> panels = null;
                if (fpPanels != null) {
                    panels = new ArrayList<>();
                    for (final String panel : fpPanels.split(",")) {
                        panels.add(panel);
                    }
                }

                // If it's the anchor snp, start the haplotype
                if (anchor == null || anchor.trim().equals("") || name.equals(anchor)) {
                    final HaplotypeBlock type = new HaplotypeBlock(maf);
                    type.addSnp(new Snp(name, chrom, pos, major, minor, maf, panels));
                    anchorToHaplotype.put(name, type);
                } else {  // Otherwise save it for later
                    final HaplotypeMapFileEntry entry = makeHaplotypeMapFileEntry(
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
            anchorToHaplotype.values().forEach(this::addHaplotype);
        } catch (IOException ioe) {
            throw new PicardException("Error parsing haplotype map.", ioe);
        } finally {
            if (in != null) {
                try {
                    in.close();
                } catch (Exception e) { /* do nothing */ }
            }
        }
    }

    private HaplotypeMapFileEntry makeHaplotypeMapFileEntry(final String chrom, final int pos, final String name,
                                                            final byte major, final byte minor, final double maf,
                                                            final String anchorSnp, final List<String> fingerprintPanels) {
        return new HaplotypeMapFileEntry(chrom, pos, name, major, minor, maf, anchorSnp, fingerprintPanels);
    }

    /** Constructs an empty HaplotypeMap using the provided SAMFileHeader's sequence dictionary. */
    public HaplotypeMap(final SAMFileHeader header) {
        initialize(header);
    }

    public HaplotypeMap(final File file) {
        if (VcfUtils.isVariantFile(file)){
            fromVcf(file);
        } else {
            fromHaplotypeDatabase(file);
        }
    }

    private void initialize(final SAMFileHeader header){
        this.header = header;
        this.intervals = new IntervalList(header);
    }

    /**
     * Adds a HaplotypeBlock to the map and updates all the relevant caches/indices.
     */
    public void addHaplotype(final HaplotypeBlock haplotypeBlock) {
        this.haplotypeBlocks.add(haplotypeBlock);

        for (final Snp snp : haplotypeBlock.getSnps()) {
            if (haplotypesBySnp.containsKey(snp)) {
                throw new IllegalStateException("Same snp name cannot be used twice" + snp.toString());
            }

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
        this.intervals = this.intervals.sorted(); // TODO: should probably do this elsewhere
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

    public void writeAsVcf(final File output, final File refFile) throws FileNotFoundException {
        ReferenceSequenceFile ref = new IndexedFastaSequenceFile(refFile);
        try (VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(output)
                .setReferenceDictionary(ref.getSequenceDictionary())
                .build()) {

            final VCFHeader vcfHeader = new VCFHeader(
                    VCFUtils.withUpdatedContigsAsLines(Collections.emptySet(), refFile, header.getSequenceDictionary(), false),
                    Collections.singleton(HET_GENOTYPE_FOR_PHASING));

            VCFUtils.withUpdatedContigsAsLines(Collections.emptySet(), refFile, header.getSequenceDictionary(), false);

            vcfHeader.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(), VCFHeaderVersion.VCF4_2.getVersionString()));
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele Frequency, for each ALT allele, in the same order as listed"));
            vcfHeader.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
            vcfHeader.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.PHASE_SET_KEY, 1, VCFHeaderLineType.String, "Phase-set identifier for phased genotypes."));
            vcfHeader.addMetaDataLine(new VCFHeaderLine(VCFHeader.SOURCE_KEY,"HaplotypeMap::writeAsVcf"));
            vcfHeader.addMetaDataLine(new VCFHeaderLine("reference","HaplotypeMap::writeAsVcf"));


          //  vcfHeader.addMetaDataLine(new VCFHeaderLine());
            writer.writeHeader(vcfHeader);
            final LinkedList<VariantContext> variants = new LinkedList<>(this.asVcf(ref));
            variants.sort(vcfHeader.getVCFRecordComparator());
            variants.forEach(writer::add);
        }
    }

    public Collection<VariantContext> asVcf(final ReferenceSequenceFile ref) {

        final List<VariantContext> entries = new ArrayList<>();
        final SortedSet<Snp> snps = new TreeSet<>(getAllSnps());
        final Map<Snp, Boolean> allele1MatchesReference = new HashMap<>(snps.size());

        for( final Snp snp : snps) {
            final ReferenceSequence seq = ref.getSubsequenceAt(snp.getChrom(), snp.getPos(), snp.getPos());
            if (seq.getBases()[0] == snp.getAllele1()) {
                allele1MatchesReference.put(snp, true);
            } else if (seq.getBases()[0] == snp.getAllele2()) {
                allele1MatchesReference.put(snp, false);
            } else {
                throw new RuntimeException("One of the two alleles should agree with the reference: " + snp.toString());
            }
        }

        for (final HaplotypeBlock block : this.getHaplotypes()) {
            Snp anchorSnp = null;
            final SortedSet<Snp> blocksSnps = new TreeSet<>(block.getSnps());

            for (final Snp snp : blocksSnps) {

                if (anchorSnp == null) {
                    anchorSnp = snp;
                }

                final String alleleString = snp.getAlleleString();
                final boolean swap = allele1MatchesReference.get(snp);
                final String reference = !swap ? alleleString.substring(0, 1) : alleleString.substring(1, 2);
                final String alternate = swap ? alleleString.substring(0, 1) : alleleString.substring(1, 2);

                final double maf = !swap ? snp.getMaf() : 1 - snp.getMaf();

                VariantContextBuilder builder = new VariantContextBuilder()
                        .chr(snp.getChrom())
                        .start(snp.getPos())
                        .stop(snp.getPos())
                        .alleles(reference, alternate)
                        .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, maf)
                        .id(snp.getName());
                GenotypeBuilder genotypeBuilder = new GenotypeBuilder(HET_GENOTYPE_FOR_PHASING);

                if (blocksSnps.size() > 1 && swap) {
                    genotypeBuilder.alleles(Arrays.asList(builder.getAlleles().get(1), builder.getAlleles().get(0)));
                } else {
                    genotypeBuilder.alleles(builder.getAlleles());
                }

                if (blocksSnps.size() > 1) {
                    genotypeBuilder.phased(true);
                    genotypeBuilder.attribute(VCFConstants.PHASE_SET_KEY, anchorSnp.getPos());
                }
                builder.genotypes(genotypeBuilder.make());

                entries.add(builder.make());
            }
        }
        return entries;
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

            final List<HaplotypeMapFileEntry> entries = new ArrayList<>();
            for (final HaplotypeBlock block : this.getHaplotypes()) {
                String anchor = null;
                final SortedSet<Snp> snps = new TreeSet<>(block.getSnps());

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
            throw new PicardException("Error writing out haplotype map to file: " + file.getAbsolutePath(), ioe);
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
            this.panels = new ArrayList<>();
            if (fingerprintPanels != null) {
                this.panels.addAll(fingerprintPanels);
                Collections.sort(this.panels);
            }            
        }

        public String getPanels() {
            if (panels == null) return "";
            final StringBuilder sb = new StringBuilder();

            for (final String panel : panels) {
                if (sb.length() > 0) sb.append(',');
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
