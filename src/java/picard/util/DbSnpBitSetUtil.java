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
package picard.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.lang.IllegalArgumentException;import java.lang.String;import java.util.BitSet;import java.util.Collection;import java.util.EnumSet;import java.util.HashMap;import java.util.Map;import java.util.Set;

/**
 * Utility class to use with DbSnp files to determine is a locus is
 * a dbSnp site.
 */
public class DbSnpBitSetUtil {

    private final Map<String, BitSet> sequenceToBitSet = new HashMap<String,BitSet>();

    /** Little tuple class to contain one bitset for SNPs and another for Indels. */
    public static class DbSnpBitSets {
        public DbSnpBitSetUtil snps;
        public DbSnpBitSetUtil indels;
    }

    /** Private empty contructor for use by factory methods only. */
    private DbSnpBitSetUtil() { }

    /** Constructor that creates a bit set with bits set to true for all variant types. */
    public DbSnpBitSetUtil(final File dbSnpFile, final SAMSequenceDictionary sequenceDictionary) {
        this(dbSnpFile, sequenceDictionary, EnumSet.noneOf(VariantType.class));
    }

    /**
     * Constructor.
     *
     * For each sequence, creates  a BitSet that denotes whether a dbSNP entry
     * is present at each base in the reference sequence.  The set is
     * reference.length() + 1 so that it can be indexed by 1-based reference base.
     * True means dbSNP present, false means no dbSNP present.
     *
     * @param dbSnpFile in VCF format.
     * @param sequenceDictionary Optionally, a sequence dictionary corresponding to the dbSnp file, else null.
     * If present, BitSets will be allocated more efficiently because the maximum size will be known.
     * @param variantsToMatch what types of variants to load.
     */
    public DbSnpBitSetUtil(final File dbSnpFile,
                           final SAMSequenceDictionary sequenceDictionary,
                           final Collection<VariantType> variantsToMatch) {

        if (dbSnpFile == null) throw new IllegalArgumentException("null dbSnpFile");
        final Map<DbSnpBitSetUtil, Set<VariantType>> tmp = new HashMap<DbSnpBitSetUtil, Set<VariantType>>();
        tmp.put(this, EnumSet.copyOf(variantsToMatch));
        loadVcf(dbSnpFile, sequenceDictionary, tmp);
    }

    /** Factory method to create both a SNP bitmask and an indel bitmask in a single pass of the VCF. */
    public static DbSnpBitSets createSnpAndIndelBitSets(final File dbSnpFile,
                                                        final SAMSequenceDictionary sequenceDictionary) {

        final DbSnpBitSets sets = new DbSnpBitSets();
        sets.snps   = new DbSnpBitSetUtil();
        sets.indels = new DbSnpBitSetUtil();

        final Map<DbSnpBitSetUtil, Set<VariantType>> map = new HashMap<DbSnpBitSetUtil, Set<VariantType>>();
        map.put(sets.snps,   EnumSet.of(VariantType.SNP));
        map.put(sets.indels, EnumSet.of(VariantType.insertion, VariantType.deletion));
        loadVcf(dbSnpFile, sequenceDictionary, map);
        return sets;
    }

    /** Private helper method to read through the VCF and create one or more bit sets. */
    private static void loadVcf(final File dbSnpFile,
                                final SAMSequenceDictionary sequenceDictionary,
                                final Map<DbSnpBitSetUtil, Set<VariantType>> bitSetsToVariantTypes) {

        final VCFFileReader variantReader = new VCFFileReader(dbSnpFile);
        final CloseableIterator<VariantContext> variantIterator = variantReader.iterator();

        while (variantIterator.hasNext()) {
            final VariantContext kv = variantIterator.next();

            for (final Map.Entry<DbSnpBitSetUtil, Set<VariantType>> tuple : bitSetsToVariantTypes.entrySet()) {
                final DbSnpBitSetUtil bitset            = tuple.getKey();
                final Set<VariantType> variantsToMatch  = tuple.getValue();

                BitSet bits = bitset.sequenceToBitSet.get(kv.getChr());
                if (bits == null) {
                    final int nBits;
                    if (sequenceDictionary == null) nBits = kv.getEnd() + 1;
                    else nBits = sequenceDictionary.getSequence(kv.getChr()).getSequenceLength() + 1;
                    bits = new BitSet(nBits);
                    bitset.sequenceToBitSet.put(kv.getChr(), bits);
                }
                if (variantsToMatch.isEmpty() ||
                        (kv.isSNP() && variantsToMatch.contains(VariantType.SNP)) ||
                        (kv.isIndel() && variantsToMatch.contains(VariantType.insertion)) ||
                        (kv.isIndel() && variantsToMatch.contains(VariantType.deletion))) {

                    for (int i = kv.getStart(); i <= kv.getEnd(); i++)  bits.set(i, true);
                }
            }
        }

        CloserUtil.close(variantIterator);
        CloserUtil.close(variantReader);
    }

    /**
     * Returns true if there is a dbSnp entry at pos in sequenceName, otherwise false
     */
    public boolean isDbSnpSite(final String sequenceName, final int pos) {
        // When we have a dbSnpFile with no sequence dictionary, this line will be necessary
        if (sequenceToBitSet.get(sequenceName) == null) {
            return false;
        }
        if (pos > sequenceToBitSet.get(sequenceName).length()) {
            return false;
        }
        return sequenceToBitSet.get(sequenceName).get(pos);
    }

}
