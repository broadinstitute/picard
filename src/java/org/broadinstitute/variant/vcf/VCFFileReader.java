package org.broadinstitute.variant.vcf;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

/**
 * Simplified interface for reading from VCF/BCF files.
 */
public class VCFFileReader implements Closeable, Iterable<VariantContext> {

	private final FeatureReader<VariantContext> reader;

	/**
	 * Returns true if the given file appears to be a BCF file.
	 */
	public static boolean isBCF(final File file) {
		return file.getAbsolutePath().endsWith(".bcf");
	}

	/**
	 * Returns the SAMSequenceDictionary from the provided VCF file.
	 */
	public static SAMSequenceDictionary getSequenceDictionary(final File file) {
		final SAMSequenceDictionary dict = new VCFFileReader(file).getFileHeader().getSequenceDictionary();
		CloserUtil.close(file);
		return dict;
	}

    /** Constructs a VCFFileReader that requires the index to be present. */
	public VCFFileReader(final File file) {
		this(file, true);
	}

    /** Allows construction of a VCFFileReader that will or will not assert the presence of an index as desired. */
	public VCFFileReader(final File file, final boolean requireIndex) {
		this.reader = AbstractFeatureReader.getFeatureReader(
						file.getAbsolutePath(),
						isBCF(file) ? new BCF2Codec() : new VCFCodec(),
						requireIndex);
	}

    /** Returns the VCFHeader associated with this VCF/BCF file. */
	public VCFHeader getFileHeader() {
		return (VCFHeader) reader.getHeader();
	}

    /** Returns an iterator over all records in this VCF/BCF file. */
	public CloseableIterator<VariantContext> iterator() {
		try { return reader.iterator(); }
        catch (final IOException ioe) {
			throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
		}
	}

    /** Queries for records within the region specified. */
    public CloseableIterator<VariantContext> query(final String chrom, final int start, final int end) {
        try { return reader.query(chrom, start, end); }
        catch (final IOException ioe) {
            throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
        }
    }

	public void close() {
		try { this.reader.close(); }
        catch (final IOException ioe) {
			throw new TribbleException("Could not close a variant context feature reader.", ioe);
		}
	}
}
