package net.sf.picard.vcf;

import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.RuntimeIOException;
import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;

public class VariantContextIterator implements CloseableIterator<VariantContext> {

	private final VCFCodec vcfCodec = new VCFCodec();

	private final VCFHeader vcfHeader;

	private final AsciiLineReader reader;

	private String line = null;

	public VariantContextIterator(final File vcfFile) {
		this.reader = new AsciiLineReader(IoUtil.openFileForReading(vcfFile));
		final Object header = vcfCodec.readHeader(reader);
		if ( ! (header instanceof VCFHeader)) {
			throw new IllegalArgumentException("The file " + vcfFile.getAbsolutePath() + " did not have a VCF header");
		}
		this.vcfHeader = (VCFHeader) header;
	}

	// TODO: Add a c'tor that reads intervals.

	@Override
	public void close() {
		CloserUtil.close(reader);
	}

	public VCFHeader getHeader() {
		return this.vcfHeader;
	}

	@Override
	public boolean hasNext() {
		try {
			if (line == null) line = reader.readLine();
		} catch (IOException e) {
			throw new RuntimeIOException(e);
		}
		return line != null;
	}

	@Override
	public VariantContext next() {
		if ( ! this.hasNext()) throw new NoSuchElementException("Called next() on an exhausted VariantContextIterator");
		final String tmp = line;
		line = null;
		return vcfCodec.decode(tmp);
	}

	/**
	 * Unsupported.
	 */
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
