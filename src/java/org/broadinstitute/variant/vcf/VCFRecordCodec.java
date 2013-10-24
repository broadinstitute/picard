package org.broadinstitute.variant.vcf;

import net.sf.samtools.util.SortingCollection;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * Writes VariantContext instances to an OutputStream without headers or metadata. For use
 * with SortingCollection ONLY.
 */
public class VCFRecordCodec implements SortingCollection.Codec<VariantContext> {

	final VCFCodec vcfDecoder = new VCFCodec();

	final VCFEncoder vcfEncoder;

	private PrintStream outputStream = null;

	private BufferedReader inputReader = null;

	public VCFRecordCodec(final VCFHeader header) {
		this.vcfEncoder = new VCFEncoder(header, false);
		// Explicitly set the version because it's not available in the header itself.
		this.vcfDecoder.setVCFHeader(header, VCFHeaderVersion.VCF4_1);
	}

	@Override
	public void setOutputStream(final OutputStream stream) {
		this.outputStream = new PrintStream(stream);
	}

	@Override
	public void setInputStream(final InputStream stream) {
		this.inputReader = new BufferedReader(new InputStreamReader(stream));
	}

	@Override
	public void encode(final VariantContext context) {
		this.outputStream.println(this.vcfEncoder.encode(context));
	}

	@Override
	public VariantContext decode() {
		try {
			final String line;
			return ((line = inputReader.readLine()) != null) ? this.vcfDecoder.decode(line) : null;
		} catch (final IOException ioe) {
			throw new RuntimeException("Could not decode/read a VCF record for a sorting collection: " + ioe.getMessage(), ioe);
		}
	}

	@Override
	public VCFRecordCodec clone() {
		return new VCFRecordCodec(this.vcfEncoder.getVCFHeader());
	}
}

