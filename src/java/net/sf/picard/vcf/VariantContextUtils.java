package net.sf.picard.vcf;

import net.sf.picard.PicardException;
import net.sf.samtools.Defaults;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.EnumSet;
import java.util.zip.GZIPOutputStream;

public final class VariantContextUtils {

	/**
	 * Create a VariantContextWriter for the given output file. If the output file has a .gz extension
	 * a GZIPOutputStream is used to compress the data on the fly, otherwise an "ordinary", non-
	 * compressing VariantContextWriter is returned. If compressed, the returned ...Writer will create
	 * .gz a file that conforms to (and can be decompressed by) ordinary gzip. No sequence dictionaries
	 * are used.
	 *
	 * Default compression level for compressed files is 5: it seems to be a good tradeoff between
	 * compression ratio and time.
	 */
	public static VariantContextWriter getConditionallyCompressingWriter(final File output, final EnumSet<Options> options) {
		return output.getName().endsWith(".gz")
				? getCompressingWriter(output, options)
				: VariantContextWriterFactory.create(output, null, options);
	}

	/**
	 * Create a compressing VariantContextWriter for the given File, even if the extension on the File
	 * is not .gz. The returned ...Writer will create a file that conforms to (and can be decompressed by)
	 * ordinary gzip. No sequence dictionaries are used.
	 *
	 * Default compression level for compressed files is 5: it seems to be a good tradeoff between
	 * compression ratio and time.
	 */
	public static VariantContextWriter getCompressingWriter(final File output, final EnumSet<Options> options) {
		try {
			final GZIPOutputStream gzipOutputStream = new GZIPOutputStream(new FileOutputStream(output)) {{
				def.setLevel(Defaults.COMPRESSION_LEVEL);
			}};
			final OutputStream outputStream = new BufferedOutputStream(gzipOutputStream);
			return VariantContextWriterFactory.create(output, outputStream, null, options);

		} catch (final Exception e) {
			throw new PicardException("Could not create a compressed output stream for the VCF writer: " + e.getMessage(), e);
		}
	}
}
