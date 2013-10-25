package org.broadinstitute.variant.vcf;

import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.LazyGenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.IntGenotypeFieldAccessors;

import java.lang.reflect.Array;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Functions specific to encoding VCF records.
 */
public class VCFEncoder {

    /**
     * The encoding used for VCF files: ISO-8859-1
     */
    public static final Charset VCF_CHARSET = Charset.forName("ISO-8859-1");
    private static final String QUAL_FORMAT_STRING = "%.2f";
	private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

	private final IntGenotypeFieldAccessors GENOTYPE_FIELD_ACCESSORS = new IntGenotypeFieldAccessors();

	private VCFHeader header;

	private boolean allowMissingFieldsInHeader = false;

	/**
	 * Prepare a VCFEncoder that will encode records appropriate to the given VCF header, optionally
	 * allowing missing fields in the header.
	 */
	public VCFEncoder(final VCFHeader header, final boolean allowMissingFieldsInHeader) {
		if (header == null) throw new NullPointerException("The VCF header must not be null.");
		this.header = header;
		this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
	}

	/**
	 * Please see the notes in the default constructor
	 */
	@Deprecated
	public void setVCFHeader(final VCFHeader header) {
		this.header = header;
	}

	/**
	 * Please see the notes in the default constructor
	 */
	@Deprecated
	public void setAllowMissingFieldsInHeader(final boolean allow) {
		this.allowMissingFieldsInHeader = allow;
	}

	public String encode(final VariantContext context) {
		if (this.header == null) {
			throw new NullPointerException("The header field must be set on the VCFEncoder before encoding records.");
		}

		final StringBuilder stringBuilder = new StringBuilder();

		// CHROM
		stringBuilder.append(context.getChr()).append(VCFConstants.FIELD_SEPARATOR);

		// POS
		stringBuilder.append(String.valueOf(context.getStart())).append(VCFConstants.FIELD_SEPARATOR);

		// ID
		stringBuilder.append(context.getID()).append(VCFConstants.FIELD_SEPARATOR);

		// REF
		stringBuilder.append(context.getReference().getDisplayString()).append(VCFConstants.FIELD_SEPARATOR);

		// ALT
		if ( context.isVariant() ) {
			Allele altAllele = context.getAlternateAllele(0);
			String alt = altAllele.getDisplayString();
			stringBuilder.append(alt);

			for (int i = 1; i < context.getAlternateAlleles().size(); i++) {
				altAllele = context.getAlternateAllele(i);
				alt = altAllele.getDisplayString();
				stringBuilder.append(",");
				stringBuilder.append(alt);
			}
		} else {
			stringBuilder.append(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
		}

		stringBuilder.append(VCFConstants.FIELD_SEPARATOR);

		// QUAL
		if ( ! context.hasLog10PError()) stringBuilder.append(VCFConstants.MISSING_VALUE_v4);
		else stringBuilder.append(formatQualValue(context.getPhredScaledQual()));
		stringBuilder.append(VCFConstants.FIELD_SEPARATOR);

		// FILTER
		stringBuilder.append(getFilterString(context)).append(VCFConstants.FIELD_SEPARATOR);

		// INFO
		final Map<String, String> infoFields = new TreeMap<String, String>();
		for (final Map.Entry<String, Object> field : context.getAttributes().entrySet() ) {
			if ( ! this.header.hasInfoLine(field.getKey())) fieldIsMissingFromHeaderError(context, field.getKey(), "INFO");

			final String outputValue = formatVCFField(field.getValue());
			if (outputValue != null) infoFields.put(field.getKey(), outputValue);
		}
		writeInfoString(infoFields, stringBuilder);

		// FORMAT
		final GenotypesContext gc = context.getGenotypes();
		if (gc.isLazyWithData() && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String) {
			stringBuilder.append(VCFConstants.FIELD_SEPARATOR);
			stringBuilder.append(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
		} else {
			final List<String> genotypeAttributeKeys = context.calcVCFGenotypeKeys(this.header);
			if ( ! genotypeAttributeKeys.isEmpty()) {
				for (final String format : genotypeAttributeKeys)
					if ( ! this.header.hasFormatLine(format))
						fieldIsMissingFromHeaderError(context, format, "FORMAT");

				final String genotypeFormatString = ParsingUtils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);

				stringBuilder.append(VCFConstants.FIELD_SEPARATOR);
				stringBuilder.append(genotypeFormatString);

				final Map<Allele, String> alleleStrings = buildAlleleStrings(context);
				addGenotypeData(context, alleleStrings, genotypeAttributeKeys, stringBuilder);
			}
		}

		return stringBuilder.toString();
	}

	VCFHeader getVCFHeader() {
		return this.header;
	}

	boolean getAllowMissingFieldsInHeader() {
		return this.allowMissingFieldsInHeader;
	}

	private String getFilterString(final VariantContext vc) {
		if (vc.isFiltered()) {
			for (final String filter : vc.getFilters()) {
				if ( ! this.header.hasFilterLine(filter)) fieldIsMissingFromHeaderError(vc, filter, "FILTER");
			}

			return ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters()));
		}
		else if (vc.filtersWereApplied()) return VCFConstants.PASSES_FILTERS_v4;
		else return VCFConstants.UNFILTERED;
	}

	private String formatQualValue(final double qual) {
		String s = String.format(QUAL_FORMAT_STRING, qual);
		if ( s.endsWith(QUAL_FORMAT_EXTENSION_TO_TRIM) )
			s = s.substring(0, s.length() - QUAL_FORMAT_EXTENSION_TO_TRIM.length());
		return s;
	}

	private void fieldIsMissingFromHeaderError(final VariantContext vc, final String id, final String field) {
		if ( ! allowMissingFieldsInHeader)
			throw new IllegalStateException("Key " + id + " found in VariantContext field " + field
					+ " at " + vc.getChr() + ":" + vc.getStart()
					+ " but this key isn't defined in the VCFHeader.  We require all VCFs to have"
					+ " complete VCF headers by default.");
	}

	String formatVCFField(final Object val) {
		final String result;
		if ( val == null )
			result = VCFConstants.MISSING_VALUE_v4;
		else if ( val instanceof Double )
			result = formatVCFDouble((Double) val);
		else if ( val instanceof Boolean )
			result = (Boolean)val ? "" : null; // empty string for true, null for false
		else if ( val instanceof List ) {
			result = formatVCFField(((List)val).toArray());
		} else if ( val.getClass().isArray() ) {
			final int length = Array.getLength(val);
			if ( length == 0 )
				return formatVCFField(null);
			final StringBuilder sb = new StringBuilder(formatVCFField(Array.get(val, 0)));
			for ( int i = 1; i < length; i++) {
				sb.append(",");
				sb.append(formatVCFField(Array.get(val, i)));
			}
			result = sb.toString();
		} else
			result = val.toString();

		return result;
	}

	/**
	 * Takes a double value and pretty prints it to a String for display
	 *
	 * Large doubles => gets %.2f style formatting
	 * Doubles < 1 / 10 but > 1/100 </>=> get %.3f style formatting
	 * Double < 1/100 => %.3e formatting
	 * @param d
	 * @return
	 */
	public static String formatVCFDouble(final double d) {
		final String format;
		if ( d < 1 ) {
			if ( d < 0.01 ) {
				if ( Math.abs(d) >= 1e-20 )
					format = "%.3e";
				else {
					// return a zero format
					return "0.00";
				}
			} else {
				format = "%.3f";
			}
		} else {
			format = "%.2f";
		}

		return String.format(format, d);
	}

	static int countOccurrences(final char c, final String s) {
		int count = 0;
		for (int i = 0; i < s.length(); i++) {
			count += s.charAt(i) == c ? 1 : 0;
		}
		return count;
	}

	static boolean isMissingValue(final String s) {
		// we need to deal with the case that it's a list of missing values
		return (countOccurrences(VCFConstants.MISSING_VALUE_v4.charAt(0), s) + countOccurrences(',', s) == s.length());
	}

	/*
	 * Add the genotype data
	 */
	private void addGenotypeData(final VariantContext vc, final Map<Allele, String> alleleMap, final List<String> genotypeFormatKeys, final StringBuilder builder) {
		final int ploidy = vc.getMaxPloidy(2);

		for (final String sample : this.header.getGenotypeSamples()) {
			builder.append(VCFConstants.FIELD_SEPARATOR);

			Genotype g = vc.getGenotype(sample);
			if (g == null) g = GenotypeBuilder.createMissing(sample, ploidy);

			final List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
			for (final String field : genotypeFormatKeys) {
				if (field.equals(VCFConstants.GENOTYPE_KEY)) {
					if ( ! g.isAvailable()) {
						throw new IllegalStateException("GTs cannot be missing for some samples if they are available for others in the record");
					}

					writeAllele(g.getAllele(0), alleleMap, builder);
					for (int i = 1; i < g.getPloidy(); i++) {
						builder.append(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
						writeAllele(g.getAllele(i), alleleMap, builder);
					}
					continue;

				} else {
					final String outputValue;
					if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY ) ) {
						outputValue = g.isFiltered() ? g.getFilters() : VCFConstants.PASSES_FILTERS_v4;
					} else {
						final IntGenotypeFieldAccessors.Accessor accessor = GENOTYPE_FIELD_ACCESSORS.getAccessor(field);
						if ( accessor != null ) {
							final int[] intValues = accessor.getValues(g);
							if ( intValues == null )
								outputValue = VCFConstants.MISSING_VALUE_v4;
							else if ( intValues.length == 1 ) // fast path
								outputValue = Integer.toString(intValues[0]);
							else {
								final StringBuilder sb = new StringBuilder();
								sb.append(intValues[0]);
								for ( int i = 1; i < intValues.length; i++) {
									sb.append(",");
									sb.append(intValues[i]);
								}
								outputValue = sb.toString();
							}
						} else {
							Object val = g.hasExtendedAttribute(field) ? g.getExtendedAttribute(field) : VCFConstants.MISSING_VALUE_v4;

							final VCFFormatHeaderLine metaData = this.header.getFormatHeaderLine(field);
							if ( metaData != null ) {
								final int numInFormatField = metaData.getCount(vc);
								if ( numInFormatField > 1 && val.equals(VCFConstants.MISSING_VALUE_v4) ) {
									// If we have a missing field but multiple values are expected, we need to construct a new string with all fields.
									// For example, if Number=2, the string has to be ".,."
									final StringBuilder sb = new StringBuilder(VCFConstants.MISSING_VALUE_v4);
									for ( int i = 1; i < numInFormatField; i++ ) {
										sb.append(",");
										sb.append(VCFConstants.MISSING_VALUE_v4);
									}
									val = sb.toString();
								}
							}

							// assume that if key is absent, then the given string encoding suffices
							outputValue = formatVCFField(val);
						}
					}

					if ( outputValue != null )
						attrs.add(outputValue);
				}
			}

			// strip off trailing missing values
			for (int i = attrs.size()-1; i >= 0; i--) {
				if ( isMissingValue(attrs.get(i))) attrs.remove(i);
				else break;
			}

			for (int i = 0; i < attrs.size(); i++) {
				if ( i > 0 || genotypeFormatKeys.contains(VCFConstants.GENOTYPE_KEY)) {
					builder.append(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
				}
				builder.append(attrs.get(i));
			}
		}
	}

	/*
	 * Create the info string; assumes that no values are null
	 */
	private void writeInfoString(final Map<String, String> infoFields, final StringBuilder builder) {
		if ( infoFields.isEmpty() ) {
			builder.append(VCFConstants.EMPTY_INFO_FIELD);
			return;
		}

		boolean isFirst = true;
		for (final Map.Entry<String, String> entry : infoFields.entrySet()) {
			if (isFirst) isFirst = false;
			else builder.append(VCFConstants.INFO_FIELD_SEPARATOR);

			builder.append(entry.getKey());

			if ( ! entry.getValue().equals("")) {
				final VCFInfoHeaderLine metaData = this.header.getInfoHeaderLine(entry.getKey());
				if ( metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0 ) {
					builder.append("=");
					builder.append(entry.getValue());
				}
			}
		}
	}

	private Map<Allele, String> buildAlleleStrings(final VariantContext vc) {
		final Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size()+1);
		alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

		final List<Allele> alleles = vc.getAlleles();
		for ( int i = 0; i < alleles.size(); i++ ) {
			alleleMap.put(alleles.get(i), String.valueOf(i));
		}

		return alleleMap;
	}

	private void writeAllele(final Allele allele, final Map<Allele, String> alleleMap, final StringBuilder builder) {
		final String encoding = alleleMap.get(allele);
		if ( encoding == null )
			throw new RuntimeException("Allele " + allele + " is not an allele in the variant context");
		builder.append(encoding);
	}
}
