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

import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.PicardException;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * class to represent a genetic fingerprint as a set of HaplotypeProbabilities
 * objects that give the relative probabilities of each of the possible haplotypes
 * at a locus.
 *
 * @author Tim Fennell
 */
public class Fingerprint extends TreeMap<HaplotypeBlock, HaplotypeProbabilities> {
    private final String sample;
    private final Path source;
    private final String info;

    public Fingerprint(final String sample, final Path source, final String info) {
        this.sample = sample;
        this.source = source;
        this.info = info;
    }

    public String getSample() { return sample; }

    public Path getSource() { return source; }

    public String getInfo() { return info; }

    public String getPrintableId() {
        return getSample() + "@" + (source == null ? "" : source.toUri().toString()) + (info == null ? "" : (":" + info));
    }

    public void add(final HaplotypeProbabilities h) {
        put(h.getHaplotype(), h);
    }

    /**
     * Merges the likelihoods from the supplied Fingerprint into the likelihoods for this fingerprint.
     */
    public Fingerprint merge(final Fingerprint other) {
        final Set<HaplotypeBlock> haps = new HashSet<>();
        haps.addAll(keySet());
        haps.addAll(other.keySet());

        for (final HaplotypeBlock haplotype : haps) {
            HaplotypeProbabilities probabilities = get(haplotype);
            final HaplotypeProbabilities otherProbabilities = other.get(haplotype);
            if (probabilities == null) {
                probabilities = otherProbabilities.deepCopy();
                put(haplotype, probabilities);
            } else if (otherProbabilities != null) {
                probabilities.merge(otherProbabilities);
            }
        }
        return this;
    }

    public static Function<FingerprintIdDetails, String> getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType CROSSCHECK_BY) {
        final Function<FingerprintIdDetails, String> groupByTemp;
        switch (CROSSCHECK_BY) {
            case READGROUP:
                groupByTemp = details -> details.platformUnit;
                break;
            case LIBRARY:
                groupByTemp = details -> details.sample + "::" + details.library;
                break;
            case FILE:
                groupByTemp = details -> details.file + "::" + details.sample;
                break;
            case SAMPLE:
                groupByTemp = details -> details.sample;
                break;
            default:
                throw new PicardException("unpossible");
        }

        // if the groupBy string is null (e.g. a vcf file has no read group info) then the hashcode is
        // used intending to be unique per object (ignoring possible collisions)
        return key -> {
            final String temp = groupByTemp.apply(key);
            return temp == null ? Integer.toString(key.hashCode()) : temp;
        };
    }

    public static Map<FingerprintIdDetails, Fingerprint> mergeFingerprintsBy(
            final Map<FingerprintIdDetails, Fingerprint> fingerprints,
            final Function<FingerprintIdDetails, String> by) {

        // collect the various entries according to the grouping "by"

        final Map<String, List<Map.Entry<FingerprintIdDetails, Fingerprint>>> collection =
                fingerprints.entrySet()
                        .stream()
                        .collect(Collectors.groupingBy(entry -> by.apply(entry.getKey())));

        return collection.entrySet().stream()
                .collect(Collectors.toMap(
                        entry -> {
                            // merge the keys (unequal values are eliminated by merge).
                            final List<Map.Entry<FingerprintIdDetails, Fingerprint>> entryList = entry.getValue();
                            final FingerprintIdDetails finalId;
                            if (entryList.size() == 1) {
                                finalId = entryList.get(0).getKey();
                            } else {
                                finalId = new FingerprintIdDetails();
                                entryList.forEach(id -> finalId.merge(id.getKey()));
                            }
                            finalId.group = entry.getKey();
                            return finalId;

                        }, entry -> {
                            // merge the values by merging the fingerprints.

                            //use the "by" function to determine the "info" part of the fingerprint
                            final Set<Fingerprint> fingerprintsSet = entry.getValue().stream().map(Map.Entry::getValue).collect(Collectors.toSet());
                            if (fingerprintsSet.size() == 1) {
                                return fingerprintsSet.iterator().next();
                            }
                            final FingerprintIdDetails firstDetail = entry.getValue().get(0).getKey();
                            final Fingerprint mergedFp = new Fingerprint(firstDetail.sample, null, by.apply(firstDetail));

                            fingerprintsSet.forEach(mergedFp::merge);
                            return mergedFp;
                        }));
    }

    enum CrosscheckMode implements CommandLineParser.ClpEnum {
        CHECK_SAME_SAMPLE {
            @Override
            public String getHelpDoc() {
                return "In this mode, each sample in INPUT will only be checked against a single corresponding sample in SECOND_INPUT. " +
                        "If a corresponding sample cannot be found, the program will proceed, but report the missing samples" +
                        " and return the value specified in EXIT_CODE_WHEN_MISMATCH. The corresponding samples are those that equal each other, after possible renaming " +
                        "via INPUT_SAMPLE_MAP and SECOND_INPUT_SAMPLE_MAP. In this mode CROSSCHECK_BY must be SAMPLE.";
            }
        },
        CHECK_ALL_OTHERS {
            @Override
            public String getHelpDoc() {
                return "In this mode, each sample in INPUT will be checked against all the samples in SECOND_INPUT.";
            }
        }
    }
}
