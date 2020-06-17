/*
 * The MIT License
 *
 * Copyright (c) 2020 The Broad Institute
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
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;

/**
 * Liftover SNPs in HaplotypeMaps from one reference to another
 *
 */
@CommandLineProgramProperties(
        summary = "Lifts over a haplotype database from one reference to another. Based on UCSC liftOver.\n" +
                " Uses a UCSC chain file to guide the liftOver.",
        oneLineSummary = "Lifts over a haplotype database from one reference to another",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class LiftOverHaplotypeMap extends CommandLineProgram {

    @Argument(doc = "Haplotype database to be lifted over.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Where to write lifted-over haplotype database.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "Sequence dictionary to write into the output haplotype database. (Any file from which a dictionary is extractable.)",
            shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME)
    public File SEQUENCE_DICTIONARY;

    @Argument(doc = "Chain file that guides LiftOver. (UCSC format)")
    public File CHAIN;

    public static final int LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS = 101;

    private static final Log log = Log.getInstance(LiftOverHaplotypeMap.class);

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsReadable(CHAIN);
        IOUtil.assertFileIsWritable(OUTPUT);

        final LiftOver liftOver = new LiftOver(CHAIN);

        final HaplotypeMap fromMap = new HaplotypeMap(INPUT);

        final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY);

        liftOver.validateToSequences(dict);
        final SAMFileHeader newHeader = new SAMFileHeader(dict);
        final HaplotypeMap toMap = new HaplotypeMap(newHeader);

        boolean anyFailed = false;

        for (final HaplotypeBlock fromHaplotypeBlock : fromMap.getHaplotypes()) {
            final HaplotypeBlock toHaplotypeBlock = new HaplotypeBlock(fromHaplotypeBlock.getMaf());

            for (final Snp snp : fromHaplotypeBlock.getSnps()) {
                final Interval interval = liftOver.liftOver(
                        new Interval(snp.getChrom(), snp.getPos(), snp.getPos()));
                if (interval != null) {
                    toHaplotypeBlock.addSnp(new Snp(snp.getName(), interval.getContig(), interval.getStart(),
                            snp.getAllele1(), snp.getAllele2(), snp.getMaf(),
                            snp.getFingerprintPanels()));
                } else {
                    anyFailed = true;
                    log.warn("Liftover failed for ", snp.getName(), "(", snp.getChrom(), ":", snp.getPos(), ")");
                }
            }

            toMap.addHaplotype(toHaplotypeBlock);
        }

        toMap.writeToFile(OUTPUT);
        return anyFailed ? LIFTOVER_FAILED_FOR_ONE_OR_MORE_SNPS : 0;
    }
}
