/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.pedigree.PedFile;
import picard.pedigree.Sex;

import java.io.File;
import java.util.Map;
import java.util.Set;

/**
 * Attempts to use a set of genome data samples to determine sex within the cohort provided.
 *
 * @author Tim Fennell
 * @author Jonathan Barlev
 * @author Yossi Farjoun
 */

// Abstract class that provides the outline of a sex inferencing CLP.
// The implementing class will have to implement getSexInferencer()
public abstract class InferSex extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output pedigree file name.")
    public File OUTPUT;

    @Argument(doc="List of possible names for male sex chromosome(s)")
    public Set<String> MALE_CHROMS = CollectionUtil.makeSet("Y", "chrY");

    @Argument(doc="List of possible names for female sex chromosome(s)")
    public Set<String> FEMALE_CHROMS = CollectionUtil.makeSet("X", "chrX");

    final Log log = Log.getInstance(InferSex.class);

    abstract SexInferenceEngine getSexInference();

    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        // Determine gender for everyone
        log.info("Determining sample sexes");
        final Map<String,Sex> sampleSexes =  getSexInference().determineSexes();
        PedFile.fromSexMap(sampleSexes).write(OUTPUT);

        return 0;
    }
}
