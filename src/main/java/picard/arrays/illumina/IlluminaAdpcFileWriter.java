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

package picard.arrays.illumina;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * A class to encompass writing an Illumina adpc.bin file.
 *
 * <p> This file is used as input to verifyIDIntensity, a contamination checking tool
 * for Illumina Genotyping Arrays.</p>
 *
 * <p>Here is the format of the file</p>
 *
 * <p>The file size is (16 byte offset) + (18 bytes) * # INDS * SNP</p>
 *
 * <p>Note that I do not know what the header (16 bytes) should contain.
 * verifyIDIntensity, does not care, so we are putting garbage in there now.
 * I presume it *should* contain the number of probes (at a minimum)</p>
 *
 * <p>Each genotype is ordered in the following way.
 *           (ind1-snp1) - (ind1-snp2) - (ind1-snp3) ... (ind1-snpN) (ind2-snp1)</p>
 *
 * <p>The 18 bytes are composed of the following information.</p>
 * <pre>
 *   2-short - A intensity
 *   2-short - B intensity
 *   4-float - A normalized intensity
 *   4-float - B normalized intensity
 *   4-float - GC score : clustering confidence
 *   2-short - genotype value : 0 (AA) 1 (AB) 2 (BB) 3 (NN)
 * </pre>
 */

public class IlluminaAdpcFileWriter implements AutoCloseable {
    private final DataOutputStream outputStream;

    public IlluminaAdpcFileWriter(final File adpcFile) throws IOException {
        outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(adpcFile)));
        writeHeaderData();
    }

    private void writeHeaderData() throws IOException {
        outputStream.write("1234567890123456".getBytes());
    }

    public void write(Iterable<Record> illuminaAdpcRecords) throws IOException {
        for (Record illuminaAdpcRecord : illuminaAdpcRecords) {
            illuminaAdpcRecord.write(outputStream);
        }
    }

    @Override
    public void close() throws Exception {
        outputStream.close();
    }

    public static class Record {
        final short aIntensity;
        final short bIntensity;
        final float aNormalizedIntensity;
        final float bNormalizedIntensity;
        final float gcScore;
        final IlluminaGenotype genotype;

        public Record(short aIntensity, short bIntensity, float aNormalizedIntensity, float bNormalizedIntensity, float gcScore, IlluminaGenotype genotype) {
            this.aIntensity = aIntensity;
            this.bIntensity = bIntensity;
            this.aNormalizedIntensity = aNormalizedIntensity;
            this.bNormalizedIntensity = bNormalizedIntensity;
            this.gcScore = gcScore;
            this.genotype = genotype;
        }

        public void write(final DataOutputStream outputStream) throws IOException {
            InfiniumDataFile.writeShort(outputStream, aIntensity);
            InfiniumDataFile.writeShort(outputStream, bIntensity);
            InfiniumDataFile.writeFloat(outputStream, aNormalizedIntensity);
            InfiniumDataFile.writeFloat(outputStream, bNormalizedIntensity);
            InfiniumDataFile.writeFloat(outputStream, gcScore);
            InfiniumDataFile.writeShort(outputStream, genotype.value);
        }
    }
}
