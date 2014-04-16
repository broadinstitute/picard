/*
* Copyright (c) 2014 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.variantcontext.writer;

import java.io.*;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import net.sf.samtools.Defaults;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.BlockCompressedOutputStream;
import org.broad.tribble.AbstractFeatureReader;
import org.broadinstitute.variant.VariantBaseTest;

public class VariantContextWriterBuilderUnitTest extends VariantBaseTest {
    private SAMSequenceDictionary dictionary;

    private final File vcf = new File("test.vcf");
    private final File vcfIdx = new File("test.vcf.idx");
    private final File vcfMD5 = new File("test.vcf.md5");
    private final File bcf = new File("test.bcf");
    private final File bcfIdx = new File("test.bcf.idx");
    private final File unknown = new File("test.unknown");

    @BeforeSuite
    public void before() throws IOException {
        dictionary = createArtificialSequenceDictionary();
        vcf.deleteOnExit();
        vcfIdx.deleteOnExit();
        vcfMD5.deleteOnExit();
        bcf.deleteOnExit();
        bcfIdx.deleteOnExit();
        unknown.deleteOnExit();
    }

    @Test
    public void testSetOutputFile() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary);

        VariantContextWriter writer = builder.setOutputFile(vcf.getAbsolutePath()).build();
        Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputFile VCF String");
        Assert.assertFalse(((VCFWriter) writer).getOutputStream() instanceof BlockCompressedOutputStream, "testSetOutputFile VCF String was compressed");

        writer = builder.setOutputFile(vcf).build();
        Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputFile VCF File");
        Assert.assertFalse(((VCFWriter)writer).getOutputStream() instanceof BlockCompressedOutputStream, "testSetOutputFile VCF File was compressed");

        for (final String extension : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS) {
            final String filename = "test.vcf" + extension;
            final File file = new File(filename);
            file.deleteOnExit();

            writer = builder.setOutputFile(filename).build();
            Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputFile " + extension + " String");
            Assert.assertTrue(((VCFWriter) writer).getOutputStream() instanceof BlockCompressedOutputStream, "testSetOutputFile " + extension + " String was not compressed");

            writer = builder.setOutputFile(file).build();
            Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputFile " + extension + " File");
            Assert.assertTrue(((VCFWriter) writer).getOutputStream() instanceof BlockCompressedOutputStream, "testSetOutputFile " + extension + " File was not compressed");
        }

        writer = builder.setOutputFile(bcf).build();
        Assert.assertTrue(writer instanceof BCF2Writer, "testSetOutputFile BCF String");

        writer = builder.setOutputFile(bcf.getAbsolutePath()).build();
        Assert.assertTrue(writer instanceof BCF2Writer, "testSetOutputFile BCF File");
    }

    @Test
    public void testSetOutputFileType() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputFile(unknown);

        VariantContextWriter writer = builder.setOutputFileType(VariantContextWriterBuilder.OutputType.VCF).build();
        Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputFileType VCF");
        Assert.assertFalse(((VCFWriter) writer).getOutputStream() instanceof BlockCompressedOutputStream, "testSetOutputFileType VCF was compressed");

        writer = builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).build();
        Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputFile BLOCK_COMPRESSED_VCF");
        Assert.assertTrue(((VCFWriter) writer).getOutputStream() instanceof BlockCompressedOutputStream, "testSetOutputFileType BLOCK_COMPRESSED_VCF was not compressed");

        writer = builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BCF).build();
        Assert.assertTrue(writer instanceof BCF2Writer, "testSetOutputFileType BCF");
    }

    @Test
    public void testSetOutputStream() {
        final OutputStream stream = new ByteArrayOutputStream();

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .setOutputStream(stream);

        VariantContextWriter writer = builder.build();
        Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputStream default");

        writer = builder.setOption(Options.FORCE_BCF).build();
        Assert.assertTrue(writer instanceof BCF2Writer, "testSetOutputStream FORCE_BCF");

        writer = builder.setOutputVCFStream(stream).build();
        Assert.assertTrue(writer instanceof BCF2Writer, "testSetOutputStream FORCE_BCF 2");

        writer = builder.unsetOption(Options.FORCE_BCF).build();
        Assert.assertTrue(writer instanceof VCFWriter, "testSetOutputStream VCF");

        writer = builder.setOutputBCFStream(stream).build();
        Assert.assertTrue(writer instanceof BCF2Writer, "testSetOutputStream BCF");
    }

    @Test
    public void testAsync() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputFile(vcf);

        VariantContextWriter writer = builder.build();
        Assert.assertEquals(writer instanceof AsyncVariantContextWriter, Defaults.USE_ASYNC_IO, "testAsync default");

        writer = builder.setOption(Options.USE_ASYNC_IO).build();
        Assert.assertTrue(writer instanceof AsyncVariantContextWriter, "testAsync option=set");

        writer = builder.unsetOption(Options.USE_ASYNC_IO).build();
        Assert.assertFalse(writer instanceof AsyncVariantContextWriter, "testAsync option=unset");
    }

    @Test
    public void testBuffering() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputFile(vcf)
                .unsetOption(Options.INDEX_ON_THE_FLY);     // so the potential BufferedOutputStream is not wrapped in a PositionalOutputStream

        VariantContextWriter writer = builder.build();
        Assert.assertTrue(((VCFWriter)writer).getOutputStream() instanceof BufferedOutputStream, "testBuffering was not buffered by default");

        writer = builder.unsetBuffering().build();
        Assert.assertFalse(((VCFWriter) writer).getOutputStream() instanceof BufferedOutputStream, "testBuffering was buffered when unset");

        writer = builder.setBuffer(8192).build();
        Assert.assertTrue(((VCFWriter) writer).getOutputStream() instanceof BufferedOutputStream, "testBuffering was not buffered when set");
    }

    @Test
    public void testMD5() {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputFile(vcf);

        VariantContextWriter writer = builder.build();
        writer.close();
        Assert.assertEquals(vcfMD5.exists(), Defaults.CREATE_MD5, "MD5 default setting not respected");

        if (vcfMD5.exists())
            vcfMD5.delete();

        writer = builder.setCreateMD5().build();
        writer.close();
        Assert.assertTrue(vcfMD5.exists(), "MD5 not created when requested");
        vcfMD5.delete();

        writer = builder.unsetCreateMD5().build();
        writer.close();
        Assert.assertFalse(vcfMD5.exists(), "MD5 created when not requested");

        writer = builder.setCreateMD5(false).build();
        writer.close();
        Assert.assertFalse(vcfMD5.exists(), "MD5 created when not requested via boolean parameter");

        writer = builder.setCreateMD5(true).build();
        writer.close();
        Assert.assertTrue(vcfMD5.exists(), "MD5 not created when requested via boolean parameter");
        vcfMD5.delete();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidImplicitFileType() {
        new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputFile("test.bam")
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSetInvalidFileType() {
        new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputFile("test.bam")
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF_STREAM)
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidSetFileTypeForStream() {
        new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputStream(new ByteArrayOutputStream())
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnsupportedIndexOnTheFlyForStreaming() {
        new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputStream(new ByteArrayOutputStream())
                .setOption(Options.INDEX_ON_THE_FLY)
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUnsupportedDefaultIndexOnTheFlyForStreaming() {
        new VariantContextWriterBuilder()
                .setReferenceDictionary(dictionary)
                .setOutputStream(new ByteArrayOutputStream())
                .build();
    }

}