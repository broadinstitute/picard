/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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
package picard.vcf;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.time.Instant;

/**
 * Tests for {@link LiftoverVcf#formatProgressEvent(Instant, String, long, long, String, long)} and
 * {@link LiftoverVcf#escapeJsonString(String)}, the helpers backing LiftoverVcf's PROGRESS_JSON
 * NDJSON sidecar output.
 */
public class LiftoverVcfProgressJsonTest {

    private static final Instant FIXED = Instant.parse("2026-05-16T12:34:56Z");

    @Test
    public void testFormatProgressEventReadStage() {
        final String line = LiftoverVcf.formatProgressEvent(
                FIXED, "read", 1_000_000L, 12_345L, "chr3:148946422", 60L);
        Assert.assertEquals(line,
                "{\"timestamp\":\"2026-05-16T12:34:56Z\","
                        + "\"stage\":\"read\","
                        + "\"records_processed\":1000000,"
                        + "\"records_rejected\":12345,"
                        + "\"last_position\":\"chr3:148946422\","
                        + "\"elapsed_seconds\":60}");
    }

    @Test
    public void testFormatProgressEventNullLastPositionSerializesAsJsonNull() {
        final String line = LiftoverVcf.formatProgressEvent(
                FIXED, "read_complete", 0L, 0L, null, 0L);
        Assert.assertTrue(line.contains("\"last_position\":null"),
                "Expected JSON null for last_position when null was passed in, got: " + line);
        Assert.assertFalse(line.contains("\"last_position\":\"null\""),
                "Did not expect last_position to be the string \"null\": " + line);
    }

    @Test
    public void testFormatProgressEventLineHasNoNewline() {
        // NDJSON / JSON Lines: each event must be exactly one line. Newlines come from the
        // calling println(), not the event content itself.
        final String line = LiftoverVcf.formatProgressEvent(
                FIXED, "read", 100L, 0L, "chr1:1", 1L);
        Assert.assertFalse(line.contains("\n"), "Event line must not contain newline: " + line);
        Assert.assertFalse(line.contains("\r"), "Event line must not contain CR: " + line);
    }

    @Test
    public void testFormatProgressEventTimestampTruncatedToSeconds() {
        // Sub-second precision is dropped so output is deterministic across platforms.
        final Instant subSecond = Instant.parse("2026-05-16T12:34:56.789Z");
        final String line = LiftoverVcf.formatProgressEvent(
                subSecond, "read", 1L, 0L, "chr1:1", 0L);
        Assert.assertTrue(line.contains("\"timestamp\":\"2026-05-16T12:34:56Z\""),
                "Expected timestamp truncated to second precision, got: " + line);
    }

    @DataProvider(name = "jsonEscapeCases")
    public Object[][] jsonEscapeCases() {
        return new Object[][]{
                {"plain", "plain"},
                {"with \"quotes\"", "with \\\"quotes\\\""},
                {"back\\slash", "back\\\\slash"},
                {"new\nline", "new\\nline"},
                {"tab\there", "tab\\there"},
                {"ctrlchar", "ctrl\\u0001char"},
                {"", ""},
        };
    }

    @Test(dataProvider = "jsonEscapeCases")
    public void testEscapeJsonString(final String input, final String expected) {
        Assert.assertEquals(LiftoverVcf.escapeJsonString(input), expected);
    }

    @Test
    public void testFormatProgressEventEscapesStageAndPosition() {
        final String line = LiftoverVcf.formatProgressEvent(
                FIXED, "stage\"with\"quotes", 1L, 0L, "contig\\with\\back:1", 0L);
        Assert.assertTrue(line.contains("\"stage\":\"stage\\\"with\\\"quotes\""),
                "Stage should be escaped, got: " + line);
        Assert.assertTrue(line.contains("\"last_position\":\"contig\\\\with\\\\back:1\""),
                "Position should be escaped, got: " + line);
    }
}
