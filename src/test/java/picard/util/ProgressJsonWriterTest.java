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
package picard.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.time.Instant;
import java.util.List;

/**
 * Tests for {@link ProgressJsonWriter}: the NDJSON event formatter, the JSON string escaper, and
 * the on-disk encoding/termination of the sidecar.
 */
public class ProgressJsonWriterTest {

    private static final Instant FIXED = Instant.parse("2026-05-16T12:34:56Z");

    @Test
    public void testFormatEventReadStage() {
        final String line = ProgressJsonWriter.formatEvent(
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
    public void testFormatEventNullLastPositionSerializesAsJsonNull() {
        final String line = ProgressJsonWriter.formatEvent(
                FIXED, "read_complete", 0L, 0L, null, 0L);
        Assert.assertTrue(line.contains("\"last_position\":null"),
                "Expected JSON null for last_position when null was passed in, got: " + line);
        Assert.assertFalse(line.contains("\"last_position\":\"null\""),
                "Did not expect last_position to be the string \"null\": " + line);
    }

    @Test
    public void testFormatEventLineHasNoNewline() {
        // NDJSON: each event must be exactly one line. The line terminator is added by the writer,
        // not by the event content itself.
        final String line = ProgressJsonWriter.formatEvent(
                FIXED, "read", 100L, 0L, "chr1:1", 1L);
        Assert.assertFalse(line.contains("\n"), "Event line must not contain newline: " + line);
        Assert.assertFalse(line.contains("\r"), "Event line must not contain CR: " + line);
    }

    @Test
    public void testFormatEventTimestampTruncatedToSeconds() {
        // Sub-second precision is dropped so output is deterministic across platforms.
        final Instant subSecond = Instant.parse("2026-05-16T12:34:56.789Z");
        final String line = ProgressJsonWriter.formatEvent(
                subSecond, "read", 1L, 0L, "chr1:1", 0L);
        Assert.assertTrue(line.contains("\"timestamp\":\"2026-05-16T12:34:56Z\""),
                "Expected timestamp truncated to second precision, got: " + line);
    }

    @DataProvider(name = "jsonEscapeCases")
    public Object[][] jsonEscapeCases() {
        // (char) 1 is the U+0001 control character, built at runtime so the source stays ASCII.
        final String withControlChar = "ctrl" + ((char) 1) + "char";
        return new Object[][]{
                {"plain", "plain"},
                {"with \"quotes\"", "with \\\"quotes\\\""},
                {"back\\slash", "back\\\\slash"},
                {"new\nline", "new\\nline"},
                {"tab\there", "tab\\there"},
                {withControlChar, "ctrl\\u0001char"},
                {"", ""},
        };
    }

    @Test(dataProvider = "jsonEscapeCases")
    public void testEscapeJsonString(final String input, final String expected) {
        Assert.assertEquals(ProgressJsonWriter.escapeJsonString(input), expected);
    }

    @Test
    public void testFormatEventEscapesStageAndPosition() {
        final String line = ProgressJsonWriter.formatEvent(
                FIXED, "stage\"with\"quotes", 1L, 0L, "contig\\with\\back:1", 0L);
        Assert.assertTrue(line.contains("\"stage\":\"stage\\\"with\\\"quotes\""),
                "Stage should be escaped, got: " + line);
        Assert.assertTrue(line.contains("\"last_position\":\"contig\\\\with\\\\back:1\""),
                "Position should be escaped, got: " + line);
    }

    @Test
    public void testEmitWritesLfTerminatedLinesUtf8() throws IOException {
        final File f = File.createTempFile("progress", ".ndjson");
        f.deleteOnExit();
        // A non-ASCII code point (U+00E9, e-acute) must round-trip as UTF-8, not the platform
        // charset. Built from the code point so this source file stays pure ASCII.
        final String accentedPosition = "chr" + new String(Character.toChars(0x00E9)) + ":2";
        try (final ProgressJsonWriter w = new ProgressJsonWriter(f)) {
            w.emit("read", 1L, 0L, "chr1:1", 0L);
            w.emit("read_complete", 2L, 0L, accentedPosition, 1L);
        }

        final byte[] bytes = Files.readAllBytes(f.toPath());
        final String content = new String(bytes, StandardCharsets.UTF_8);
        Assert.assertFalse(content.contains("\r"), "NDJSON must use LF, not CRLF: " + content);

        final List<String> lines = Files.readAllLines(f.toPath(), StandardCharsets.UTF_8);
        Assert.assertEquals(lines.size(), 2, "Each event should be its own line");
        Assert.assertTrue(lines.get(0).contains("\"stage\":\"read\""));
        Assert.assertTrue(lines.get(1).contains(accentedPosition),
                "Non-ASCII position should round-trip as UTF-8: " + lines.get(1));
        // File must end with a single trailing LF (last line fully terminated).
        Assert.assertEquals(bytes[bytes.length - 1], (byte) '\n', "File should end with LF");
    }
}
