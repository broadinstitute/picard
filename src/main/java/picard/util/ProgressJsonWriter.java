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

import picard.PicardException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.time.Instant;
import java.time.temporal.ChronoUnit;
import java.util.Locale;

/**
 * Writes progress events as JSON Lines (NDJSON): one UTF-8, LF-terminated JSON object per event.
 * Gives downstream tooling a stable parsing target instead of regex-scanning human-readable logs.
 * Each event is flushed immediately so the file can be tailed while the producing tool runs.
 */
public final class ProgressJsonWriter implements AutoCloseable {

    private final File file;
    private final Writer writer;

    /** Opens {@code file} for writing, truncating any existing content. */
    public ProgressJsonWriter(final File file) {
        this.file = file;
        try {
            // Plain UTF-8 stream (not IOUtil's gzip-by-extension wrapping) so the sidecar can be
            // tailed line-by-line while the tool runs.
            this.writer = new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(file), StandardCharsets.UTF_8));
        } catch (final IOException e) {
            throw new PicardException("Could not open progress JSON file: " + file, e);
        }
    }

    /** Appends one event as a single NDJSON line, stamped with the current time, and flushes. */
    public void emit(final String stage, final long recordsProcessed, final long recordsRejected,
                     final String lastPosition, final long elapsedSeconds) {
        try {
            writer.write(formatEvent(Instant.now(), stage, recordsProcessed, recordsRejected,
                    lastPosition, elapsedSeconds));
            writer.write('\n');
            writer.flush();
        } catch (final IOException e) {
            throw new PicardException("Could not write to progress JSON file: " + file, e);
        }
    }

    @Override
    public void close() {
        try {
            writer.close();
        } catch (final IOException e) {
            throw new PicardException("Could not close progress JSON file: " + file, e);
        }
    }

    /**
     * Builds one NDJSON event line. Field order is fixed so downstream consumers can rely on a
     * stable schema. Package-private and static for unit testing.
     */
    static String formatEvent(final Instant timestamp, final String stage,
                              final long recordsProcessed, final long recordsRejected,
                              final String lastPosition, final long elapsedSeconds) {
        final StringBuilder sb = new StringBuilder(160);
        sb.append('{');
        sb.append("\"timestamp\":\"").append(timestamp.truncatedTo(ChronoUnit.SECONDS)).append("\",");
        sb.append("\"stage\":\"").append(escapeJsonString(stage)).append("\",");
        sb.append("\"records_processed\":").append(recordsProcessed).append(',');
        sb.append("\"records_rejected\":").append(recordsRejected).append(',');
        sb.append("\"last_position\":");
        if (lastPosition == null) {
            sb.append("null");
        } else {
            sb.append('"').append(escapeJsonString(lastPosition)).append('"');
        }
        sb.append(',');
        sb.append("\"elapsed_seconds\":").append(elapsedSeconds);
        sb.append('}');
        return sb.toString();
    }

    /** Escapes a string for embedding in a JSON double-quoted value. */
    static String escapeJsonString(final String s) {
        final StringBuilder sb = new StringBuilder(s.length() + 4);
        for (int i = 0; i < s.length(); i++) {
            final char c = s.charAt(i);
            switch (c) {
                case '"':  sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\b': sb.append("\\b"); break;
                case '\f': sb.append("\\f"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                default:
                    if (c < 0x20) {
                        sb.append(String.format(Locale.ROOT, "\\u%04x", (int) c));
                    } else {
                        sb.append(c);
                    }
            }
        }
        return sb.toString();
    }
}
