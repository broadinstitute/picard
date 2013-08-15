package org.broad.tribble.readers;

import net.sf.samtools.Defaults;
import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.TribbleException;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * A collection of factories for generating {@link LineReader}s.
 * 
 * @author mccowan
 */
public class LineReaderUtil {
    public enum LineReaderOption {
        ASYNCHRONOUS, SYNCHRONOUS
    }

    /**
     * Like {@link #fromBufferedStream(java.io.InputStream, org.broad.tribble.readers.LineReaderUtil.LineReaderOption)}, but the synchronicity
     * option is determined by {@link net.sf.samtools.Defaults}: if asynchronous I/O is enabled, an asynchronous line reader will be
     * returned.
     */
    public static LineReader fromBufferedStream(final InputStream stream) {
        return fromBufferedStream(stream, Defaults.USE_ASYNC_IO ? LineReaderOption.ASYNCHRONOUS : LineReaderOption.SYNCHRONOUS);
    }

    /**
     * Convenience factory for composing a LineReader from an InputStream.
     */
    public static LineReader fromBufferedStream(final InputStream bufferedStream, final LineReaderOption option) {
        final InputStreamReader bufferedInputStreamReader = new InputStreamReader(bufferedStream);
        switch (option) {
            case ASYNCHRONOUS:
                return new AsynchronousLineReader(bufferedInputStreamReader);
            case SYNCHRONOUS:
                return new LineReader() {
                    final LongLineBufferedReader reader = new LongLineBufferedReader(bufferedInputStreamReader);

                    @Override
                    public String readLine() {
                        try {
                            return reader.readLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    }

                    @Override
                    public void close() {
                        CloserUtil.close(reader);
                    }
                };
            default:
                throw new TribbleException(String.format("Unrecognized LineReaderUtil option: %s.", option));
        }
    }

}
