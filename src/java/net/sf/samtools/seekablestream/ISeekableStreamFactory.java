package net.sf.samtools.seekablestream;

import java.io.IOException;
import java.net.URL;

/**
 * Factory for creating {@link SeekableStream}s based on URLs/paths.
 * Implementations can be set as the default with {@link SeekableStreamFactory#setInstance(ISeekableStreamFactory)}
 * @author jacob
 * @date 2013-Oct-24
 */
public interface ISeekableStreamFactory {

    public SeekableStream getStreamFor(URL url) throws IOException;

    public SeekableStream getStreamFor(String path) throws IOException;
}
