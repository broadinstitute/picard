package org.broad.tribble;

import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.readers.LocationAware;
import org.broad.tribble.readers.PositionalBufferedStream;

import java.io.IOException;
import java.io.InputStream;

/**
 * Implements common methods of {@link FeatureCodec}s that read from {@link PositionalBufferedStream}s.
 * @author mccowan
 */
abstract public class BinaryFeatureCodec<T extends Feature> implements FeatureCodec<T, PositionalBufferedStream> {
    @Override
    public PositionalBufferedStream makeSourceFromStream(final InputStream bufferedInputStream) {
        if (bufferedInputStream instanceof PositionalBufferedStream)
            return (PositionalBufferedStream) bufferedInputStream;
        else
            return new PositionalBufferedStream(bufferedInputStream);
    }

    /** {@link PositionalBufferedStream} is already {@link LocationAware}. */
    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        return makeSourceFromStream(bufferedInputStream);
    }

    @Override
    public void close(final PositionalBufferedStream source) {
        CloserUtil.close(source);
    }

    @Override
    public boolean isDone(final PositionalBufferedStream source) {
        try {
            return source.isDone();
        } catch (IOException e) {
            throw new RuntimeException("Failure reading from stream.", e);
        }
    }
}
