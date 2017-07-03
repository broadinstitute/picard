package picard.annotation;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * Interface for receiving records from file as stream
 */
public interface Reader<T> {

    Path path();

    Function<String, Optional<T>> lineMapper();

    static <T> Reader<T> of(final Path path, final Function<String, Optional<T>> mapper) {
        return new Reader<T>() {
            @Override
            public Path path() {
                return path;
            }

            @Override
            public Function<String, Optional<T>> lineMapper() {
                return mapper;
            }
        };
    }

     default Stream<T> records() throws IOException {
        return Files.lines(path())
                .map(lineMapper())
                .filter(record -> {
                    final boolean present = record.isPresent();
                    if (!present)
                        throw new AnnotationException(String.format("Unable to produce record from file %s:", path()));
                    return true;
                })
                .map(Optional::get);
    }
}
