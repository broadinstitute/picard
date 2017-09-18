package picard.util;

import org.broadinstitute.barclay.utils.Utils;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * Utility for loading properties files from resources.
 */
public class PropertyUtils {

    /**
     * Attempt to load a Properties object from an on-disk properties file.
     * @param propertyFilePath name of the properties file to load. Must have a .properties extension
     * @param clazz class used to obtain a class loader to use to locate the properties file
     * @return null if the files doesn't exist or isn't readable, otherwise a Properties object
     */
    public static Properties loadPropertiesFile(final String propertyFilePath, final Class<?> clazz) {
        Utils.nonNull(propertyFilePath);

        try (final InputStream inputStream = clazz.getClassLoader().getResourceAsStream(propertyFilePath)) {
            if (inputStream != null) {
                final Properties properties = new Properties();
                properties.load(inputStream);
                return properties;
            } else {
                return null;
            }
        } catch (IOException ex) {
            throw new RuntimeException(String.format("IOException loading properties file %s", propertyFilePath), ex);
        }
    }

}
