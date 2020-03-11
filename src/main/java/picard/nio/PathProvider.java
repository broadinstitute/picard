package picard.nio;

/**
 * A class whose purpose is to initialize the various plugins that provide Path support.
 * <p>
 * At the moment only google-GCS support is available, but it should be easy to follow the pattern and add support for
 * other Path providers.
 */
public enum PathProvider {
    GCS {
        @Override
        public boolean initialize() {
            try {
                GoogleStorageUtils.initialize();
                return true;
            } catch (NoClassDefFoundError e) {
                return false;
            }
        }
    };

    public final boolean isAvailable;

    protected abstract boolean initialize();

    PathProvider() {
        //noinspection AbstractMethodCallInConstructor
        isAvailable = initialize();
    }
}
