package picard.util.help;

import org.broadinstitute.barclay.help.DefaultDocWorkUnitHandler;
import org.broadinstitute.barclay.help.DocWorkUnit;

import org.broadinstitute.barclay.help.HelpDoclet;

/**
 * The Picard Documentation work unit handler class that is the companion to PicardHelpDoclet.
 *
 * NOTE: Methods in this class are intended to be called by Gradle/Javadoc only, and should not be called
 * by methods that are used by the Picard runtime, as this class assumes a dependency on com.sun.javadoc classes
 * which may not be present.
 */
public class PicardHelpDocWorkUnitHandler extends DefaultDocWorkUnitHandler {

    private final static String PICARD_JAVADOC_TAG_PREFIX = "picard"; // prefix for custom javadoc tags used by Picard

    private final static String PICARD_FREEMARKER_TEMPLATE_NAME = "generic.template.html";

    public PicardHelpDocWorkUnitHandler(final HelpDoclet doclet) {
        super(doclet);
    }
    /**
     * @return Prefix for custom picard tags that should be lifted from the javadoc and stored in the
     * FreeMarker map. These will be available in the template returned by {@link #getTemplateName}.
     */
    @Override
    protected String getTagFilterPrefix() { return PICARD_JAVADOC_TAG_PREFIX; }

    /**
     * @param workUnit the classdoc object being processed
     * @return the name of a the freemarker template to be used for the class being documented.
     * Must reside in the folder passed to the Barclay Doclet via the "-settings-dir" parameter to
     * Javadoc.
     */
    @Override
    public String getTemplateName(final DocWorkUnit workUnit) { return PICARD_FREEMARKER_TEMPLATE_NAME; }
}
