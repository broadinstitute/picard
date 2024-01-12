package picard.cmdline.argumentcollections;

import com.google.common.base.Strings;
import org.broadinstitute.barclay.argparser.Argument;

/**
 * Argument collection to encapsulate special behavior of the REQUESTER_PAYS_PROJECT argument
 */
public class RequesterPaysArgumentCollection {
    /** The System property which acts as a default for REQUESTER_PAYS_PROJECT */
    public static final String PROPERTY_NAME = "picard.googleProjectForRequesterPays";
    public static final String REQUESTER_PAYS_PROJECT_FULL_NAME = "REQUESTER_PAYS_PROJECT";

    @Argument(doc="Google project for access to 'requester pays' buckets and objects.  " +
            "If this is not specified then value of the system property " + PROPERTY_NAME + " acts as the default.",
            common = true, optional = true,
    fullName = REQUESTER_PAYS_PROJECT_FULL_NAME)
    public String requesterPaysProject = null;

    public String getProjectForRequesterPays() {
        final String value = ! Strings.isNullOrEmpty(requesterPaysProject)
                ? requesterPaysProject
                : getSystemProperty();
        return Strings.isNullOrEmpty(value) ? null : value; // "" -> null
    }

    private String getSystemProperty() {
        return System.getProperty(PROPERTY_NAME);
    }

    public String getDescription(){
         final String value = getProjectForRequesterPays();
         if(!Strings.isNullOrEmpty(requesterPaysProject)){
             return "Requester Pays Project set by argument: " + value;
         } else if( !Strings.isNullOrEmpty(getSystemProperty())){
             return "Requester Pays Project set by system property: " + value;
         } else {
             return "Requester Pays Project not set.";
         }
    }

}
