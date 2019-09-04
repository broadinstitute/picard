package picard.annotation;

import java.util.*;

public class Gff3Codec extends GffCodec {
    public final static Set<String> GFF3_FILE_EXTENSIONS = new HashSet<>(Arrays.asList("gff", "gff3"));

    public Gff3Codec() {
        super(GFF3_FILE_EXTENSIONS);
    }

    @Override
    public Map<String,String> parseAttributes(final String attributesString) {
        final Map<String, String> attributes = new LinkedHashMap<>();
        final String[] splitLine = attributesString.split(";");
        for(String attribute : splitLine) {
            final String[] key_value = attribute.split("=");
            if (key_value.length<2) {
                continue;
            }
            attributes.put(key_value[0], key_value[1]);
        }
        return attributes;
    }

    @Override
    public String getFirstLineStart() {
        return "##gff-version 3";
    }
}