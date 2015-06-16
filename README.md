[![Build Status](https://travis-ci.org/broadinstitute/picard.svg?branch=master)](https://travis-ci.org/broadinstitute/picard)

A set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats.  

Picard is implemented using the HTSJDK Java library[HTSJDK][1], supporting
accessing of common file formats, such as [SAM][2] and [VCF][3], used for high-throughput
sequencing data.  

It's also possible to build a version of Picard that supports reading from
GA4GH API, e.g. Google Genomics:
1.Fetch gatk-tools-java library from https://github.com/iliat/gatk-tools-java
2.Build: ant gatk-tools-java-picard-jar
3.Copy the resulting jar to lib/gatk-tools-java/ under Picard folder, e.g.
cp dist/gatk-tools-java-picard-1.0.jar ../picard/lib/gatk-tools-java/
4.Build Picard: ant -lib lib/ant package-commands-ga4gh
5.You can now run 
java -jar dist/picard.jar ViewSam \
INPUT=https://www.googleapis.com/genomics/v1beta2/readgroupsets/CK256frpGBD44IWHwLP22R4/ \
GA4GH_CLIENT_SECRETS=../client_secrets.json

Please see the [Picard Documentation](http://broadinstitute.github.io/picard) for more information.

[1]: http://github.com/samtools/htsjdk
[2]: http://samtools.sourceforge.net
[3]: http://vcftools.sourceforge.net/specs.html
