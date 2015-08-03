[![Build Status](https://travis-ci.org/broadinstitute/picard.svg?branch=master)](https://travis-ci.org/broadinstitute/picard)

A set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats.  

Picard is implemented using the HTSJDK Java library[HTSJDK][1], supporting
accessing of common file formats, such as [SAM][2] and [VCF][3], used for high-throughput
sequencing data.  

It's also possible to build a version of Picard that supports reading from
GA4GH API, e.g. Google Genomics:

* Fetch [gatk-tools-java](https://github.com/gatk-tools-java) 
 
```git clone https://github.com/gatk-tools-java```

* Build gatk-tools-java: 

```gatk-tools-java$ mvn compile package```

* Copy the resulting jar into Picard ```lib/``` folder:
```
gatk-tools-java$ mkdir ../picard/lib/gatk-tools-java
gatk-tools-java$ cp target/gatk-tools-java*minimized.jar ../picard/lib/gatk-tools-java/
```

* Build Picard version with GA4GH support: 

```picard$ ant -lib lib/ant package-commands-ga4gh```

* If you have not yet worked with Google Genomics API and need to set up authentication, please follow the instructions [here](https://cloud.google.com/genomics/install-genomics-tools#authenticate) to set up credentials and obtain ```client_secrets.json``` file.


* You can now run 
```
java -jar dist/picard.jar ViewSam \
INPUT=https://www.googleapis.com/genomics/v1beta2/readgroupsets/CK256frpGBD44IWHwLP22R4/ \
GA4GH_CLIENT_SECRETS=../client_secrets.json
```

* To run using GRPC as opposed to REST Genomics API implementation (which is much faster) use the following command that utilizes ALPN jars that come with gatk-tools-java and enables GRPC support:
```
java \
-Xbootclasspath/p:../gatk-tools-java/lib/alpn-boot-8.1.3.v20150130.jar \
-Dga4gh.using_grpc=true \
-jar dist/picard.jar ViewSam \
INPUT=https://www.googleapis.com/genomics/v1beta2/readgroupsets/CK256frpGBD44IWHwLP22R4/ \
GA4GH_CLIENT_SECRETS=../client_secrets.json
```
For Java 7 (as opposed to 8) use ```alpn-boot-7.1.3.v20150130.jar```.


Please see the [Picard Documentation](http://broadinstitute.github.io/picard) for more information.

[1]: http://github.com/samtools/htsjdk
[2]: http://samtools.sourceforge.net
[3]: http://vcftools.sourceforge.net/specs.html