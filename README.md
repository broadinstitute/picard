[![Coverage Status](https://coveralls.io/repos/github/broadinstitute/picard/badge.svg?branch=master)](https://coveralls.io/github/broadinstitute/picard?branch=master)
[![Build Status](https://travis-ci.org/broadinstitute/picard.svg?branch=master)](https://travis-ci.org/broadinstitute/picard)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/broadinstitute/picard/blob/master/LICENSE.txt)

A set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats.  

Picard is implemented using the HTSJDK Java library [HTSJDK][1] to support
accessing file formats that are commonly used for high-throughput
sequencing data such as [SAM][2] and [VCF][3].  

As of version 2.0.1 (Nov. 2015) Picard requires Java 1.8 (jdk8u66). The last version to support Java 1.7 was release 1.141.

#### Building Picard

* First, clone the repo:
```
    git clone https://github.com/broadinstitute/picard.git
    cd picard/
```

* Picard is now built using [gradle](http://gradle.org/). A wrapper script (`gradlew`) is included which will download the appropriate version of gradle on the first invocation.
    
* To build a fully-packaged, runnable Picard jar with all dependencies included, run:
```
    ./gradlew shadowJar
```

* The resulting jar will be in `build/libs`. To run it, the command is:
```
    java -jar build/libs/picard.jar
    
    or
    
    java -jar build/libs/picard-<VERSION>-all.jar 
```    

    
* To build a jar containing only Picard classes (without its dependencies), run:
```
    ./gradlew jar
```    
    
* To clean the build directory, run:
```
    ./gradlew clean
```

#### Running Tests

* To run all tests, the command is:
```
    ./gradlew test
```

* To run a specific test, the command is:
```
    ./gradlew test -Dtest.single=TestClassName 
```

#### Changing the released version of HTSJDK that Picard depends on

To switch Picard's HTSJDK dependency to a different released version:

* Open `build.gradle`
* Edit VERSION in the following line to be a different released version of HTSJDK. HTSJDK releases are listed [here](https://github.com/samtools/htsjdk/releases)
```
    final htsjdkVersion = System.getProperty('htsjdk.version', 'VERSION')`
```
* Open a pull request with this change

#### Building Picard with a Custom Version of HTSJDK

During development in Picard, it is sometimes necessary to build locally against an unreleased version or branch of HTSJDK. 

* To build against an unreleased version of HTSJDK's master branch:
    * Go to the [Broad artifactory](https://artifactory.broadinstitute.org/artifactory/simple/libs-snapshot-local/com/github/samtools/htsjdk/), where continuous snapshots of HTSJDK's master branch are published, and select the version you want to use. For example, `2.5.1-9-g5740ca1-SNAPSHOT`. You can search by tag or short git commit hash.
    * In your Picard clone, run `./gradlew shadowJar -Dhtsjdk.version=VERSION`, where VERSION is the version of the HTSJDK master branch snapshot you want to use.
    
* To build against a version of HTSJDK that has *not* yet been merged into HTSJDK's master branch:
    * Clone [HTSJDK](https://github.com/samtools/htsjdk/), and in your clone check out the tag or branch you want to build Picard with.
    * Run `./gradlew install printVersion` in your htsjdk clone to install that version to your local maven repository. Take note of the version number that gets printed at the end.
    * Switch back to your Picard clone, and run `./gradlew shadowJar -Dhtsjdk.version=VERSION`, where VERSION is the version of HTSJDK you installed to your local maven repository.

#### Releasing Picard

Full instructions on how to create a new release of 
Picard are [here](https://github.com/broadinstitute/picard/wiki/How-to-release-Picard)

#### Path providers

Picard has limited support for reading from Path providers. 
Currently only google's api is supported, and only a few tools support this.
To run with this support you need to compile the cloudJar target with gradle:
```bash
./gradlew cloudJar

```
then run picard as follows:

```bash
java -jar build/lib/picardcloud.jar <Picard arguments starting from program>
```
For example:

```bash 
java -jar build/lib/picardcloud.jar CrosscheckFingerprints \
   I=gs://sample1.vcf \
   I=gs://sample2.vcf \
   CROSSCHECK_BY=FILE \
   H=Haplotype_db.txt \
   O=crosscheck.out
```

Alternatively, you can run the tool via the [GATK](https://software.broadinstitute.org/gatk/download/) which bundles the Google-Cloud
jar, and should thus "Just Work".

#### GA4GH API

It's also possible to build a version of Picard that supports reading from
GA4GH API, e.g. Google Genomics:

* Fetch [gatk-tools-java](https://github.com/googlegenomics/gatk-tools-java) 
 
```git clone https://github.com/googlegenomics/gatk-tools-java```

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

----

Picard is migrating to semantic versioning (http://semver.org/). We will eventually adhere to it strictly and bump our major version whenever there are breaking changes to our API, but until we more clearly define what constitutes our official API, clients should assume that every release potentially contains at least minor changes to public methods.

Please see the [Picard Documentation](http://broadinstitute.github.io/picard) for more information.

[1]: http://github.com/samtools/htsjdk
[2]: http://samtools.sourceforge.net
[3]: http://vcftools.sourceforge.net/specs.html
