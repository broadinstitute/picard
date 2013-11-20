#! /bin/bash
#
# The MIT License
#
# Copyright (c) 2013 The Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

# Build libIntelDeflater.so, the JNI library that wraps Intel IPP compression library.
# Note that this is not built as part of standard release process.  Rather, it is built manually and then
# copied to Picard-public/lib/jni.

# Assumes OpenJDK exists at $OPENJDK.  I used openjdk-7-fcs-src-b147-27_jun_2011.zip
# Assumes that Picard-public java sources have been compiled
# Assumes IPP8_CODE_SAMPLES_DIR points to Intel IPP sample code built with -fPIC
set -e

if [ "$OPENJDK" = "" ]
then echo "ERROR: OPENJDK environment variable not defined." >&2
     exit 1
fi

if [ "$IPP8_CODE_SAMPLES_DIR" = "" ]
then echo "ERROR: IPP8_CODE_SAMPLES_DIR environment variable not defined." >&2
     exit 1
fi

rootdir=$(dirname $(dirname $(dirname $0)))


builddir=$rootdir/lib_build
rm -rf $builddir
mkdir -p $builddir

# Create JNI C header file
javah -jni -classpath $rootdir/classes -d $builddir net.sf.samtools.util.zip.IntelDeflater

# Compile source and create library.
gcc -I$builddir -I$JAVA_HOME/include/ -I$JAVA_HOME/include/linux/ -I$OPENJDK/jdk/src/share/native/common/ \
-I$OPENJDK/jdk/src/solaris/native/common/ -c -O3 -fPIC IntelDeflater.c
gcc  -shared -o $builddir/libIntelDeflater.so IntelDeflater.o  -L${IPP8_CODE_SAMPLES_DIR}/__cmake/data-compression.intel64.make.static.release/__lib/release \
-lzlib  -lstdc++ -Wl,-Bstatic  -lbfp754  -ldecimal  -liomp5  -liompstubs5  -lipgo  -lippac  -lippcc  -lippch  -lippcv  \
-lippdc  -lippdi  -lippgen  -lippi  -lippj  -lippm  -lippr  -lippsc  -lippvc  -lippvm  -lirng  -lmatmul  -lpdbx  \
-lpdbxinst  -lsvml  -lipps  -limf  -lirc  -lirc_s  -lippcore -Wl,-Bdynamic



