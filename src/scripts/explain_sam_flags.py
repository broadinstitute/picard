#!/usr/bin/env python

# The MIT License
#
# Copyright (c) $today.year The Broad Institute
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
#
# $Header$

"""usage %prog decimal-flag [decimal-flag...]

Explain each flag on the command line in plain English
"""

from __future__ import division
import sys

lstFlags = [
    ("read paired", 0x1),
    ("read mapped in proper pair", 0x2),
    ("read unmapped", 0x4),
    ("mate unmapped", 0x8),
    ("read reverse strand", 0x10),
    ("mate reverse strand", 0x20),
    ("first in pair", 0x40),
    ("second in pair", 0x80),
    ("not primary alignment", 0x100),
    ("read fails platform/vendor quality checks", 0x200),
    ("read is PCR or optical duplicate", 0x400),
    ("supplementary alignment", 0x800)
    ]
    

def explain_sam_flags(iFlags):
    print iFlags, ":"
    for strFlagName, iMask in lstFlags:
        if iFlags & iMask:
            print "\t" + strFlagName

def main(argv=None):
    if argv is None:
        argv = sys.argv

    for strArg in argv[1:]:
        explain_sam_flags(int(strArg))

if __name__ == "__main__":
    sys.exit(main())
    
