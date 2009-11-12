#!/usr/bin/env python
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2008 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
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
    ("read is PCR or optical duplicate", 0x400)
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
    
