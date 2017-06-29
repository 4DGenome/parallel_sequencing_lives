#!/usr/bin/env python

# Import python modules and functions
import os, sys

# read resolution in bp passed as a script parameter
resolution_bp = int(sys.argv[1])

# converto to Kb or Mb
def nice(reso):

    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)

    return '%dkb' % (reso / 1000)


print nice(resolution_bp)