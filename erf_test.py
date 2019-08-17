#!/usr/bin/python

from scipy.special import erf

#-------------------------------------------------------------------------
# open outputfile
outfile = open("test.out", mode='w')

a = erf(1)
outfile.write(str(a))

outfile.close
