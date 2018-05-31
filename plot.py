#!/usr/bin/python
##########################################################################
#                       Plot ELDEST results and make video               #
##########################################################################
# Purpose:                                                               #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import glob
import os

filenames = glob.glob('*.dat')
n_files = len(filenames)
times = []
for i in range(0,n_files):
    temp = filenames[i]
    nodat = temp.replace('.dat','')
    if (nodat[0] == 'm'):
        nodat = nodat.replace('m','-')
    fl = float(nodat)
    times.append(fl)

times.sort()

fn_red = []
filenames.sort()
for i in range(0,n_files):
    temp = filenames[i]
    nodat = temp.replace('.dat','')
    fn_red.append(nodat)

#print filenames
    infile = open('mk_images.gnu', 'r')
    outfile = open('tmpfile', 'w')

    replacements = {'number':str(i), 'filename':filenames[i] }

    for line in infile:
       if 'FloatBarrier' not in line:
           for src, target in replacements.iteritems():
               line = line.replace(src, target)
           outfile.write(line)

    infile.close()
    outfile.close()

    os.system('gnuplot tmpfile')
    os.remove('tmpfile')

os.system("ffmpeg -f image2 -pattern_type glob -r 10.0 -i 'image.*.png' -qscale 1 out.mp4")
