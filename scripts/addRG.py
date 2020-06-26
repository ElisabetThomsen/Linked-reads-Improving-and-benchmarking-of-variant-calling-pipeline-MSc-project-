#!/usr/bin/env python3

"""
Author: Elisabet Thomsen
Date: 19-03-2020
Description: Program that adds lane information to SAM file.
Handy for e.g. EMA aligner.
In the header it adds what you input on the commandline.
In the reads, it finds the lanenumber as the number that comes after the third ":" in the readname.
Input: SAM file. On commandline: lanes to add and were to add them.
Separate lanes with "," (no whitespace!) and the same for where to add them.
Output: SAM file with added lane information.
Usage: ./<scriptname> <infile> <outfile> <lanes to add to RG> <which RGs to add lanes to>
e.g.
./addRG.py rgtest.sam rgtest_out.sam :4,:5 1,3
where the original RG is:
ID:K00126 PL:illumina PU:@K00126_334_HGKNJBBXX SM:NA12878
and we want it to add ":4" and ":5" to ID and PU. Result:
ID:K00126:4 PL:illumina PU:@K00126_334_HGKNJBBXX:4 SM:NA17878
ID:K00126:5 PL:illumina PU:@K00126_334_HGKNJBBXX:5 SM:NA12878
"""

import sys, copy

# Check command line
if len(sys.argv) != 5:
    print('Usage: ./<scriptname> <infile> <outfile> <lanes to add to RG> <which RGs to add lanes to>')
    sys.exit(1)

infilename = sys.argv[1]
outfilename = sys.argv[2]

# Open files
try:
    infile = open(infilename, 'r')
    outfile = open(outfilename, 'w')
except IOError as err:
    print('Cant open files:', str(err))
    sys.exit(1)

# Initialize
lanes = sys.argv[3].split(',')
addto = sys.argv[4].split(',')
printchunk = ''
chunksize = 100000

# Convert into integer
for i in range(len(addto)):
    addto[i] = int(addto[i])

# Iterate through file
for line in infile:
    # Find RG line and replace it with a new RG line per added lane
    if line.startswith('@RG'):
        RGline = line.split()
        for lane in lanes:
            RGline2 = copy.deepcopy(RGline)
            for RG in addto:
                RGline2[RG] += lane
            RGline2 = '\t'.join(RGline2)
            printchunk += RGline2 + '\n'
    # Find lines that not are header
    elif not line.startswith('@'):
        linelist = line.split()
        # Find lane number for this read
        lane = linelist[0].split(':')[3]
        count = 0
        # Find RG tag in this line
        for tag in reversed(linelist):
            count += 1
            if tag.startswith('RG'):
                # Add lane number to RG tag
                tag += ":" + lane
                linelist[-count] = tag
                break
        line = '\t'.join(linelist)
        printchunk += line + '\n'

    # Else just add line to printchunk
    else:
        printchunk += line

    # When chunk is chunksize, print it
    if len(printchunk) > chunksize:
        outfile.write(printchunk)
        printchunk = ''

# Print last chunk
outfile.write(printchunk)
printchunk = ''

# Close files
infile.close()
outfile.close()
