#!/usr/bin/env python3

"""
Description: Program that split TP and FP variants into seperate files.
It expect input from gatk VariantsToTable, and the Truth info to be in the 11th colum.
Usage: ./scripname.py <infilename>
"""

import sys

infilename = sys.argv[1]
TPfilename = infilename + '.TP'
FPfilename = infilename + '.FP'

infile = open(infilename, 'rb')
TPfile = open(TPfilename, 'wb')
FPfile = open(FPfilename, 'wb')

chunksize = 1024*1024
FP_chunk = b''
TP_chunk = b''

# Get the header
for line in infile:
    FP_chunk += line
    TP_chunk += line
    break

# Iterate through lines, and print TP and FP into seperate files
for line in infile:
    if line.split()[11] == b'NA':
        FP_chunk += line
        if len(FP_chunk) > chunksize:
            FPfile.write(FP_chunk)
            FP_chunk = b''
    else:
        TP_chunk += line
        if len(TP_chunk) > chunksize:
            TPfile.write(TP_chunk)
            TP_chunk = b''

# Print the last chunks
FPfile.write(FP_chunk)
TPfile.write(TP_chunk)

# Close files
infile.close()
TPfile.close()
FPfile.close()
