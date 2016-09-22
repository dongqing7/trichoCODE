#!/usr/bin/python
import sys

if len(sys.argv) != 4:
    print "usage: "+sys.argv[0] + " start_mark file outfile"
    sys.exit(255)

mark = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

wh = open(outfile, 'w')

for line in open(infile):
    if line.startswith(mark):
        continue

    wh.write(line)

wh.close()

