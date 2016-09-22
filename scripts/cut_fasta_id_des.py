#!/usr/bin/python
import sys

def cut(f_fasta):

    for line in open(f_fasta):
        line = line.strip()
        if line.startswith('>'):
            splits = line.split(' ')
            print splits[0]
        else:
            print line

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "usage: %s pasa_update.faa " % (sys.argv[1])
        exit(255)

    infile = sys.argv[1]

    cut(infile)

