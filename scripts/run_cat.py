#!/usr/bin/python

import commands
import os, sys

if len(sys.argv) != 2:
    print "usage: " + sys.argv[0] + " list"

file_list = sys.argv[1]

for line in open(file_list):
    line = line.strip()
    if line == "":
        continue
    splits = line.split(' ')
    for contig in splits:
        cmd = "cat " + contig + "/evm.out.gff3 >> evm.gff3"
        result = os.system(cmd)
        if result != 0:
            print commands.getstatusoutput(cmd)
            exit(255)
        else:
            print cmd

