#!/usr/bin/python

import sys
import re
"""
"""

if len(sys.argv) != 4:
    print "usage: " + sys.argv[0] + " species_diff.gff3 species_diff.csv species_diff_fit.gff3"
    sys.exit(255)

file_gff3 = sys.argv[1]
file_csv = sys.argv[2]
file_out = sys.argv[3]

class Gene(object):
    def __init__(self):
        self.contig = ""
        self.start = None
        self.stop = None
        self.strand = "."

def write_file(outfile, list_lines):
    outfile_handler = open(outfile, 'w')

    for line in list_lines:
        outfile_handler.write(line + '\n')
    print("write file " + outfile + " done!")

list_diff_gid = list()

for line in open(file_csv):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue
    if line.startswith("not") == True:
        break

    splits = line.split('\t')

    gid = ""
    # g4401.t1 > g4401
    if splits[0].find('AUGUSTUS') != -1:
        gid = splits[0].split(' ')[0].split('.')[0]
    # mRNA.scaffold_10-snap.228 > scaffold_10-snap.228
    elif splits[0].find('SNAP') != -1:
        gid = splits[0].split(' ')[0][5:]
    # 795_t > 795_g
    elif splits[0].find('GeneMark.hmm') != -1:
        gid = splits[0].split(' ')[0].split('_')[0] + '_g'

    list_diff_gid.append(gid)

list_extract = list()

# record_mark: if True, control the mRNA exon CDS stop start BE extracted
record_mark = False
for line in open(file_gff3):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue

    splits = line.split('\t')
       
    if splits[2] == "gene":
        record_mark = False 
        gene_id = ""
        if splits[1] == "AUGUSTUS":
            gene_id = splits[8].split(';')[0].split('=')[1]
        elif splits[1] == "SNAP":
            gene_id = splits[8].split(';')[0].split('=')[1]
        elif splits[1] == "GeneMark.hmm":
            gene_id = splits[8].split(';')[0].split('=')[1].split('_')[0] + "_g"

        if gene_id in list_diff_gid:
            record_mark = True
            # record this line
            # gene
            list_extract.append(line)
            continue

    if record_mark == True:
        list_extract.append(line)

print "process extract meaningful proteins done!"

write_file(file_out, list_extract)
