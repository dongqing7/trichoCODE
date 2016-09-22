#!/usr/bin/python

import sys
import re
"""
"""

if len(sys.argv) != 4:
    print "usage: " + sys.argv[0] + " augustus.gff3 snap.gff3 genemark.gff3"
    sys.exit(255)

file_augs = sys.argv[1]
file_snap = sys.argv[2]
file_gmak = sys.argv[3]

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

# read evm results to list_evm
list_diff = list()

list_augs = list()
for line in open(file_augs):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue

    if line.find("gene") == -1:
        continue
    
    splits = line.split('\t')
    gene_id = splits[8].split(';')[0].split('=')[1]
    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3])
    gene.stop = int(splits[4])
    gene.strand = splits[6]

    list_augs.append(gene_id)
    list_diff.append(gene)
print "process Augustus done, different count from evm: " + str(len(list_augs))

wlist_augs = list()
for line in open(file_augs):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue
    
    for gene_id in list_augs:
        if line.find(gene_id+';') != -1 or line.find(gene_id+'.') != -1:
            wlist_augs.append(line)
            break 
write_file("augustus_diff_2.gff3", wlist_augs)

list_snap = list()
for line in open(file_snap):
    line = line.strip()
    if line == "":
        continue

    if line.find("gene") == -1:
        continue
#   contig_10   SNAP    gene    90626   91842   .   -   .   ID=contig_10-snap.27;
    splits = line.split('\t')
    gene_id = splits[8].split(';')[0].split('=')[1]
    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3])
    gene.stop = int(splits[4])
    gene.strand = splits[6]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    for gene_diff in list_diff:
        if gene.contig == gene_diff.contig:
            if gene.start > gene_diff.stop or gene.stop < gene_diff.start:
                pass
            else:
                found = True
                break
    if found == False:
        list_snap.append(gene_id)
        list_diff.append(gene)
print "process SNAP done, different count from evm: " + str(len(list_snap))
wlist_snap = list()
for line in open(file_snap):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue
    
    for gene_id in list_snap:
        if line.find(gene_id+';') != -1 or line.find(gene_id+'.') != -1:
            wlist_snap.append(line)
            break
write_file("snap_diff_2.gff3", wlist_snap)

list_gmak = list()
for line in open(file_gmak):
    line = line.strip()
    if line == "":
        continue

    if line.find("gene") == -1:
        continue

    splits = line.split('\t')
    gene_id = splits[8].split(';')[0].split('=')[1].split('_')[0] + "_"
    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3])
    gene.stop = int(splits[4])
    gene.strand = splits[6]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    for gene_diff in list_diff:
        if gene.contig == gene_diff.contig:
            if gene.start > gene_diff.stop or gene.stop < gene_diff.start:
                pass
            else:
                found = True
                break
    if found == False:
        list_gmak.append(gene_id)
        list_diff.append(gene)
print "process GeneMark-ES done, different count from evm: " + str(len(list_gmak))
wlist_gmak = list()
for line in open(file_gmak):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue
    
    for gene_id in list_gmak:
        if line.find('='+gene_id) != -1: #or line.find(':'+gene_id) != -1:
            wlist_gmak.append(line)
            break
write_file("genemark_diff_2.gff3", wlist_gmak)

