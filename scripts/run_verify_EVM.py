#!/usr/bin/python

import sys
import re
"""
"""

if len(sys.argv) != 7:
    print "usage: " + sys.argv[0] + " evm.gff3 augustus.gff3 snap.gff3 genemark.gff3 assembly.bed exonerate.out"
    sys.exit(255)

file_evm = sys.argv[1]
file_augs = sys.argv[2]
file_snap = sys.argv[3]
file_gmak = sys.argv[4]
file_assm = sys.argv[5]
file_exon = sys.argv[6]

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
list_evm = list()
for line in open(file_evm):
    line = line.strip()
    if line == "":
        continue

    if line.find("gene") == -1:
        continue

    splits = line.split('\t')
    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3])
    gene.stop = int(splits[4])
    gene.strand = splits[6]

    list_evm.append(gene)
print "process EVM done, evm has genes: " + str(len(list_evm))

list_augs = list()
for line in open(file_augs):
    line = line.strip()
    if line == "" or line.startswith('#'):
        continue

    if line.find("gene") == -1:
        continue
    
    splits = line.split('\t')
    gene_id = splits[8].split(';')[0].split('=')[1]
    contig = splits[0]
    start = int(splits[3])
    stop = int(splits[4])
    strand = splits[6]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    # check the gene wether in list_evm
    for gene in list_evm:
        if contig == gene.contig:
            # different gene location
            if start > gene.stop or stop < gene.start:
                pass
            else:
                found = True
                break
    
    if found == False:
        list_augs.append(gene_id)
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
write_file("augustus_diff.gff3", wlist_augs)

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
    contig = splits[0]
    start = int(splits[3])
    stop = int(splits[4])
    strand = splits[6]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    for gene in list_evm:
        if contig == gene.contig:
            if start > gene.stop or stop < gene.start:
                pass
            else:
                found = True
                break
    if found == False:
        list_snap.append(gene_id)
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
write_file("snap_diff.gff3", wlist_snap)

list_gmak = list()
for line in open(file_gmak):
    line = line.strip()
    if line == "":
        continue

    if line.find("gene") == -1:
        continue

    splits = line.split('\t')
    gene_id = splits[8].split(';')[0].split('=')[1].split('_')[0] + "_"
    contig = splits[0]
    start = int(splits[3])
    stop = int(splits[4])
    strand = splits[6]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    for gene in list_evm:
        if contig == gene.contig:
            if start > gene.stop or stop < gene.start:
                pass
            else:
                found = True
                break
    if found == False:
        list_gmak.append(gene_id)
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
write_file("genemark_diff.gff3", wlist_gmak)

list_assm = list()
for line in open(file_assm):
    line = line.strip()
    if line == "":
        continue

    splits = line.split('\t')

    contig = splits[0]
    start = int(splits[1])
    stop = int(splits[2])
    strand = splits[5]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    for gene in list_evm:
        if contig == gene.contig:
            if start > gene.stop or stop < gene.start:
                continue
            else:
                found = True
                break
    if found == False:
        list_assm.append(line)
print "process PASA Assembly done, different count from evm: " + str(len(list_assm))
write_file("pasa_diff.gff3", list_assm)

list_exon = list()
for line in open(file_exon):
    line = line.strip()
    if line == "" or not line.startswith("vulgar"):
        continue

    splits = line.split()

    contig = splits[5]
    start = int(splits[6]) if int(splits[6]) < int(splits[7]) else int(splits[7])
    stop = int(splits[7]) if int(splits[6]) < int(splits[7]) else int(splits[6])
    strand = splits[8]

    found = False
    # don't care about the stranded, it is considered by the EVM,
    # we can judge them through the genome brwose.
    for gene in list_evm:
        if contig == gene.contig:
            if start > gene.stop or stop < gene.start:
                continue
            else:
                found = True
                break
    if found == False:
        list_exon.append(line)
print "process Exonerate done, different count from evm: " + str(len(list_exon))
write_file("exone_diff.gff3", list_exon)

