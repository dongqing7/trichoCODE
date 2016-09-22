#!/usr/bin/python
import sys

"""
check EVM results consistency with nonRNA, repeats(repeatMasker, repeatModeler),
PSI, write conflicts to outfile
"""

if len(sys.argv) != 7:
    print "usage: " + sys.argv[0] + " rRNA tRNA rep psi evm.gff3 outfile"
    sys.exit(255)

file_rRNA = sys.argv[1]
file_tRNA = sys.argv[2]
#file_rep_mask = sys.argv[3]
file_rep_modl = sys.argv[3]
file_psi = sys.argv[4]
file_evm = sys.argv[5]
file_out = sys.argv[6]

class Gene(object):
    def __init__(self):
        self.contig = ""
        self.start = None
        self.stop = None
        self.strand = "."

# check rRNA
list_rRNA = list()
for line in open(file_rRNA):
    line = line.strip()
    if line.startswith('#') or line == "":
        continue

    splits = line.split()

    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3]) if int(splits[3]) < int(splits[4]) else int(splits[4])
    gene.stop = int(splits[4]) if int(splits[3]) < int(splits[4]) else int(splits[3])
    gene.strand = splits[6]

    list_rRNA.append(gene)

# check tRNA
list_tRNA = list()
for line in open(file_tRNA):
    if line.startswith('#') or line == "":
        continue
    line = line.strip().replace(' ', '')
    splits = line.split()
    
    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[2]) if int(splits[2]) < int(splits[3]) else int(splits[3])
    gene.stop = int(splits[3]) if int(splits[2]) < int(splits[3]) else int(splits[2])
    gene.strand = "+" if int(splits[2]) < int(splits[3]) else "-"

    list_tRNA.append(gene)

# check repeatModeler
list_repMod = list()
for line in open(file_rep_modl):
    if line.startswith('#') or line == "":
        continue
    line = line.strip()
    splits = line.split('\t')

    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3])
    gene.stop = int(splits[4])
    gene.strand = splits[6]

    list_repMod.append(gene)

# check psi
list_psi = list()
for line in open(file_psi):
    if line.startswith('#') or line == "":
        continue
    line = line.strip()
    splits = line.split('\t')

    gene = Gene()
    gene.contig = splits[0]
    gene.start = int(splits[3])
    gene.stop = int(splits[4])
    gene.strand = splits[6]

    list_psi.append(gene)

results_handler = open(file_out, 'w')
for line in open(file_evm):
    line = line.strip()
    if line.startswith('#') or line == "":
        continue
    if line.find("gene") == -1:
        continue

    splits = line.split('\t')
    contig = splits[0]
    start = int(splits[3])
    stop = int(splits[4])
    strand = splits[6]

    # check rRNA
    mark_rRNA = False
    for gene in list_rRNA:
        if gene.contig == contig:
            if start > gene.stop or stop < gene.start:
                pass
            else:
                mark_rRNA = True
                break

    # check tRNA
    mark_tRNA = False
    for gene in list_tRNA:
        if gene.contig == contig:
            if start > gene.stop or stop < gene.start:
                pass
            else:
                mark_tRNA = True
                break

    # check repeat Mod
    mark_repMod = False
    for gene in list_repMod:
        if gene.contig == contig:
            if start > gene.stop or stop < gene.start:
                pass
            else:
                mark_repMod = True
                break

    # check PSI
    mark_psi = False
    for gene in list_psi:
        if gene.contig == contig:
            if start > gene.stop or stop < gene.start:
                pass
            else:
                mark_psi = True
                break

    if mark_rRNA:
        line += "\trRNA"
    else:
        line += "\t-"

    if mark_tRNA:
        line += "\ttRNA"
    else:
        line += "\t-"
    
    if mark_repMod:
        line += "\tnew repeat"
    else:
        line += "\t-"
    
    if mark_psi:
        line += "\tTransposonPSI"
    else:
        line += "\t-"
    
    results_handler.write(line + '\n')

results_handler.close()

