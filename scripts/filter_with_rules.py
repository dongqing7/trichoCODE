#!/usr/bin/python
import sys

class Gene(object):
    def __init__(self):
        self.gid = 0
        self.start = 0
        self.stop = 0
        self.exons = list()
        self.introns = list()

class Exon(object):
    def __init__(self):
        self.start = 0
        self.stop = 0

class Intron(object):
    def __init__(self):
        self.splice5 = ''
        self.splice3 = ''
        self.length = 0

def filter_file(file_gff3, max_intron, min_cds):
    gene_list = list()

    exon_end = 0
    cds_length = 0
    intron_bigger = 0

    for line in open(file_gff3):
        line = line.strip()
        if line == '':
            continue
        
        # record the gene
        gene_list.append(line)

        splits = line.split('\t')
        if splits[2] == 'gene':
            gene = Gene()
            gene.gid = splits[8].split(';')[0].split('=')[1]
            
            # for the first gene, then go to the next line
            if cds_length == 0 and intron_bigger == 0:
                continue

            if cds_length/3 < min_cds or intron_bigger > max_intron:
                #print '########' + gene_list[0]
                #print cds_length, intron_bigger
                pass
            else:
                for l in gene_list[0:(len(gene_list)-1)]:
                    print l
                print ''

            # prepare for next gene
            gene_list = list()
            # add this gene line
            gene_list.append(line)
            cds_length = 0
            intron_bigger =0
            exon_end = 0
            
        # find the biggest intron in a gene
        if splits[2] == 'exon':
            if splits[6] == '-':
                if exon_end == 0:
                    exon_end = int(splits[3])
                else:
                    intron = exon_end - int(splits[4])
                    exon_end = int(splits[3])
                    if intron_bigger < intron:
                        intron_bigger = intron
            else:
                if exon_end == 0:
                    exon_end = int(splits[4])
                else:
                    intron = int(splits[3]) - exon_end
                    exon_end = int(splits[4])
                    if intron_bigger < intron:
                        intron_bigger = intron
        # accumulate cds length 
        if splits[2] == 'CDS':
            cds_length += int(splits[4]) - int(splits[3])


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "usage: %s file_gff3 intron_max aa_min " % (sys.argv[1])
        exit(255)

    infile = sys.argv[1]
    max_intron = int(sys.argv[2])
    min_aa = int(sys.argv[3])

    filter_file(infile, max_intron, min_aa)

