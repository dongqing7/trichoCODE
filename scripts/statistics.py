#!/usr/bin/python
from __future__ import division
import sys

def count_matched(f_in):
    dict_gene_support = dict()
    count_transcripts = 0
    count_exonerate = 0
    count_transposon = 0

    for line in open(f_in):
        line = line.strip()

        splits = line.split('\t')

        groups = splits[8].split('; ')
        gene_id = groups[0].split(' ')[1].strip('"')

        if not dict_gene_support.has_key(gene_id):
            dict_gene_support[gene_id] = 0

            if splits[9] != '-':
                count_transcripts += 1
            if splits[10] != '-':
                count_exonerate += 1
            if splits[11] != '-':
                count_transposon += 1

    print 'Statistics of training set'
    print 'Supported by transcript data:   ' + str(count_transcripts)
    print 'Supported by homology proteins: ' + str(count_exonerate)
    print 'Found Transposon:               ' + str(count_transposon)

def count_gff(f_gff):
    count_genes = 0
    count_exons = 0
    # count for genes
    exon_end = 0
    next_exon_start = 0

    # count for intron
    count_introns = 0
    length_total_introns = 0

    for line in open(f_gff):
        line = line.strip()
        if line.startswith('#') or line == '':
            continue

        splits = line.split('\t')
        if splits[2] == 'gene':
            count_genes += 1
            exon_end = 0

        # count for exons
        # total num, ave exons
        if splits[2] == 'exon':
            count_exons += 1

            # count for introns
            # total num, total length, ave length
            if exon_end == 0:
                if splits[6] == '+':
                    exon_end = int(splits[4])
                elif splits[6] == '-':
                    exon_end = int(splits[3])
            else:
                tmp_exon_end = 0
                if splits[6] == '+':
                    next_exon_start = int(splits[3])
                    tmp_exon_end = int(splits[4])
                elif splits[6] == '-':
                    next_exon_start = int(splits[4])
                    tmp_exon_end = int(splits[3])
                length_total_introns += abs(next_exon_start - exon_end)
                exon_end = tmp_exon_end
                count_introns += 1

    print 'Statistics of gff'
    print 'Total genes:            ' + str(count_genes)
    print 'Exons'
    print '  Total number:         ' + str(count_exons)
    print '  Ave Exons per gene:   ' + str(round(count_exons/count_genes, 1))
    print 'Introns'
    print '  Total number:         ' + str(count_introns)
    print '  Total length:         ' + str(length_total_introns)
    print '  Ave Introns per gene: ' + str(round(count_introns/count_genes, 1))
    print '  Ave length:           ' + str(round(length_total_introns/count_introns, 2))

def count_pasa(f_gff):
    count_5_utr = 0 
    count_3_utr = 0 
    count_alter_splicing = 0
    count_exon_correction = 0

    for line in open(f_gff):
        line = line.strip()
        if line.startswith('#PROT'):
            continue

        if line.startswith('#'):
            # count PASA update information
            if line.find('single gene model update') != -1:
                count_exon_correction += 1
            elif line.find('alt-splice') != -1:
                count_alter_splicing += 1
        else:
            # count for UTR
            if line.find('three_prime_UTR') != -1:
                count_3_utr += 1
            elif line.find('five_prime_UTR') != -1:
                count_5_utr += 1
    
    print 'Statistics of PASA update'
    print '5 prime UTR:              ' + str(count_5_utr)
    print '3 prime UTR:              ' + str(count_3_utr)
    print 'alternative splicing:     ' + str(count_alter_splicing)
    print 'exon boundary correction: ' + str(count_exon_correction)
     

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage: %s matched.gtf type(match|gff|pasa)" % (sys.argv[0])
        exit(255)

    f_in = sys.argv[1]
    t = sys.argv[2]

    if t == 'match':    count_matched(f_in)
    if t == 'gff':  count_gff(f_in)
    if t == 'pasa': count_pasa(f_in)
