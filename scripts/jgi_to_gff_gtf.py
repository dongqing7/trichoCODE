#!/usr/bin/env python

__author__ = 'biomobil'

import sys
#from collections import Counter


class Features(object):
    def __init__(self):
        self.feature = ''
        self.start = -1
        self.end = -1
        self.score = '.'
        self.frame = '.'

class Gene(object):
    def __init__(self, geneid):
        self.contig = ''
        self.source = ''
        self.strand = ''
        self.start = ''
        self.stop = ''
        self.transid = geneid
        self.geneid = geneid
        self.features = list()

dict_genes = dict()
list_genes_order = list()
def read_jgi_gtf(f_jgi):
    ph_geneid = ''
    ct_features = 0
    ct_genes = 0
    for line in open(f_jgi):
        line = line.strip()
        if line == '' or line.find('_codon') != -1 :
            continue

        # partial gene, don't consider this for now

        # complete gene (start_codon, stop_codon)
        # first circle
        splits = line.split('\t')
        group_info = splits[8].split(';')
        ph_geneid = group_info[1].strip().split(' ')[1]
        if dict_genes.has_key(ph_geneid):
            gene = dict_genes[ph_geneid]
            gene.start = gene.start if int(gene.start) < int(splits[3]) else splits[3]
            gene.stop = gene.stop if int(gene.stop) > int(splits[4]) else splits[4]
            feature = Features()
            feature.feature = splits[2]
            feature.start = splits[3]
            feature.end = splits[4]
            feature.frame = splits[7]
            if gene.strand == '-':
                gene.features.insert(ct_features, feature)
                ct_features += 1
            else:
                gene.features.append(feature)
        else:
            ct_genes += 1
            gene = Gene(ph_geneid)
            gene.contig = splits[0]
            gene.source = splits[1]
            gene.strand = splits[6]
            gene.start = splits[3]
            gene.stop = splits[4]
            feature = Features()
            feature.feature = splits[2]
            feature.start = splits[3]
            feature.end = splits[4]
            feature.frame = splits[7]
            dict_genes[ph_geneid] = gene
            gene.features.append(feature)
            list_genes_order.append(ph_geneid)
            ct_features += 1

   # print 'genes count is: ' + str(ct_genes)

def write_output(format):
    for geneid in list_genes_order:
        if format == 'gtf':
            gene = dict_genes[geneid]
            for feature in gene.features:
                print gene.contig + '\t' + gene.source + '\t' + feature.feature + '\t' + feature.start + '\t' + feature.end + '\t' + '.' + '\t' + gene.strand + '\t' + feature.frame + '\t' + 'gene_id \"' + gene.geneid + '\"; ' + 'transcript_id \"' + gene.transid + '\"'
        elif format == 'gff':
            gene = dict_genes[geneid]
            gene.geneid = gene.geneid.strip('"')
            gene.transid = gene.transid.strip('"')
            # gene
            print gene.contig + '\t' + gene.source + '\t' + 'gene' + '\t' + gene.start + '\t' + gene.stop + '\t.\t' + gene.strand + '\t.\t' + 'ID=' + gene.geneid
            print gene.contig + '\t' + gene.source + '\t' + 'mRNA' + '\t' + gene.start + '\t' + gene.stop + '\t.\t' + gene.strand + '\t.\t' + 'ID=model.' + gene.transid + ';Parent=' + gene.geneid
            i = 1
            j = 1
            for feature in gene.features:
                if feature.feature == 'exon':
                    print gene.contig + '\t' + gene.source + '\t' + feature.feature + '\t' + feature.start + '\t' + feature.end + '\t' + '.' + '\t' + gene.strand + '\t' + feature.frame + '\t' + 'ID=exon' + str(i) + '.' + gene.geneid + ';Parent=' + gene.geneid
                    i += 1
                if feature.feature == 'CDS':
                    print gene.contig + '\t' + gene.source + '\t' + feature.feature + '\t' + feature.start + '\t' + feature.end + '\t' + '.' + '\t' + gene.strand + '\t' + feature.frame + '\t' + 'ID=cds' + str(j) + '.' + gene.geneid + ';Parent=' + gene.geneid
                    j += 1

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print "usage: " + sys.argv[0] + " jgi.gff outtype(gff|gtf)"
        sys.exit(255)
    f_jgi = sys.argv[1]
    format = sys.argv[2]
    read_jgi_gtf(f_jgi)
    write_output(format)

