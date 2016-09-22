#!/usr/bin/python
import sys
import argparse

"""
mark the genes/exons with RNASeq
"""

dict_cuff = dict()
dict_cuff_exon = dict()
exon_support_count = 0
gene_support_count = 0

class Gene(object):
    def __init__(self):
        self.contig = ""
        self.start = None
        self.stop = None
        self.strand = "."

class Exon(object):
    def __init__(self):
        self.contig = ""
        self.start = None
        self.stop = None
        self.strand = "."

def read_genes(file_cuff):
    global dict_cuff

    # alignment from STAR & cufflinks cuffmerge
    # scaffold_1  Cufflinks   exon    6052    7419    .   +   .   gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "6"; oId "CUFF.1.1"; tss_id "TSS1";
    for line in open(file_cuff):
        line = line.strip()
        if line == "":
            continue

        splits = line.split('\t')
        gene_id = splits[8].split(';')[0].split(' ')[1]
        if dict_cuff.has_key(gene_id):
            # update value of start & stop
            dict_cuff[gene_id].start = dict_cuff[gene_id].start if dict_cuff[gene_id].start < int(splits[3]) else int(splits[3])
            dict_cuff[gene_id].stop = dict_cuff[gene_id].stop if dict_cuff[gene_id].stop > int(splits[4]) else int(splits[4])

        else:
            gene = Gene()
            gene.contig = splits[0]
            gene.strand = splits[6]
            gene.start = int(splits[3])
            gene.stop = int(splits[4])

            dict_cuff[gene_id] = gene
    print "process genes from cufflinks done"

def read_exons(file_cuff):
    global dict_cuff_exon
    
    for line in open(file_cuff):
        line = line.strip()
        if line == "":
            continue

        splits = line.split('\t')
        gene_id = splits[8].split('; ')[0].split(' ')[1]
        exon_order = splits[8].split('; ')[2].split(' ')[1]
        exon_id = gene_id + exon_order
        exon = Exon()
        exon.contig = splits[0]
        exon.start = splits[3]
        exon.stop = splits[4]
        exon.strand = splits[6]

        dict_cuff_exon[exon_id] = exon
    print "process exons from cufflinks done"

def find_match(file_gm, file_out_training):
    global exon_support_count
    global gene_support_count

    dict_ges = dict()
    results_handler = open(file_out_training, 'w')
    list_lines = list()
    tmp_geneid = ""
    for line in open(file_gm):
        line = line.strip()
        if line == "" or line.find('CDS') != -1:
            continue
        list_lines.append(line)

        splits = line.split('\t')
        gene_id = splits[8].split(';')[0].split(' ')[1]

        mark_exon = False
        for exon_item in dict_cuff_exon.itervalues():
            if splits[0] == exon_item.contig:
                if splits[6] == exon_item.strand:
                    if splits[4] < exon_item.start or splits[3] > exon_item.stop:
                        pass
                    else:
                        exon_support_count += 1
                        mark_exon = True
                        break

        if dict_ges.has_key(gene_id):
            # update value of start & stop
            dict_ges[gene_id].start = dict_ges[gene_id].start if dict_ges[gene_id].start < int(splits[3]) else int(splits[3])
            dict_ges[gene_id].stop = dict_ges[gene_id].stop if dict_ges[gene_id].stop > int(splits[4]) else int(splits[4])
        else:
            if tmp_geneid == "":
                gene = Gene()
                gene.contig = splits[0]
                gene.strand = splits[6]
                gene.start = int(splits[3])
                gene.stop = int(splits[4])
                dict_ges[gene_id] = gene
                tmp_geneid = gene_id
                continue

            start = dict_ges[tmp_geneid].start
            stop = dict_ges[tmp_geneid].stop
            strand = dict_ges[tmp_geneid].strand
            contig = dict_ges[tmp_geneid].contig
            
            mark_trans_found = False
            # find match with cufflinks&merge alignment
            for gene_item in dict_cuff.itervalues():
                if contig == gene_item.contig:
                    if strand == gene_item.strand or gene_item.strand == '.':
                        if stop < gene_item.start or start > gene_item.stop:
                            pass
                        else:
                            gene_support_count += 1
                            mark_trans_found = True
                            break

            for line in list_lines:
                if mark_trans_found:
                    line += "\ttranscript"
                else:
                    line += "\t-"
                results_handler.write(line + '\n')
            list_lines = list()
            
            gene = Gene()
            gene.contig = splits[0]
            gene.strand = splits[6]
            gene.start = int(splits[3])
            gene.stop = int(splits[4])
            dict_ges[gene_id] = gene
            
            tmp_geneid = gene_id
    results_handler.close()
    print "gene supported: " + str(gene_support_count)
    print "exon supported: " + str(exon_support_count)
    print "process done"

if __name__ == '__main__':

    parser_arg = argparse.ArgumentParser(description='find gene models with matched evidence')
    parser_arg.add_argument('--transcriptalg', action='store', type=str, dest='gtf_rnaseq', help='Path to gff file of transcript alignment')
    parser_arg.add_argument('--gtf', action='store', dest='gtf_in', help='Path to genemark-et output file')
    parser_arg.add_argument('--out', action='store', dest='gtf_out', help='Name of output file')
    parser_arg.add_argument('--version', action='version', version='1.0')

    args = parser_arg.parse_args()

    if not len(sys.argv) > 1:
        print parser_arg.parse_args(['-h'])
        sys.exit(255)

    file_gtf_rnaseq = args.gtf_rnaseq
    file_gtf_in = args.gtf_in
    file_gtf_out = args.gtf_out

    read_genes(file_gtf_rnaseq)
    read_exons(file_gtf_rnaseq)
    
    find_match(file_gtf_in, file_gtf_out)

