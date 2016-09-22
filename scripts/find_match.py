#!/usr/bin/python
import sys
import argparse

"""
the subset of GeneMark-ES predictions that have structures 
consistent with transcript alignments(PASA) and Exonerate(GeneWise)
predictions are extracted to further supplement the Training Set.

Also predicted proteins found to exhibit homology to known repeats
(ascertained by usnig TransposonPSI) are excluded from this collection.
"""

# read from TransposonPSI protein
dict_psi = dict()

# read from PASA alignment
dict_pasa = dict()

dict_cuff = dict()
# read from Exonerate homology
list_exr = list()

class Gene(object):
    def __init__(self):
        self.contig = ""
        self.start = None
        self.stop = None
        self.strand = "."


def read_psi(file_psi):
    global dict_psi
    for line in open(file_psi):
        line = line.strip()
        if line == "":
            continue

        splits = line.split('\t')
        psi_id = splits[8].split(';')[0].split('=')[1]

        gene = Gene()
        gene.contig = splits[0]
        gene.start = int(splits[3])
        gene.stop = int(splits[4])
        gene.strand = splits[6]

        dict_psi[psi_id] = gene
    print "process psi done"

def read_pasa(file_pasa):
    global dict_pasa
    if file_pasa != "":
        for line in open(file_pasa):
            line = line.strip()
            if line == "":
                continue

            splits = line.split('\t')
            # gff3
            gene_id = splits[8].split('|')[1].split(';')[0]
            # gtf
            # gene_id = splits[8].split(';')[0].split(' ')[1]

            if splits[2] != "gene":
                continue

            if not dict_pasa.has_key(gene_id):
                gene = Gene()
                gene.contig = splits[0]
                gene.start = int(splits[3])
                gene.stop = int(splits[4])
                gene.strand = splits[6]

                dict_pasa[gene_id] = gene
        print "process pasa done"

def read_cuff(file_cuff):
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
            gene.stop =   int(splits[4])

            dict_cuff[gene_id] = gene
    print "process cufflinks&merge done"

def read_exonerate(file_exon):
    global list_exr
    for line in open(file_exon):
        line = line.strip()
        if line == "":
            continue

        if line.startswith("vulgar"):
            splits = line.split(' ')
            gene = Gene()
            gene.contig = splits[5]
            gene.start = int(splits[6]) if int(splits[6]) < int(splits[7]) else int(splits[7])
            gene.stop = int(splits[7]) if int(splits[6]) < int(splits[7]) else int(splits[6])
            gene.strand = splits[8]

            list_exr.append(gene)
    print "process Exonerate done"

def find_match(type_transcript, file_gm, file_out_training):
    # build GeneMark-ES prediction dictionary.
    dict_ges = dict()
    results_handler = open(file_out_training, 'w')
    list_lines = list()
    # this tmp used to remember last gene
    tmp_geneid = ''
    for line in open(file_gm):
        line = line.strip()
        if line == "": # or line.find('_codon') != -1:
            continue

        splits = line.split('\t')
        gene_id = splits[8].split('; ')[0].split(' ')[1]

        if dict_ges.has_key(gene_id):
            # update value of start & stop
            dict_ges[gene_id].start = dict_ges[gene_id].start if dict_ges[gene_id].start < int(splits[3]) else int(splits[3])
            dict_ges[gene_id].stop = dict_ges[gene_id].stop if dict_ges[gene_id].stop > int(splits[4]) else int(splits[4])
            list_lines.append(line)
        else:
            # when read the 1st gene
            if tmp_geneid == '':
                gene = Gene()
                gene.contig = splits[0]
                gene.strand = splits[6]
                gene.start = int(splits[3])
                gene.stop = int(splits[4])
                dict_ges[gene_id] = gene
                # remember this gene
                tmp_geneid = gene_id
                list_lines.append(line)
                continue

            # when arrived here, means the last gene is compeletely read done
            start = dict_ges[tmp_geneid].start
            stop = dict_ges[tmp_geneid].stop
            strand = dict_ges[tmp_geneid].strand
            contig = dict_ges[tmp_geneid].contig

            mark_trans_found = False
            if type_transcript == 'pasa':
                # find match with PASA alignment

                for gene_item in dict_pasa.itervalues():
                    if contig == gene_item.contig:
                        if strand == gene_item.strand:
                            if stop < gene_item.start or start > gene_item.stop:
                                pass
                            else:
                                mark_trans_found = True
                                break
            elif type_transcript == 'cuff':
                # find match with cufflinks&merge alignment
                for gene_item in dict_cuff.itervalues():
                    if contig == gene_item.contig:
                        if strand == gene_item.strand or gene_item.strand == '.':
                            if stop < gene_item.start or start > gene_item.stop:
                                pass
                            else:
                                mark_trans_found = True
                                break

            # find matches with Exonerate homology
            mark_exon_found = False
            for gene_item in list_exr:
                if contig == gene_item.contig:
                    if strand == gene_item.strand:
                        if stop < gene_item.start or start > gene_item.stop:
                            pass
                        else:
                            mark_exon_found = True
                            break

            # find matches with TransposonPSI
            mark_psi_found = False
            for gene_item in dict_psi.itervalues():
                if contig == gene_item.contig:
                    if stop < gene_item.start or start > gene_item.stop:
                        pass
                    else:
                        mark_psi_found = True
                        break

            for line in list_lines:
                if mark_trans_found:
                    line += "\ttranscript"
                else:
                    line += "\t-"
                if mark_exon_found:
                    line += "\tExonerate"
                else:
                    line += "\t-"
                if mark_psi_found:
                    line += "\tTransposonPSI"
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
    print "process done"

if __name__ == '__main__':

    parser_arg = argparse.ArgumentParser(description='find gene models with matched evidence')
    parser_arg.add_argument('--transposon', action='store', type=str, dest='gff_psi', help='Path to transposonPSI output file')
    parser_arg.add_argument('--transcriptalg', action='store', type=str, dest='alg_gff', help='Path to gff file of transcript alignment')
    parser_arg.add_argument('--transcripttype', action='store', type=str, dest='alg_tool', help='Type of transcript alignment: \'cuff\' or \'pasa\'')
    parser_arg.add_argument('--exonerate', action='store', dest='vul_exonerate', help='Path to exonerate output file')
    parser_arg.add_argument('--genemark-et', action='store', dest='gtf_genemark', help='Path to genemark-et output file')
    parser_arg.add_argument('--out', action='store', dest='gtf_matched_set', help='Name of output file')
    parser_arg.add_argument('--version', action='version', version='1.0')

    args = parser_arg.parse_args()

    if not len(sys.argv) > 1:
        print parser_arg.parse_args(['-h'])
        sys.exit(255)

    file_psi = args.gff_psi
    file_alg_gff = args.alg_gff
    file_exon = args.vul_exonerate
    file_gm = args.gtf_genemark
    file_out_training = args.gtf_matched_set
    type_transcript = args.alg_tool
    
    read_psi(file_psi)
    read_exonerate(file_exon)

    if type_transcript != None and type_transcript.strip() != '':
        if type_transcript == 'cuff':
            read_cuff(file_alg_gff)
        elif type_transcript == 'pasa':
            read_pasa(file_alg_gff)
        else:
            print sys.stderr, '--transcripttype could be recognised'
            exit(1)
    else:
        type_transcript = ''
    
    find_match(type_transcript, file_gm, file_out_training)

