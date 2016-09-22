#!/usr/bin/python

import sys
import argparse
"""
making the training set
"""

def prepare_training_set(file_matched, file_for_training, filter_rule):
    results_handler = open(file_for_training, 'w')
    for line in open(file_matched):
        line = line.strip()
        if line == "":
            continue

        splits = line.split("\t")

        if len(splits) == 10:
            if splits[9] != '-':
                    groups = list()
                    for part in splits[8].strip(';').split('; '):
                        #groups.append(part.strip().split(' ')[0] + ' ' + '\"' + part.strip().split(' ')[1] +'\"')
                        groups.append(part.strip().split(' ')[0] + ' ' + part.strip().split(' ')[1])

                    new_line = '\t'.join(splits[:8]) + '\t' + '; '.join(groups)
                    results_handler.write(new_line + '\n')
            continue 
        # if not transposon
        if splits[11] == '-':
            if filter_rule == 'loose':
                print 'l'
                # criteria is loose: just demaing either of homology proteins mapping or transcripts alignment
                if splits[9] != '-' or splits[10] != '-':
                    groups = list()
                    for part in splits[8].strip(';').split('; '):
                        #groups.append(part.strip().split(' ')[0] + ' ' + '\"' + part.strip().split(' ')[1] +'\"')
                        groups.append(part.strip().split(' ')[0] + ' ' + part.strip().split(' ')[1])

                    new_line = '\t'.join(splits[:8]) + '\t' + '; '.join(groups)
                    results_handler.write(new_line + '\n')
            if filter_rule == 'strict':
                # criteria is strict: demading both of homology proteins mapping and transcripts alignment
                if splits[9] != '-' and splits[10] != '-':
                    groups = list()
                    for part in splits[8].strip(';').split('; '):
                        #groups.append(part.strip().split(' ')[0] + ' ' + '\"' + part.strip().split(' ')[1] +'\"')
                        groups.append(part.strip().split(' ')[0] + ' ' + part.strip().split(' ')[1])

                    new_line = '\t'.join(splits[:8]) + '\t' + '; '.join(groups)
                    results_handler.write(new_line + '\n')
            if filter_rule == 'bias_RNASeq':
                # criteria bias on RNA-Seq
                if splits[9] != '-':
                    # gff like
                    groups = list()
                    for part in splits[8].strip(';').split('; '):
                        #groups.append(part.strip().split(' ')[0] + ' ' + '\"' + part.strip().split(' ')[1] +'\"')
                        groups.append(part.strip().split(' ')[0] + ' ' + part.strip().split(' ')[1])

                    new_line = '\t'.join(splits[:8]) + '\t' + '; '.join(groups)
                    results_handler.write(new_line + '\n')

    results_handler.close()

def prepare_training_set_gli(file_for_training, file_for_glimmerhmm):
    gene_id_fix = ""
    minus_mark = False
    list_genes = list()
    results_handler = open(file_for_glimmerhmm, 'w')
    for line in open(file_for_training):
        line = line.strip()
        if line == "":
            continue

        splits = line.split('\t')
        gene_id = splits[8].split(';')[0].split(' ')[1]

        new_line = ""

        # checking for mark
        if gene_id == gene_id_fix:
            if splits[2] == "CDS":
               if splits[6] == "+":
                    new_line = splits[0] + ' ' + splits[3] + ' ' + splits[4]
                    results_handler.write(new_line + '\n')
               else:
                    minus_mark = True
                    new_line = splits[0] + ' ' + splits[4] + ' ' + splits[3]
                    list_genes.append(new_line)
        else:
            # a new line for a new gene
            gene_id_fix = gene_id
            if minus_mark:
                list_genes.reverse()
                for line in list_genes:
                    results_handler.write(line + '\n')

                minus_mark = False
                list_genes = list()

            results_handler.write('\n')

    results_handler.close()

#
# gene_id_fix = ""
# minus_mark = False
# list_genes = list()
# results_handler = open("d:\\filter_for_train_snap.txt", 'w')
# for line in open("d:\\test.csv"):
#     line = line.strip()
#     if line == "":
#         continue
#
#     splits = line.split('\t')
#     gene_id = splits[8].split(';')[0].split(' ')[1]
#
#     # checking for mark
#     if splits[11] == '-':
#         if splits[9] != '-' or splits[10] != '-':
#             if gene_id == gene_id_fix:
#                 if splits[2] == "CDS":
#                     if splits[6] == "+":
#                         new_line = splits[3] + ' ' + splits[4] + ' ' + gene_id
#                         list_genes.append(new_line)
#                     else:
#                         minus_mark = True
#                         new_line = splits[4] + ' ' + splits[3] + ' ' + gene_id
#                         list_genes.append(new_line)
#             else:
#                 # a new line for a new gene
#                 gene_id_fix = gene_id
#
#                 if  len(list_genes) == 0:
#                     continue
#                 elif len(list_genes) == 1:
#                     results_handler.write('>' + splits[0] + '\n')
#                     results_handler.write("Esngl " + list_genes[0] + '\n')
#                 elif len(list_genes) == 2:
#                     results_handler.write('>' + splits[0] + '\n')
#                     if minus_mark == True:
#                         results_handler.write("Einit " + list_genes[1] + '\n')
#                         results_handler.write("Eterm " + list_genes[0] + '\n')
#                     else:
#                         results_handler.write("Einit " + list_genes[0] + '\n')
#                         results_handler.write("Eterm " + list_genes[1] + '\n')
#                 elif len(list_genes) > 2:
#                     results_handler.write('>' + splits[0] + '\n')
#                     if minus_mark == True:
#                         list_genes.reverse()
#
#                         results_handler.write("Einit " + list_genes[0] + '\n')
#                         for line in list_genes[1:len(list_genes) - 1]:
#                             results_handler.write("Exon " + line + '\n')
#                         results_handler.write("Eterm " + list_genes[len(list_genes) - 1] + '\n')
#                     else:
#                         results_handler.write("Einit " + list_genes[0] + '\n')
#                         for line in list_genes[1:len(list_genes) - 1]:
#                             results_handler.write("Exon " + line + '\n')
#                         results_handler.write("Eterm " + list_genes[len(list_genes) - 1] + '\n')
#
#                 list_genes = list()
#                 minus_mark = False
#
# results_handler.close()
#
# w_h_genome = open("d:\\4742_genome.dna", 'w')
#
# list_zff_seqid = list()
# list_glm_seqid = list()
# for line in open("d:\\filter_for_train_snap.txt"):
#     line = line.strip()
#     if line == "":
#         continue
#
#     if line.startswith('>'):
#         if line[1:] not in list_zff_seqid:
#             list_zff_seqid.append(line[1:])
#         if line[1:] not in list_glm_seqid:
#             list_glm_seqid.append(line[1:])
#
# print str(len(list_zff_seqid)), str(len(list_glm_seqid))
#
# mark_in = False
# for line in open("d:\\NJAU-4742.fna"):
#     line = line.strip()
#     if line == "":
#         continue
#
#     if line.startswith('>'):
#         if line[1:] in list_zff_seqid:
#             w_h_genome.write(line + '\n')
#             mark_in = True
#         else:
#             mark_in = False
#     else:
#         if mark_in:
#             w_h_genome.write(line + '\n')

if __name__ == '__main__':

    parser_arg = argparse.ArgumentParser(description='Using matched_set.gtf to prepare training set')
    parser_arg.add_argument('--matched_set', dest='matched_set', help='Path to matched set file')
    parser_arg.add_argument('--filter_rule', dest='filter_rule', help='three filtering rule: loose, strict, bias_RNASeq')
    parser_arg.add_argument('--training_set', dest='training_set', help='Path to the output file: traning set')
    parser_arg.add_argument('--training_set_gli', dest='training_set_gli', help='Path to the output file: glimmer training set')

    args = parser_arg.parse_args()

    if not len(sys.argv) > 1:
        print parser_arg.parse_args(['-h'])
        sys.exit(255)

    file_matched = args.matched_set
    filter_rule = args.filter_rule
    file_for_training = args.training_set
    file_for_glimmerhmm = args.training_set_gli

    prepare_training_set(file_matched, file_for_training, filter_rule)
    if file_for_glimmerhmm != None:
        prepare_training_set_gli(file_for_training, file_for_glimmerhmm)

