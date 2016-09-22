#!/usr/bin/python

import sys

if len(sys.argv) != 2:
    print "usage: %s splice_junc_list " % (sys.argv[0])
    sys.exit(255)

file_splict_junc_list = sys.argv[1]

class SJ:
    def __init__(self, chrom, start, end, strand, motif, anno, num_uniq, num_cros, overhang):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.motif = motif
        self.anno = anno
        self.num_uniq = num_uniq
        self.num_cros = num_cros
        self.overhang = overhang
    
    def __repr__(self):
        return self.chrom + '\t' + self.start + '\t' + self.end + '\t' + self.strand + '\t' + self.motif + '\t' + self.anno + '\t' + str(self.num_uniq) + '\t' + str(self.num_cros) + '\t' + self.overhang

dict_records = dict()
for line in open(file_splict_junc_list):
    line = line.strip()

    if len(dict_records) == 0:
        for line in open(line):
            line = line.strip()
            splits = line.split('\t')
            record = SJ(splits[0], splits[1], splits[2], splits[3], splits[4], splits[5], int(splits[6]), int(splits[7]), splits[8])
            dict_records[splits[0] + '_' + splits[1] + '_' + splits[2]] = record
            # for testing
#            if splits[0] == 'scaffold_26' and splits[1] == '356102':
#                print record

    else:
        for line in open(line):
            line = line.strip()
            splits = line.split('\t')
            location = splits[0] + '_' + splits[1] + '_' + splits[2]
            if dict_records.has_key(location):
                dict_records[location].num_uniq += int(splits[6])
                dict_records[location].num_cros += int(splits[7])
                dict_records[location].overhang += '_' + splits[8]
            else:
                record = SJ(splits[0], splits[1], splits[2], splits[3], splits[4], splits[5], int(splits[6]), int(splits[7]), splits[8])
                dict_records[location] = record

for record in sorted(dict_records):
    print dict_records[record]
