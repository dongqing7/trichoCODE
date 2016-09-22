#!/usr/bin/python
import sys

if len(sys.argv) != 4:
    print "usage: %s start_sj.tab threshold gmet|aug" % (sys.argv[0])
    sys.exit(255)

file_sj = sys.argv[1]
threshold = sys.argv[2]
mode = sys.argv[3]

list_records = list()
count = 0

# scaffold_10 1001995 1002284 2   2   0   3   0   48_12
for line in open(file_sj):
    line = line.strip()
    
    splits = line.split('\t')
    if int(splits[6]) < int(threshold):
        continue

    strand = '+' if splits[3] == '1' else '-'
    start = splits[1] if int(splits[1]) < int(splits[2]) else splits[2]
    stop = splits[1] if int(splits[1]) > int(splits[2]) else splits[2]
    if mode == 'gmet':
        new_line = splits[0] + '\t' + 'STAR240f1' + '\t' + 'intron' + '\t'  + \
                    start + '\t' + stop + '\t' + splits[6] + '\t' + strand + '\t.\t.'
    elif mode == 'aug':
        new_line = splits[0] + '\t' + 'b2h' + '\t' + 'intron' + '\t'  + \
                    start + '\t' + stop + '\t' + splits[6] + '\t' + strand + '\t.\tmult='+splits[6]+';pri=4;src=E'

    print new_line
