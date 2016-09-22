#!/usr/bin/python

import commands
import os,sys

if len(sys.argv) != 6:
    print "usage: "+sys.argv[0]+" weight gene_predictions.gff3 protein_alg transcript genome"
    exit(255)

os.system("rm evm.gff3")

weight = sys.argv[1]
gene_predictions = sys.argv[2]
protein_alg = sys.argv[3]
transcript = sys.argv[4]
genome = sys.argv[5]

cmd_part = "../EvmUtils/partition_EVM_inputs.pl --genome " + genome + " --gene_predictions " + gene_predictions + " --protein_alignments " + protein_alg +" --transcript_alignments " + transcript + " --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out" 

res = os.system(cmd_part)
if res != 0:
    print commands.getstatusoutput(cmd_part)
    exit(255)

cmd_comm = "../EvmUtils/write_EVM_commands.pl --genome " + genome + " --weights `pwd`/" + weight + " --gene_predictions " + gene_predictions + " --protein_alignments " + protein_alg+ " --transcript_alignments " + transcript + " --output_file_name evm.out  --partitions partitions_list.out >  commands.list"

res = os.system(cmd_comm)
if res != 0:
    print commands.getstatusoutput(cmd_comm)
    exit(255)

cmd_exe = "../EvmUtils/execute_EVM_commands.pl commands.list | tee run.log;../EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out;../EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome " + genome

res = os.system(cmd_exe)
if res != 0:
    print commands.getstatusoutput(cmd_exe)
    exit(255)

print "please run: \" echo scaffold_* > list;run_cat.py list \""
