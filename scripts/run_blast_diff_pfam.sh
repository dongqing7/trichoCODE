#!/bin/sh

if [ $# != 6 ]
then
    echo "USAGE: $0 species_diff.gff3 genome species_diff pfam_db trichCode_dir evm_dir"
    exit 255;
fi

echo "extract protein sequences ..."
$6/EvmUtils/gff3_file_to_proteins.pl $1 $2 > $3.faa

echo "processing blastping ..."
#blastp -query $3.faa -db ~/db/nr/nr -outfmt 5 -max_target_seqs 1 -out $3.xml -num_threads 2 -evalue 1e-7
rpsblast -query $3.faa -db $4 -out $3.pfam -outfmt 5 -max_target_seqs 3 -num_threads 4 -evalue 1e-10 

echo "process rpblast results ..."
#$5/scripts/process_blast_result.py $3.pfam b
$5/scripts/process_rpsblast_pfam.py $3.pfam

echo "extract fit gene (gff3)"
$5/scripts/run_extract_fit.py $1 $3.csv $3_fit.gff3

echo "translate into faa"
$6/EvmUtils/gff3_file_to_proteins.pl $3_fit.gff3 $2 > $3_fit.faa

