#!/bin/sh

if [ $# != 3 ]
then
    echo "USAGE: $0 maker2zff genome_fas tag"
    exit 255;
fi

$1 training_set.gff3;grep '^>' genome.ann | tr -d '>' > genome.keep; fasta_sort.pl genome.keep < $2 > genome.dna; fathom genome.ann genome.dna -validate; fathom genome.ann genome.dna -categorize 1000; fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log;rm params -rf; mkdir params;cd params/; forge ../export.ann ../export.dna; cd ..; hmm-assembler.pl $3 params/ > $3.hmm; snap $3.hmm genome.dna > snap_$3.zff

