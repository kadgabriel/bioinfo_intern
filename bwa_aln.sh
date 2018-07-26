#!/usr/bin/env bash
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#bwa alignment
cd $dir/bwa
pwd
bwa aln $1 $dir/reads/$2_read1.fastq > $2_aln_sa1.sai
bwa aln $1 $dir/reads/$2_read2.fastq > $2_aln_sa2.sai

bwa sampe $1 $dir/bwa/$2_aln_sa1.sai $dir/bwa/$2_aln_sa2.sai $dir/reads/$2_read1.fastq $dir/reads/$2_read2.fastq > $dir/bwa/aligned_pairs_$2.sam