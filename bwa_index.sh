#!/usr/bin/env bash
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#bwa alignment
rm -rf $dir/bwa
mkdir $dir/bwa
cd "$dir/bwa"
bwa index -a bwtsw $dir/$1 -p $2
