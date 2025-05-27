#!/bin/bash 

CHRS=($(seq 1 19) "X" "Y")

mkdir mus38_split

for chr in "${CHRS[@]}"
do
  echo "extracting chr${chr}"
	faidx --regex ".*chr${chr}.*" Mmus38.geve.nt_v1.fa -o "mus38_split/chr${chr}.fa"
done
