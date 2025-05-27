#!/bin/bash

# CHRS=($(seq 1 19) "X" "Y")
CHRS=("Y")
for chr in "${CHRS[@]}"
do
  echo "extracting .fa from chr${chr}.gtf"
	agat_sp_extract_sequences.pl -g ./repeatmasker/chr${chr}.gtf -t mRNA --merge 0\
    -f ../single_cell/mouse/GRCm39.primary_assembly.genome.fa -o ./repeatmasker/chr${chr}.fa
done

