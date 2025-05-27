#!/bin/bash

CHRS=($(seq 1 19) "X" "Y")
chr="Y"
mkdir ./blastdb/rm_mus39_chr${chr}
makeblastdb -in repeatmasker/chr${chr}.fa -title rm_mus39_chr${chr} \
      -dbtype nucl -parse_seqids -out ./blastdb/rm_mus39_chr${chr}/rm_mus39_chr${chr}
