#!/bin/bash

# A major issue is that AGAT, the programme that we use to 
# turn GTF files into FA files, require that each sequence has 
# an ID and be marked as mRNA (another cheap hack that I used 
# to get things to work). Consequently, you will see some cheap 
# hacks when I create these columns

# CHRS=($(seq 1 19) "X" "Y")
CHRS=("Y")
for chr in "${CHRS[@]}"
do
  echo "extracting .fa from chr${chr}.gtf"
	agat_sp_extract_sequences.pl -g ./repeatmasker/chr${chr}.gtf -t mRNA --merge 0\
    -f ../single_cell/mouse/GRCm39.primary_assembly.genome.fa -o ./repeatmasker/chr${chr}.fa
done

# another verison
cd /PHShome/jn22/siyi_summer2023/erv/00repeatmasker
CHRS=`ls gtfs_for_fasta_needed_by_salmon/*.gtf`
for chr in $CHRS
do
  echo "extracting .fa from ${chr}"
  echo "output file is ${chr%.*}.fa"
	agat_sp_extract_sequences.pl -g ${chr} -t mRNA --merge 0\
    -f ./GRCm39.primary_assembly.genome.fa -o ${chr%.*}.fa
done



