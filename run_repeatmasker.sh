#!/bin/bash

# running repeatmasker on the mouse chr automatically
# bsub -q interactive -Is -R "rusage[mem=32000]" -n 24 repeatmasker1
# parallel = threads / 2 = 24//2

chrDir=/PHShome/jn22/Downloads/chrs
outDir=/PHShome/jn22/Downloads/repeatmasker

for chr in "$@"
do
    mkdir $outDir/${chr}
    RepeatMasker -engine hmmer -parallel 12 -species "Mus musculus"\
      -gccalc -dir $outDir/${chr} -gff $chrDir/${chr}.fa
done

#  egrep -v "Simple|Satellite" my_data.out > filtered.out
