#!/bin/bash

# modifying the erv stuff 
# need gene_id & gene_name for gtf, column 1 

awk -F"\t" '$3==BED_feature' # I forgot this command :(
sed -i 's/BED_feature/gene/g' erv.gff3
sed -i 's/ID=\([A-Za-z0-9]\+\);/gene_id "\1"; /g' erv.gff3
sed -i 's/Name=\([A-Za-z0-9\:]\+\)/gene_name "\1"/g' erv.gff3
awk -F"\t" '$3=exon' erv.gff3


