#!/bin/bash

#awk '/^>/ {out = substr($1, 2) ".fasta"; print > out} !/^>/ {print >> out}' blastextended2/BLASTResults/LD_region_PacBios_extended2.fasta

filespath=blastextended2/BLASTResults/splitted_regions/*

for LDfile in ${filespath}
do
outname=$(basename ${LDfile})

echo "mapping ${LDfile}"
minimap2 -ax asm5 -t 4 --eqx ../TAIR_LD_region.fna ${LDfile} | samtools sort -O BAM - > minimap_results2/${outname}_aligned.bam && samtools index  minimap_results2/${outname}_aligned.bam

done