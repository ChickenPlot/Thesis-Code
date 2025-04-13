#!/bin/bash

bampath=minimap_results2/*.fasta_aligned.bam 
bamlist=($(ls ${bampath} -1 | sort))

fastapath=blastextended2/BLASTResults/splitted_regions
fastalist=($(ls ${fastapath} -1 | sort))

for i in {0..17..1}
do
fastafile=${fastalist[i]}
bamfile=${bamlist[i]}
outname=$(basename ${fastafile})

echo ""
echo ""

echo "PROCESSING ${outname}"

echo ""

syri -c ${bamfile} -r ../TAIR_LD_region.fna -q blastextended2/BLASTResults/splitted_regions/${fastafile} -F B --dir syri_out --prefix "${outname}_LD_syri."

done