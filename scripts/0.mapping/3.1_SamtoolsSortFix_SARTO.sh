#!/bin/bash

cd ../output/assembly_genomeref
#refs_dir=../../arabidopsis_reference/transcriptome
#bam_dir=../output/assembly_genomeref

echo "Sorting .bam files by NAME"
for bamfile in *.bam
do
    echo ${bamfile}
    # set a basename for output file
    bam_name=$(basename ${bamfile%%.sortedByCoord.out.bam})
    samtools sort -n ${bamfile} -@ 8 -o "${bam_name}.sortedByName.out.bam" 

done
echo 
echo "Done sorting"
cd ../../scripts