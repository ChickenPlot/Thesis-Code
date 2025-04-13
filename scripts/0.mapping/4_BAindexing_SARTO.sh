#!/bin/bash

# define paths

out_dir=../output
assembly_dir=${out_dir}/assembly_genomeref

#####################################################
# indexing bam files sorted by coordinate
echo "start indexing BAMs sorted by coord"

for bam_file in ${assembly_dir}/*.sortedByCoord.out.bam
do 

echo 'Processing file'

outname=$(basename ${bam_file})

samtools index ${bam_file} "${assembly_dir}/${outname}.bai"

echo 'File processed'
echo

done

echo "All file processed"


#############################################
# indexing bam files sorted by Name

echo "start indexing BAMs sorted by Name"

for bam_file in ${assembly_dir}/*.sortedByName.out.bam
do 

echo 'Processing file'

outname=$(basename ${bam_file})

samtools index ${bam_file} "${assembly_dir}/${outname}.bai"

echo 'File processed'
echo

done

echo "All file processed"
echo 
echo "Indexing Done"
