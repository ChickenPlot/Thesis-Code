#!/bin/bash
#conda activate bioTools
# define paths

#out_dir=../output
assembly_dir=assemblies
#mkdir assemblies/linear
#mkdir assemblies/linear/reheaded

#####################################################
# indexing bam files sorted by coordinate
#echo "start faidindexing fastas"

#for assembly_file in ${assembly_dir}/*.fna
#do 

#outname=$(basename ${assembly_file})
#echo 'Processing file'

#outname=$(basename ${bam_file})

#echo 'linearizing'
#seqkit seq -w 0 ${assembly_file} > "${assembly_file}_linear.fna"

#done


for assembly_file in ${assembly_dir}/*_linear.fna
do

echo 'reheading FASTA'
sed -E '/^>/s/.*chromosome:\s+([1-5]).*/>chromosome_\1/' ${assembly_file} > "${assembly_file}_reheaded.fna"

done


for assembly_file in ${assembly_dir}/*_reheaded.fna
do

echo 'indexing'
samtools faidx ${assembly_file} 

echo 'extrancting LD block'
samtools faidx ${assembly_file} chromosome_2:17126759-17134769 > "${assembly_file}_EXCTRACTED.fna"

echo 'File processed'

done

echo "All file processed"


#sed -E '/^>/s/.*chromosome\s+([1-5]).*/>chromosome_\1/' GCA_028009825.2_Col-CC_genomic.fna > Col0_reheaded.fasta
