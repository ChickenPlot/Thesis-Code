#!/bin/bash

out_dir=../output
# mkdir ../output/assembly_genomeref
assembly_dir=${out_dir}/assembly_genomeref

refs_dir=../../arabidopsis_reference/transcriptome


#######################################
### Extract counts from STAR output ###
#######################################

echo "Extracting counts from STAR"

for i in ${assembly_dir}/*ReadsPerGene.out.tab
do echo $i

cut -f1,2 $i | grep -v "^N_[anmu]" > ${assembly_dir}/`basename $i .out.tab`_unstranded.txt #with STAR data in this htseq-format you can now use HTseq import in DESeq2

done

echo "STAR generated counts processed"



####################################
#### Make count table from HTSeq ###
####################################

echo "Generating counts with HTSeq"

#pigz -d ${refs_dir}/TAIR10_GFF3_genes.gff.gz
pigz -d ${refs_dir}/TAIR10_GFF3_genes.gtf.gz

for b in ${assembly_dir}/*.sortedByCoord.out.bam
do

bname=$(basename ${b%%.out.bam})

echo "HTseq-count of ${bname%%.sortedByCoord.out.bam}" 

##  can be used with .gtf or .gff file   ##
#htseq-count -f bam -r pos -s no -t exon -i gene_id -m union --add-chromosome-info -n 8 ${b} ${refs_dir}/TAIR10_GFF3_genes.gff > ${b}.counts.txt
htseq-count -f bam -r pos -s no -t exon -i gene_id -m union --add-chromosome-info -n 8 ${b} ${refs_dir}/TAIR10_GFF3_genes.gtf > ${bname}.HTSeq_counts.txt

done

#pigz ${refs_dir}/TAIR10_GFF3_genes.gff.gz
pigz ${refs_dir}/TAIR10_GFF3_genes.gtf



