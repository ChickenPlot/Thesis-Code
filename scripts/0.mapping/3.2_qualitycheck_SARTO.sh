#!/bin/bash

###########################################
#### 1. fastqc ##### 
#####################

echo "RUNNING FastQC"
out_dir=../output
rawdata_dir=../X204SC22090037-Z01-F002/01.RawData

fastqc ${rawdata_dir}/*.fq.gz -o ${out_dir} -noextract

echo
echo "Done running FastQC"
echo

##########################################
#### 2. qualimap ####
#####################

echo "RUNNING QualiMap" 
mkdir ../output/assembly_genomeref/qualimap

bam_dir=../output/assembly_genomeref
refs_dir=../../arabidopsis_reference/transcriptome

echo "--unzipping .gtf"
pigz -d ${refs_dir}/TAIR10_GFF3_genes.gtf.gz

for bamfile in ${bam_dir}/*.sortedByName.out.bam
do
    # set a basename for output file
    qlmp_out=$(basename ${bamfile%%.sortedByName.out.bam})
    
    # run qualimap ## memory increased but still doesn't work
    # NB. Memory allocated: 36GB
    qualimap rnaseq -bam $bamfile -gtf ${refs_dir}/TAIR10_GFF3_genes.gtf -outdir ${bam_dir}/qualimap -outfile ${qlmp_out} -pe -s --java-mem-size=36000M

    mv rnaseq_qc_result.txt "${qlmp_out}_rnaseq_qc_result.txt"
#############################################################
################### rename txt file #########################
#############################################################
done

echo "--zipping .gtf back"
pigz ${refs_dir}/TAIR10_GFF3_genes.gtf

echo "Done running QualiMap"
echo

##########################################
#### 3. multiqc ####
#####################

echo "Running MultiQC"
multiqc -o ${bam_dir}/multiqc \
../output/*zip \
../output/assembly_genomeref/*Log.final.out \
../output/assembly_genomeref/qualimap/* \

echo "Done running MultiQC"

