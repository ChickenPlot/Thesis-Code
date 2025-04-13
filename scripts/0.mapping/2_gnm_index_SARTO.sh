#!/bin/bash


# cd genomics/arabidopsis_references/transcriptome/
refs_dir=../../arabidopsis_reference/transcriptome

pigz -d ${refs_dir}/TAIR10_GFF3_genes.gtf.gz

ref_genome=${refs_dir}/TAIR10_chr_all_nms.fas 
ref_gff=${refs_dir}/TAIR10_GFF3_genes.gtf
out_folder=${refs_dir}/indexed_genome

STAR --runMode genomeGenerate --genomeDir "${out_folder}" --genomeFastaFiles "${ref_genome}" --sjdbGTFfile "${ref_gff}"  --sjdbOverhang 100 --runThreadN 8 --genomeSAindexNbases 12

pigz ${refs_dir}/TAIR10_GFF3_genes.gtf
echo "Done generating index"