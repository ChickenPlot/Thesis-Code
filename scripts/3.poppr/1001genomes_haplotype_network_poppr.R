setwd("C:/Users/Simo/Desktop/new bea scripts")
library(tidyverse)
library(vcfR)
library(ape)
library(poppr)

# DATA
# Genomic data 1001 genomes
# vcf1001 <- read.PLINK(
#   "1001genomes/1001genomes_snp_only_ACGTN_filtered_pruned_chr2_genotypes.raw",
#   map.file = NULL,
#   quiet = FALSE,
#   chunkSize = 100)
# 
# saveRDS(vcf1001, file="1001genomes/1001genomes_snp_only_ACGTN_filtered_pruned_chr2_genotypes.rds")

#getwd()

vcf1001 <- readRDS("1001genomes_snp_only_ACGTN_filtered_pruned_chr2_genotypes.rds")
vcf1001

vcf1001 <- as.snpclone(vcf1001)
mlg.filter(vcf1001, threads = 1L) <- 0.25 #26 mlgs
vcf1001

## Info Accessions 1001 genomes
acc_info <- read.csv('1001access.csv')
#acc_info <- acc_info[c(1:3, 11)]

eel_haplos_info <- read.csv("EEL_SNP_Haplos.csv")
eel_haplos_info <- eel_haplos_info[,-1]

colnames(acc_info)

tmp__hap <- rep("Unknown", dim(acc_info)[1])

tmp__hap[which(acc_info$pk %in% eel_haplos_info$SampleID)] <- eel_haplos_info$Haplotype

table(tmp__hap)

acc_info$Haplotype <- tmp__hap

vcf1001@strata <- acc_info
vcf1001@strata

vcf1001_2 <- vcf1001
length(acc_info[,1])

# CALCULATE NETWORK & PLOT
# vcf1001.dist <- bitwise.dist(
#   vcf1001,
#   percent = TRUE,
#   mat = TRUE,
#   missing_match = TRUE,
#   scale_missing = FALSE,
#   euclidean = FALSE,
#   differences_only = FALSE,
#   threads = 4)
# 
# saveRDS(vcf1001.dist, file="1001genomes_snp_only_ACGTN_filtered_pruned_chr2_genotypes_distMX.rds")

vcf1001.dist <- readRDS("1001genomes_snp_only_ACGTN_filtered_pruned_chr2_genotypes_distMX.rds")

setPop(vcf1001) <- ~Haplotype
vcf1001.msn <- poppr.msn(vcf1001, vcf1001.dist,
                         palette = c("#FF4040", "#00008B", "#EEEEE0"),
                         gscale = TRUE)

plot_poppr_msn(vcf1001,
               vcf1001.msn,
               inds = "none",
               mlg = TRUE,
               gadj = 3,
               wscale = FALSE,
               nodescale = 5,
               palette = c("darkblue", "red","gray"),
               quantiles = FALSE,
               beforecut = FALSE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_with_graphopt)


#
setPop(vcf1001_2) <- ~country
vcf1001_2.msn <- poppr.msn(vcf1001_2, vcf1001.dist,
                         palette = c("#FF4040", "#00008B", "#EEEEE0"),
                         gscale = TRUE)

plot_poppr_msn(vcf1001_2,
               vcf1001_2.msn,
               inds = "none",
               mlg = TRUE,
               gadj = 3,
               wscale = FALSE,
               nodescale = 5,
               palette = c("darkblue", "red","gray"),
               quantiles = FALSE,
               beforecut = FALSE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_with_graphopt)

#vcf1001_2@strata$
