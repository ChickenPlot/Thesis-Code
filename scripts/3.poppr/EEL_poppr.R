library(tidyverse)
library(ape)
library(poppr)

#


setwd('C:/Users/Simo/genomics/damianoRNAseq-20042023/scripts/phylo')
eel_mft_fasta <- read.FASTA("EEL_1001_pseudogenome_TRIMMED.mft")


setwd('C:/Users/Simo/genomics/damianoRNAseq-20042023/output/GWAS/phenotypes')

acc_info <- read.csv('1001access.csv')
acc_info <- acc_info[c(1:3, 11)]

eel_haplos_info <- read.csv("../EEL_SNP_Haplos.csv")
eel_haplos_info <- eel_haplos_info[,-1]


#haplo_info <- read.csv('../EEL_SNP_Haplos.csv')
#haplo_info[c(1:8)]


colnames(acc_info)


tmp__hap <- rep("Unknown", dim(acc_info)[1])

tmp__hap[which(acc_info$pk %in% eel_haplos_info$SampleID)] <- eel_haplos_info$Haplotype

table(tmp__hap)

acc_info$Haplotype <- tmp__hap
  

setwd('C:/Users/Simo/genomics/damianoRNAseq-20042023/scripts/phylo')

# WHOLE GENE
eel_mft <- DNAbin2genind(eel_mft_fasta)
rownames(eel_mft$tab) %>%
  gsub("MPI-GMI\\|Ath-1001-Genomes\\|pseudo-genome\\|","",.) %>% 
  gsub("\\|Chr2","", .) -> tmp

table(as.numeric(tmp) %in% acc_info$pk)



length(rownames(eel_mft$tab))
dim(acc_info)
#dim(haplo_info)

#############
eel_mft@strata <- acc_info

length(acc_info[,1])
#length(eel_mft@strata)

rownames(eel_mft$tab) <- acc_info$pk

eel_mft <- missingno(eel_mft, type = "loci", cutoff = 0.05, quiet = FALSE, freq = FALSE)

setPop(eel_mft) <- ~Haplotype

eel_mft <- as.genclone(eel_mft)


imsn()
###########################à

eel_mft_fasta <-  eel_mft_fasta[-1]

#eel_haplos_info$SampleID %>% duplicated() %>% table
eel_mft <- DNAbin2genind(eel_mft_fasta)
rownames(eel_mft$tab) %>%
  gsub("MPI-GMI\\|Ath-1001-Genomes\\|pseudo-genome\\|","",.) %>% 
  gsub("\\|Chr2","", .) -> tmp


eel_mft@strata <- eel_haplos_info

rownames(eel_mft$tab) <- eel_haplos_info$SampleID #eel_haplos_info$SampleID

eel_mft <- missingno(eel_mft, type = "loci", cutoff = 0.05, quiet = FALSE, freq = FALSE)

setPop(eel_mft) <- ~Haplotype

eel_mft <- as.genclone(eel_mft)
imsn()









################################################################
eel_mft_sub <- popsub(eel_mft, exclude = character(0))
eel_mft_dist <- diss.dist(eel_mft_sub, percent = FALSE, mat = FALSE)
# eel_mft_nomiss <- missingno(eel_mft, type = 'mean')
# eel_mft_dist <- edwards.dist(eel_mft_nomiss)
min_span_net <- poppr.msn(eel_mft_sub, eel_mft_dist, showplot = FALSE, include.ties = TRUE)

min_span_net

set.seed(69)
plot_poppr_msn(eel_mft,
               min_span_net,
               inds = "none",
               mlg = TRUE,
               gadj = 3,
               wscale = FALSE,
               nodescale = 10,
               palette = c("gray", "darkblue","red"),
               cutoff = 3,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_with_graphopt)

# INTRON
aln.i <- DNAbin2genind(aln.i)

df.info <- data.frame(id = as.character(rownames(aln.i$tab))) %>% 
  mutate(id = gsub("_C.*", "", id)) %>%
  mutate(id = gsub("\\|.*", "", id)) %>%
  as.data.frame()

df.info <- df.info %>%
  left_join(.,
            hpl.flc,
            by = c("id" = "ecotype_id_flc")) %>%
  left_join(.,
            hpl.fri,
            by = c("id" = "accession_id_fri")) %>%
  as.data.frame()

aln.i@strata <- df.info
aln.i@strata$paper.simple_flc <- if_else(is.na(aln.i@strata$paper.simple_flc), "not available", aln.i@strata$paper.simple_flc)

rownames(aln.i$tab) <- df.info$id

aln.i <- missingno(aln.i, type = "loci", cutoff = 0.05, quiet = FALSE, freq = FALSE)

setPop(aln.i) <- ~paper.simple_flc

aln.i <- as.genclone(aln.i)

flc.i.mlg <- data.frame(
  mlg = as.character(aln.i@mlg@mlg$original),
  id = rownames(aln.i@tab))

flc.mlg <- full_join(flc.a.mlg %>%
                       rename(mlg.a = "mlg"),
                     flc.i.mlg %>%
                       rename(mlg.i = "mlg") %>%
                       mutate(id = gsub("ref", "Arabidopsis_thaliana", id)),
                     by = "id") %>%
  select(c(id, mlg.a, mlg.i)) %>%
  as.data.frame()

imsn()

aln.i_sub <- popsub(aln.i, exclude = character(0))
aln.i_dist <- diss.dist(aln.i_sub, percent = FALSE, mat = FALSE)
min_span_neti <- poppr.msn(aln.i_sub, aln.i_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
plot_poppr_msn(aln.i,
               min_span_neti,
               inds = "none",
               mlg = TRUE,
               gadj = 3,
               wscale = FALSE,
               nodescale = 10,
               palette = c("gray", "darkblue", "red"),
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = TRUE,
               scale.leg = TRUE,
               layfun = igraph::layout_with_graphopt)

write.csv(flc.mlg, "./output/flc_fri_genetrees/flc_mlg_msn.csv")
