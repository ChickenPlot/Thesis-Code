setwd("C:/Users/Simo/genomics/damianoRNAseq-20042023/scripts/lyrata_eel")

library(tidyverse)
library(ape)
library(poppr)
library(Biostrings)

#chr_4 <- read.FASTA("Arabidopsis_lyrata.v.1.0.dna.chromosome.4.fa")


#chr4 <- DNAbin2genind(chr_4)

ch4 <- readDNAStringSet("Arabidopsis_lyrata.v.1.0.dna.chromosome.4.fa")
ch4
subseq(ch4,start = (20026222-2000), end = (20027810+2000)) %>% DNAStringSet(.) %>%  writeXStringSet(.,"Lyrata_EEL.fa")


###########################

library(ggtree)
treefile <- ape::read.tree("EEL_phylotree.nwk")
nodelabels(treefile)

ggtree(treefile)
ggtree(treefile)+theme_tree2()+layout_dendrogram()

#
newtree <- drop.tip(treefile, treefile$tip.label[1])
ggtree(newtree)+geom_tiplab2()


#select_labels <- 
#selected_labels <- sort(sample(0:length(newtree$tip.label), length(newtree$tip.label)*0.6))
ggtree(newtree)+theme_tree2()+layout_circular()

old <-  newtree$tip.label
new <- lapply(strsplit(newtree$tip.label, split = "\\|"), "[", 4)

newtree$tip.label[match(old, newtree$tip.label)] <- new
ggtree(newtree)+theme_tree2()+layout_dendrogram()+
  geom_tiplab()





tree2 <- ape::read.tree("EEL_1001_pseudogenome_TRIMMED.mft.treefile")
ggtree(tree2)

tree3 <- ape::read.tree("EEL_1001_pseudogenome_TRIMMED.mft.contree")
ggtree(tree3)




# renamed tree

newlabels <- lapply(strsplit(newtree$tip.label, split = "\\|"), "[", 4)

newtree2 <- newtree 
newtree2$tip.label <- newlabels[match(newtree$tip.label,newtree$tip.label)] 
par(mfrow=c(1,2))  
ggtree(newtree)

ggtree(newtree2)
# find acc.n for the nib and subset them
main_access <- read_file("spain_relict_access_n.txt")
strsplit(main_access, "\r\n") %>% unlist -> main_access
main_access <- unlist(main_access)

ggtree(newtree2)+geom_tiplab2(aes(subset= (newtree2$tip.label%in%main_access)==T ))
