filespath = "C:/Users/Simo/genomics/damianoRNAseq-20042023/output/GWAS/phenotypes"
setwd(filespath)
getwd()

library(tidyverse)
library(readxl)
library(viridis)
library(ComplexHeatmap)

access_info <- read.csv("1001access.csv") #table from 1001
admixture_table <- read.csv("1001genomes_admixture_table_final.csv")
eel_haplos <- read.csv("../EEL_SNP_Haplos.csv")

areb_pheno <- read.csv("areb3_pheno.csv")
eel_pheno <- read.csv("eel_pheno.csv")#data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
fd_pheno <- read.csv("fd_pheno.csv")
fdp_pheno <-read.csv("fdp_pheno.csv")

gwas_phenos <- read.csv("GWAS_phenos.csv")
gwas_phenorm <- read_tsv("GWAS_pheNORM.pheno")
gwas_phenorm <- gwas_phenorm[,-2]



admixture_table$id
admixture_table$group %>% unique # 10 groups 

gwas_phenos <-  gwas_phenos[gwas_phenos$FID %in% admixture_table$id,] # 665/726
gwas_phenos<-  gwas_phenos[order(gwas_phenos$FID),]

admix_table_sub <- admixture_table[admixture_table$id %in% gwas_phenos$FID, ] 
admix_table_sub <- admix_table_sub[order(admix_table_sub$id), c(1,3,4,6,7,11,12)]

admix_table_sub[402,]
gwas_phenos[402,]


view(cbind(admix_table_sub, gwas_phenos))

access_expr_table <- cbind(admix_table_sub, gwas_phenos)
access_expr_table <- access_expr_table[,-7]

# now same with normalized
gwas_phenorm <-  gwas_phenorm[gwas_phenorm$FID %in% admixture_table$id,] # 665/726
gwas_phenorm<-  gwas_phenorm[order(gwas_phenorm$FID),]
access_NORMexpr_table <- cbind(admix_table_sub, gwas_phenorm)
access_NORMexpr_table <- access_NORMexpr_table[,-7]


# plots
#histo of distribution

access_expr_table %>% ggplot(., aes(x = group, fill = group))+
  geom_histogram(stat = "count")+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))

access_expr_table %>% ggplot(., aes(x = country, fill = group))+
  geom_histogram(stat = "count")+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = .8, 
                                   vjust = 1))


# raw pheno
access_expr_table %>% 
  ggplot(., aes(x = group, y = AT3G56850, fill = group))+
  geom_boxplot()+
  theme_minimal()+ggtitle("AREB3 expression across haplogroups")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))


access_expr_table %>% 
  ggplot(., aes(x = group, y = AT2G41070, fill = group))+
  geom_boxplot()+
  theme_minimal() +ggtitle("EEL expression across haplogroups")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))


access_expr_table %>% 
  ggplot(., aes(x = group, y = AT4G35900, fill = group))+
  geom_boxplot()+
  theme_minimal() +ggtitle("FD expression across haplogroups")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))



access_expr_table %>% 
  ggplot(., aes(x = group, y = AT2G17770, fill = group))+
  geom_boxplot()+
  theme_minimal() +ggtitle("FDP expression across haplogroups")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))


# normalized pheno
access_NORMexpr_table %>% 
  ggplot(., aes(x = group, y = AT3G56850, fill = group))+
  geom_boxplot()+
  theme_minimal()+ggtitle("AREB3 expression across haplogroups-normalized")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))


access_NORMexpr_table %>% 
  ggplot(., aes(x = group, y = AT2G41070, fill = group))+
  geom_boxplot()+
  theme_minimal() +ggtitle("EEL expression across haplogroups-normalized")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))


access_NORMexpr_table %>% 
  ggplot(., aes(x = group, y = AT4G35900, fill = group))+
  geom_boxplot()+
  theme_minimal() +ggtitle("FD expression across haplogroups-normalized")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))



access_NORMexpr_table %>% 
  ggplot(., aes(x = group, y = AT2G17770, fill = group))+
  geom_boxplot()+
  theme_minimal() +
  ggtitle("FDP expression across haplogroups-normalized")+
  #scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 20, 
                                   hjust = .8, 
                                   vjust = 1))


# show groups geographically

library("rnaturalearth")
library("rnaturalearthdata")

# one point has NA coords
access_expr_table <- access_expr_table[complete.cases(access_expr_table),]  

world <- ne_countries(scale = "medium", returnclass = "sf")

lat_lim <-  c(min(access_expr_table$latitude),max(access_expr_table$latitude))
long_lim <- c(min(access_expr_table$longitude),max(access_expr_table$longitude)) 

ggplot(data = world) +
  geom_sf() + 
  coord_sf(xlim = long_lim, ylim = lat_lim, expand = T)+
  geom_point(data = access_expr_table, aes(x = longitude , y =latitude, color = group )) + 
  theme_minimal()+scale_color_viridis_d()

""

#############################################################################################

eel_haplos <- eel_haplos[,-c(1,7)]

table(gwas_phenos$FID %in% eel_haplos$SampleID)
sub_phenos <- gwas_phenos[gwas_phenos$FID %in% eel_haplos$SampleID,]
sub_phenos$Haplotype <- eel_haplos$Haplotype






##############################################################################################




















































""