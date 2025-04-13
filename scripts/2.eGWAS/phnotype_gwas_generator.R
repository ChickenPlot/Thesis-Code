#setwd("C:/Users/Simo/Desktop/phenotypes")
setwd("C:/Users/Simo/genomics/damianoRNAseq-20042023/output/GWAS/phenotypes")

library(tidyverse)
library(viridis)
library(patchwork)



norm_and_plot <- function(x, main = 'Normalization plot') {
  
  shap <- shapiro.test(x)
  qqnorm(x,
         main=main,
         ylab='Empirical quantiles',
         xlab='Theoretical quantiles')
  qqline(x,col='red')
  mtext(paste('Shap. test p-value: ', signif(shap$p.value, digit=4)), side = 3)
}

fd <- read.csv('fd_pheno.csv')
fd <- fd[c('accession_id', 'phenotype_value')]
colnames(fd) <- c('accession_id', 'AT4G35900')

fdp <- read.csv('fdp_pheno.csv')
fdp <- fdp[c('accession_id', 'phenotype_value')]
colnames(fdp) <- c('accession_id', 'AT2G17770')

eel <- read.csv('eel_pheno.csv')
eel <- eel[c('accession_id', 'phenotype_value')]
colnames(eel) <- c('accession_id', 'AT2G41070')

areb <- read.csv('areb3_pheno.csv')
areb <- areb[c('accession_id', 'phenotype_value')]
colnames(areb) <- c('accession_id', 'AT3G56850')



##################################################################
####################### PCA ######################################


left_join(fd, fdp, by = join_by(accession_id == accession_id)) %>%
  left_join(., eel,by = join_by(accession_id == accession_id)) %>%
    left_join(., areb,by = join_by(accession_id == accession_id)) %>%
        #mutate(., IID = accession_id) %>%
          rename(FID = accession_id) -> tmp_2

tmp <- tmp_2

tmp[,c("FID", "AT4G35900")] <- tmp[,c("AT4G35900", "FID")]
colnames(tmp) <- c('FID', 'IID', 'AT2G17770', 'AT2G41070', 'AT3G56850','AT4G35900' )

#tmp_df <- cbind(tmp$FID, tmp$IID, tmp$AT2G17770)

#typeof(tmp_df)
#write_tsv(tmp, file = "GWAS_phenos.pheno")
#write.csv(tmp, file = "GWAS_phenos.csv", row.names = F)

prcomp(t(tmp[,3:6])) -> all_cov


summary(all_cov)$importance %>% t() %>% data.frame() %>%
  ggplot(., aes(x = factor(rownames(.), 
                           levels = rownames(.)), 
                y = Cumulative.Proportion, fill = Cumulative.Proportion )) + geom_bar(stat = 'identity') + scale_fill_viridis(discrete = F) + geom_hline(yintercept = .9, linetype = 2) + theme_minimal() -> pca_plot_cumulat

summary(all_cov)$importance %>% t() %>% data.frame() %>%
  ggplot(., aes(x = factor(rownames(.), 
                           levels = rownames(.)), 
                y = Proportion.of.Variance
                , fill = Proportion.of.Variance
  )) + geom_bar(stat = 'identity') + scale_fill_viridis(discrete = F, direction = -1) + geom_hline(yintercept = .9, linetype = 2) + theme_minimal() -> pca_plot_varexpl

pca_plot_varexpl + pca_plot_cumulat


covs <- all_cov$rotation[,c(1:3)] %>% data.frame()
#write_tsv(covs, 'GWAS_covs.cov')

#################################################################
################ NORMALIZATION TRIES ############################
set.seed(123)

tmp <- tmp_2

### boxCox ####
bc1 <- MASS::boxcox((tmp[,2]+.1)~1, lambda = seq(-5,5,0.5))
mylambda1 <- bc1$x[which.max(bc1$y)]
mylambda1

bc2 <- MASS::boxcox((tmp[,3]+.1)~1, lambda = seq(-5,5,0.5))
mylambda2 <- bc2$x[which.max(bc2$y)]
mylambda2

bc3 <- MASS::boxcox((tmp[,4]+.1)~1, lambda = seq(-5,5,0.5))
mylambda3 <- bc3$x[which.max(bc3$y)]
mylambda3

bc4 <- MASS::boxcox((tmp[,5]+.1)~1, lambda = seq(-5,5,0.5))
mylambda4 <- bc4$x[which.max(bc4$y)]
mylambda4
# for all 0.3535354

tmp <- tmp %>%
  #filter(xpr > 0) %>%
  mutate(AT4G35900 = ((tmp[,2]^mylambda1)-1)/mylambda1) %>%
  mutate(AT3G56850 = ((tmp[,5]^mylambda4)-1)/mylambda2) %>%
  mutate(AT2G17770 = ((tmp[,3]^mylambda2)-1)/mylambda3) %>%
  mutate(AT2G41070 = ((tmp[,4]^mylambda3)-1)/mylambda4) %>%
  as.data.frame()

pivot_longer(tmp, cols = starts_with('AT'),
                names_to = c('gene.id'), 
              values_to = 'expression_value') %>%

ggplot(., aes(x = expression_value)) + 
  geom_histogram() +
  facet_wrap(~gene.id, scales = "free") +
  theme_minimal()


par(mfrow=c(1,4))

#j<-subset(sinergon_df, (Concentrazione == 0.0), select=c(Indice_proliferazione))
#j<- as.numeric(unlist(j))
#j<- log1p(j)
shap <- shapiro.test(tmp[,2])
qqnorm(tmp[,2],
       main='Q-Q plot: Indice di proliferazione (controllo)',
       ylab='Empirical quantiles',
       xlab='Theoretical quantiles')
qqline(tmp[,2],col='red')
mtext(paste('Shap. test p-value: ', signif(shap$p.value, digit=4)), side = 3)

shap <- shapiro.test(tmp[,3])
qqnorm(tmp[,3],
       main='Q-Q plot: Indice di proliferazione (controllo)',
       ylab='Empirical quantiles',
       xlab='Theoretical quantiles')
qqline(tmp[,3],col='red')
mtext(paste('Shap. test p-value: ', signif(shap$p.value, digit=4)), side = 3)

shap <- shapiro.test(tmp[,4])
qqnorm(tmp[,4],
       main='Q-Q plot: Indice di proliferazione (controllo)',
       ylab='Empirical quantiles',
       xlab='Theoretical quantiles')
qqline(tmp[,4],col='red')
mtext(paste('Shap. test p-value: ', signif(shap$p.value, digit=4)), side = 3)

shap <- shapiro.test(tmp[,5])
qqnorm(tmp[,5],
       main='Q-Q plot: Indice di proliferazione (controllo)',
       ylab='Empirical quantiles',
       xlab='Theoretical quantiles')
qqline(tmp[,5],col='red')
mtext(paste('Shap. test p-value: ', signif(shap$p.value, digit=4)), side = 3)




#### SQRT #####

tmp <- tmp_2
tmp[,2:5] <- sqrt(tmp[,2:5])

par(mfrow = c(1,4))

norm_and_plot(tmp[,2], main = 'sqrt')
norm_and_plot(tmp[,3],  main = 'sqrt')
norm_and_plot(tmp[,4], main = 'sqrt')
norm_and_plot(tmp[,5], main = 'sqrt')


#### LOG ####

tmp <- tmp_2

tmp[,2:5] <- log1p(tmp[,2:5])

par(mfrow = c(1,4))

norm_and_plot(tmp[,2], main = 'log')
norm_and_plot(tmp[,3], main = 'log')
norm_and_plot(tmp[,4], main = 'log')
norm_and_plot(tmp[,5], main = 'log')




# 1: box
# 2: box
# 3: box
# 4: box

tmp <- tmp_2


par(mfrow = c(1,4))
#tmp[, 5] <- log1p(tmp[, 5])
tmp[,2] <- ((tmp[,2]^mylambda1)-1)/mylambda1
tmp[,3] <- ((tmp[,3]^mylambda2)-1)/mylambda2
tmp[,4] <- ((tmp[,4]^mylambda3)-1)/mylambda3
tmp[,5] <- ((tmp[,5]^mylambda4)-1)/mylambda3



norm_and_plot(tmp[,2], main = 'box')
norm_and_plot(tmp[,3], main = 'box')
norm_and_plot(tmp[,4], main = 'box')
norm_and_plot(tmp[,5], main = 'box')

tmp <- tmp %>% mutate(IID = FID)
tmp[,c("IID", "AT4G35900")] <- tmp[,c("AT4G35900", "IID")]
colnames(tmp) <- c('FID', 'IID', 'AT2G17770', 'AT2G41070', 'AT3G56850','AT4G35900' )

write_tsv(tmp, 'GWAS_pheNORM.pheno')
#########################################################
###  ##################################


read_tsv('GWAS_pheNORM.pheno') -> gwas_phenorm
read.csv("1001genomes_admixture_table_final.csv") -> admix_file


admix_file[admix_file$relict == "Y",]$id -> relict_id
admix_file[admix_file$country == "ESP",]$id -> spain_id

# subset from gwas_phenorm  only iberian and relict accessions

length(relict_id)
length(spain_id)
access_to_subset <-  c(relict_id, spain_id) %>% unique #%>% length()

access_to_subset%>% data.frame()

areb[areb$accession_id %in% access_to_subset,] -> areb_subset
eel[eel$accession_id %in% access_to_subset,] -> eel_subset
fd[fd$accession_id %in% access_to_subset,] -> fd_subset
fdp[fdp$accession_id %in% access_to_subset,] -> fdp_subset

length(access_to_subset) # 183

# all 166, we are missing  some data
length(areb_subset$accession_id)
length(eel_subset$accession_id)
length(fd_subset$accession_id)
length(fdp_subset$accession_id)

#testino -- tutto ok
table(areb_subset$accession_id %in% eel_subset$accession_id)
table(areb_subset$accession_id %in% fd_subset$accession_id)
table(areb_subset$accession_id %in% fdp_subset$accession_id)
table(eel_subset$accession_id %in% fdp_subset$accession_id)
table(eel_subset$accession_id %in% fd_subset$accession_id)
table(fd_subset$accession_id %in% fdp_subset$accession_id)






cbind(areb_subset$accession_id,
      areb_subset$accession_id, 
      fdp_subset$AT2G17770, 
      eel_subset$AT2G41070, 
      areb_subset$AT3G56850, 
      fd_subset$AT4G35900) %>% 
  data.frame() %>% 
  rename(FID = X1, 
         IID = X2, 
         AT2G17770 = X3, 
         AT2G41070 = X4,
         AT3G56850 = X5, 
         AT4G35900 = X6) -> raw_phenos




# all norma by default
norm_and_plot(raw_phenos$AT2G17770)
norm_and_plot(raw_phenos$AT2G41070)
norm_and_plot(raw_phenos$AT3G56850)
norm_and_plot(raw_phenos$AT4G35900)



write_tsv(raw_phenos, "gwas_Spain_Relict_data.pheno")































































