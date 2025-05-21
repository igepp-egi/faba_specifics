testmyY <- read.table("http://zzlab.net/GAPIT/data/mdp_traits.txt", head = TRUE)
testmyGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
testmyGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
# install.packages("BiocManager")
# BiocManager::install("Biobase")
# BiocManager::install("multtest")
# install.packages("GAPIT3")

# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(data.table)
library(vcfR)   
library(tidyverse)
library(dplyr)

pheno<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae.csv", sep = ",")        

pheno<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/emmeans_pisum.csv", sep = ",")        


pheno<- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                                      sep =",", 
                                      stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., pheno, by = "ivip")

pheno <- subset(pheno, cluster_name %in% c("North", "Winter"))
pheno <- subset(pheno, cluster_name %in% c("South", "South_America"))


pheno <- pheno[,c(1,3)]


nrow(pheno)

myGD1 <-fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column


myGD2 <-fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.indv",header = FALSE,stringsAsFactors =  FALSE)

myKI <- fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5_kinship.cXX.txt")
# Load datasets



myGM <-fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.pos") %>% rename(
    "chromosome" = "V1",
    "position" = "V2"
  )%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 

colnames(myGD2)[1]<-"ivip"
myY <- merge(myGD2, pheno, by = "ivip", all.x = TRUE)

colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)

# Add snp_ID column from dataset2 as the first column of dataset1
myKI <- cbind(snp_id = myY$ivip, myKI)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  KI = myKI,
  PCA.total=3,
  model=c("MLM", "FarmCPU", "Blink"),
  Multiple_analysis=TRUE)


myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  KI = myKI,
  PCA.total=3,
  model=c("Blink"),
  cutOff = 1e-5)

 as.data.frame(test) <- myGD$chr5_1349197397
test<-as.data.frame(myGD$chr5_1349197397)
test2 <-cbind(test,myY)
