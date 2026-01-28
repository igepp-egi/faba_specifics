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


######

pheno<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv", sep = ",")        

myGD1 <-fread(file = "../../results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column

myGD2 <-fread(file = "../../results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012.indv",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGM <-fread(file = "../../results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012.pos") %>% rename(
  "chromosome" = "V1",
  "position" = "V2"
)%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 

colnames(myGD2)[1]<-"ivip"
myY <- merge(myGD2, pheno, by = "ivip", all.x = TRUE)

colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)


myKI <- fread(file = "../../results_with_hedin/output/vfaba_GRM_fabae.cXX.txt")

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
  model=c("GLM","MLM", "MLMM", "FarmCPU", "Blink"),
  Multiple_analysis=TRUE)

# Load datasets

######
myGD1 <-fread(file = "../../files_GWAS/GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column


myGD2 <-fread(file = "../../files_GWAS/GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.indv",header = FALSE,stringsAsFactors =  FALSE)

myKI <- fread(file = "../../files_GWAS/GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5_kinship.cXX.txt")
# Load datasets



myGM <-fread(file = "../../files_GWAS/GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.pos") %>% rename(
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
  model=c("GLM","MLM", "MLMM", "FarmCPU", "Blink"),
  Multiple_analysis=TRUE)

myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  model=c("FarmCPU", "Blink"),
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


GAPIT.Power.compare(  
  GD=myGD,  
  GM=myGM,  
  nrep=100,  
  h2=0.586,  
  
  model=c("GLM",
          "MLM",
          "MLMM",
          "FarmCPU",
          "Blink"),
  NQTN=1
) 



FarmCPU.P.Threshold(
  Y=myY[,c(1,2)], #only two columns allowed, the first column is taxa name and the second is phenotype value
  GD=myGD,
  GM=myGM,
  trait="trait_name", #name of the trait, only used for the output file name
  theRep=30 #number of permutation times
)


myFarmCPU <- FarmCPU(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  p.threshold=0.5/nrow(myGM)
)
