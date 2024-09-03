#######################################################
##Author: Sandra Negro and Stephane Nicolas
##Adapted on virtual machine by Julie Ferreira de Carvalho
##ATTENTION : IMPORTANT AVANT UTILISATION
##PENSER A VERIFIER LES CHEMINS D'ACCES AUX DIFFERENTS FICHIERS/PROGRAMMES
##ET CHANGER LE NOM DE L UTILISATEUR
########################################################



########################################################
########################################################
## Senario Pipeline GWAS (mixed model)
########################################################
#########################################################
## DATE: 18 December 2018
## ATTENTION: You can find a detailed description of all inputs in the functions files.
## TESTED ON LINUX / NEED SOME ADAPTATION FOR WINDOWS
# Remove all object in R working space
rm(list=ls())



## Install and load the package required in pipeline
#install.packages("corpcor")
#install.packages("./PROGRAM/LIBRARY_R/GenABEL.data_1.0.0.tar.gz", repos = NULL, type = "source",dependencies = FALSE)
#install.packages("./PROGRAM/LIBRARY_R/GenABEL_1.8-0.tar.gz", repos = NULL, type = "source",dependencies = FALSE)

# Set working directory
#setwd("C:/Test_pipeline")
#create a folder for results:
Rep="~/GWASPipe/Result"
dir.create(Rep)
setwd(Rep)
getwd()

######################################################
## Descriptive analyses and quality control
#####################################################
### Path to the input files for the log record output file 
input.geno = "/home/gwas_user/GWASPipe/DATA/Matrix_geno_Lettres_randomized_247_5000.RData"
input.snp = "/home/gwas_user/GWASPipe/DATA/SNP_info_randomized.csv"

  ### Load the genotyping file
load(input.geno) # the R-object is "geno"
dim(geno)#247 x 5000
geno[1:5,1:5]
length(which(is.na(geno)))#49933 missing data
#colnames(geno)
#rownames(geno)
#List of genotype in genotyping matrix
unique(unlist(apply(geno,1,unique)))

  ### Load the marker file that includes detailed info
snp <- read.table(file="/home/gwas_user/GWASPipe/DATA/SNP_info_randomized.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
dim(snp)
head(snp)


	### Type of analysis to be written in the log file
Type.analysis = "Quality Control"
	### Parameters (see the description of the parameters in "Function_Quality_Control.r")
snp.dm = 0.20
maf = 0.01

snp.het = 0.15
ind.dm = 0.10
ind.het = 0.08

	### Load and call the function
source("/home/gwas_user/GWASPipe/ScriptR/Function_Quality_Control.r")
QC(geno, snp, input.geno, input.snp, Type.analysis, snp.dm, maf, snp.het, ind.dm, ind.het)

#####################################################
## BEAGLE v3; inferring missing data
#####################################################

# ATTENTION: Load BEAGLE (https://faculty.washington.edu/browning/beagle/b3.html#download) AND "beagle.jar" HAS TO BE IN THE WORKING DIRECTORY FOR BEAGLE TO WORK !!!

# Remove all objects in the working space
rm(list=ls())


	### Load the QCed matrix outputted from the "QC()" function above
load("./Results_QC/Matrix_QCed.RData") # the R-object is "geno.pruned.QC"
geno.pruned.QC[1:5,1:5]
dim(geno.pruned.QC)#247 x 4497
### number of Megabytes of memory available on your machine for the analyses
Mb = 15000 # check on your machine
# Rq: sur windows "Mb" semble etre limite a 1000, alors que sous linux on peut mettre la capacite reel de notre machine ! Demander au moulon a Olivier Langela ou Benoit si il peuvent aider?

	### Load and call the function
source("/home/gwas_user/GWASPipe/ScriptR/Function_Imputation_Beagle.r")
beagle(geno=geno.pruned.QC, Mb)

#############################################################
## Matrix conversion to plink binary files
#############################################################

### ATTENTION: the plink software should be install on your machine and for Windows, the executable file (plink.exe) should be in the working directory !!!

### Below there is a control to check if you have the same individual ID in the genotyping, kinships and phenotyping files.
### If there are differences between those files, the scenario below is extracting the same individuals for the following analyses (i.e. GWAS)

rm(list=ls())

	### Load the marker file that includes detailed info 
load(file = "./Results_QC/SNP_Info_UniquePosition_Ordered.RData")
dim(snp)#5000 x 3
head(snp)
tail(snp)
unique(snp$Chromosome)

	### Load the imputed genotyping file 
load("./Results_Beagle/Matrix_DMimputedBeagle.RData") # the R-obj is "geno"
geno[1:5,1:5]
length(which(is.na(geno)))#0 missing data
ls()

	### Nb of chromosome (= Nb of chr present in the SNP info; in maize chr=10 or chr=11 if you have unplaced SNPs put in a fictive chr11)
chr=10

	### Load the phenotyping file to extract the same ind. than in the genotyping file, this is required for the GWAS
# Load the phenotype
pheno = read.table(file="/home/gwas_user/GWASPipe/DATA/Phenotypes_247.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
dim(pheno)#247 x 17

	### Check if both genotype and phenotype matrices have the same individual ID and in the same order
com.geno.pheno = intersect(rownames(geno), pheno$variety)

if(length(com.geno.pheno) == 0){
	cat("\n")
	cat("\n")
	print("#######################################################################")
	print("The individual ID are totally different between the phenotyping and genotyping")
	print("#######################################################################")
	cat("\n")
	cat("\n")
} else {
  	### Check if there is the same nb of individuals and check if the individual ID are in the same order
 	geno=geno[rownames(geno) %in% pheno$variety, ]
	
  	pheno=pheno[match(rownames(geno), pheno$variety), ]
  
  	print("#####################################################")
  	print(paste("There are", nrow(pheno), "individual in common"))
  	print("#####################################################")
}

	### Save the phenotyping file that has the individual names in the same order than geno
write.table(pheno, file=paste0("/home/gwas_user/GWASPipe/DATA/Phenotypes_", nrow(pheno),"_ordered.csv"), sep=",", col.names=TRUE, quote=FALSE)

	### Load and call the function to create the ouput files
source("/home/gwas_user/GWASPipe/ScriptR/Function_Transfo_plink_conda.r")
Conversion(snp, geno, chr)

#############################################################
# Admixture: population structure analyses
#############################################################
## ATTENTION: ADMIXTURE DO NOT WORK ON WINDOWS, it has been written to work on Linux or MacOS
rm(list=ls())
#Rep="/home/formation/Bureau/Pipeline_WithoutResult/"
#setwd(Rep)


	### Load the genotyping matrix
#load("./Results_plink/Matrix_247_4497.RData") # the R-object is "geno" here
load("/home/gwas_user/GWASPipe/Result/Results_plink/MatrixGeno.RData") # the R-object is "geno" here
geno[1:5,1:5]

	### Path to the Admixture software. If you do not know, open a new terminal and type the command line "which admixture". You will have the path
#path = "/home/formation/Bureau/Pipeline_WithoutResult/PROGRAM/admixture" # Linux
#path = "../PROGRAM/admixture" # Linux

path = "/usr/local/env/admixture_1.3/bin/admixture"

### Serie of sub-groups to testls

Q = "1 2 3 4 5 6"

	### Last number of the serie "Q"
Qmax = 6

	### number of sub-groups we think there is for the stacked bar plot
nbQ=6
	### Ancestral fraction threshold for the stacked bar plot
Q.threshold=0.3

### Nb of dimension to include in covariate in your GWAS
nbk = 3 

	### Call and run the functions
source("/home/gwas_user/GWASPipe/ScriptR/Function_Admixture.r")

# On the screen, you will also have the values of the CV (from 1 to 10 sub-groups) that you can put on a graph: the K on the x-axis and the CV errors on the y-axis
admixture(path, Q, Qmax, geno)

# Stacked bar plot
StackedBar(geno, nbQ, Q.threshold)

##########################################################
## Kinship estimation (IBS and IBD) using GenABEL
##########################################################
rm(list=ls())

## Install and load the package
#install.packages("GenABEL")
#install.packages("/home/gwas_user/GWASPipe/PROGRAM/LIBRARY_R/GenABEL.data_1.0.0.tar.gz", repos = NULL, type = "source",dependencies = FALSE)
#install.packages("/home/gwas_user/GWASPipe/PROGRAM/LIBRARY_R/GenABEL_1.8-0.tar.gz", repos = NULL, type = "source",dependencies = FALSE)
## Beware, if you cannot install packages, use tools: install packages: install from package archive file and select both GenABEL files in folder: LIBRARY_R

library("GenABEL")

# Load the .ped .map matrix
pedfile = "./Results_plink/MatrixGeno.ped"
mapfile = "./Results_plink/MatrixGeno.map"
ped = read.table(file=pedfile, sep="", header=FALSE, stringsAsFactors=FALSE)

dir.create("/home/gwas_user/GWASPipe/Result/Results_Kinship_GenABEL")

# Load the function
source("/home/gwas_user/GWASPipe/ScriptR/Function_Kinship_GenABEL.r")

# Call the function "kin.genable()"
kin.genable(ped, pedfile, mapfile)

####################################################################################
# Estimation of kinships without the chromosome tested (Rincent et al. 2014 approach)
####################################################################################
rm(list=ls())
#Rep="/home/formation/Bureau/Pipeline_WithoutResult/"
#setwd(Rep)

# This function require a numeric genotyping matrix (Dose code: 0,1,2)
# Load Stephane function to convert a lettre matrix to a numeric matrix
source("/home/gwas_user/GWASPipe/ScriptR/Function_MatrixLetter2Numeric012NA.r")

# Load the genotyping file generated from the plink function
load("./Results_plink/MatrixGeno.RData")
ls()
geno[1:5,1:5]
# Call the function to convert the lettre to dose code
mat012=MatrixLetter2Numeric012NA(test=geno,Panel="TD",Origin="FW") 


# Load the dose code (0, 1, 2) converted genotype matrix
load(file="Matrix_Genotyping_012NA_TD_FW_.RData") # the R-obj is "test2"
test2[1:5,1:3]


# Load the SNP info file (Chr nb and position)
snp <- read.table(file="./Results_plink/SNP_Info_UniquePosition_Ordered.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
head(snp)

# Nb of chromosome (Nb of Chr present in the SNP info)
chr = 11
# choose the algorithm to be computed
algo = "Astle" # Astle and Balding 2009
#algo = "VR" # Van Raden et al. 2008

# Call and run the function
rm(mat012)
source("/home/gwas_user/GWASPipe/ScriptR/Function_Kinship_Renaud_ssChr.r")
kin(geno = test2, algo, snp, chr)

##################################################################
## Principal Coordinates Analyses (PCoA): Analysing similarity between inbred lines
##################################################################
rm(list=ls())
#Rep="/home/formation/Bureau/Pipeline_WithoutResult/"
#setwd(Rep)

# Load the package "stats"
library("stats")

	### Load the general kinship matrix
kin = read.table(file="./Results_Kinship_ssChr/IBD_General_Renaud_Astle_247.csv", sep=",", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
kin[1:5,1:5]
	### Store the Admixture output name with the path
# Without the extension
adm.name = "./Results_Admixture/MatrixGeno"

	### Nb of genetic groups you think you have
nbQ = 6 
	
### Nb of dimension to include in covariate in your GWAS
nbdim = 3 

	### Specify if the user has his/her own classification file
	user.classification = FALSE #(or TRUE)

	### If the user.classification = TRUE then load the file in "user.class" that includes the individual names in the first column and the classification in the second column 
	if(user.classification){
	user.class = read.table("./DATA/Info_Lignees.txt", sep="\t", header=TRUE, stringsAsFactors=TRUE)
	}

	### Call and run the function
source("/home/gwas_user/GWASPipe/ScriptR/Function_PCoA.r")
PCoA(kin, adm.name, nbQ, user.classification)



#########################################################
## Genome-wise association studies with FASTLMM using Rincent et al. 2014 approach
# (kinship estimated by excluding the chromosome of the tested marker) 
#########################################################
## ATTENTION: FASTLMM SHOULD BE INSTALL
## If your operating system is Windows, the executable file (FastLmmC.exe) should be in the working directory. 
rm(list=ls())
Rep2="/home/gwas_user/GWASPipe/"
setwd(Rep2)


	### Install and load the package corpcor
#install.packages("corpcor")
library("corpcor")

	### Load the input files
# Load the phenotype
# file name of pheno to be written in the log information
pheno.file = "/home/gwas_user/GWASPipe/DATA/Phenotypes_247_ordered.csv"
pheno = read.table(file=pheno.file, sep=",", header=TRUE, stringsAsFactors=FALSE)
dim(pheno)


# Load the kinships (for ex. obtained from the script: "Kinship_K_Chr_Renaud.r")
# File name with the entire path of the kinship to be written in the log information
kin.file = "/home/gwas_user/GWASPipe/Result/Results_Kinship_ssChr/IBD_Renaud_Astle_247.RData"
load(kin.file)
kin = matA1_sansChrom
kin[[1]][1:5,1:5]

# Load the structure (for ex. obtained from the script: "PCoA.r or Admixture")
# File name with the entire path of the structure to be written in the log information
struct.file = "/home/gwas_user/GWASPipe/Result/Results_PCoA/covarforFaST_LMM.txt" # PCoA Matrix
#or 
#struct.file = "/home/gwas_user/GWASPipe/Result/Results_Admixture/covarforFaST_LMM.txt" # Admixture Matrix
struct = read.table(file=struct.file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
head(struct)

# Binary file name without the chromosome number and the extension. The binary files can be obtained using the function: "Conversion()"
geno.file = "Result/Results_plink/MatrixGeno_Chr"

	### Type of analysis to be written in the log file
type.analysis = "Linear mixed model in FASTLMM for GWAS using a test data set and Rincent et al. 2014 approach"

	### Nb of chromosome (=Nb of chr present in the SNP info file file)
chr = 10
	### GWAS with fastlmm

#path ou se trouve fastlmm (linux)
path="/home/gwas_user/GWASPipe/PROGRAM/fastlmmc"
#system("chmod +x /home/gwas_user/GWASPipe/PROGRAM/fastlmmc")


#Load the gwas.fastlmm function
source("/home/gwas_user/GWASPipe/ScriptR/Function_FASTLMM_GWAS.r")
# If you don't want to use a structure matrix as covariate, type struct=NULL
#struct=NULL
gwas.fastlmm(pheno, pheno.file, kin, kin.file, struct, struct.file, geno.file, chr, type.analysis, path)

#################################################################
## GNPasso template: GWAS output template from FASTLMM results ##
#################################################################
rm(list=ls()) # remove all object in R
library("stats")
Rep3="~/GWASPipe/Result"
setwd(Rep3)
getwd()

### Load the inputs
# Load the phenotyping file
pheno = read.table(file="../DATA/Phenotypes_247_ordered.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
head(pheno)

# Load the genotyping file that contains all the markers tested
load("Results_plink/MatrixGeno.RData") # the R-object is "geno"
geno[1:5,1:3]

# Bim file containing the SNPs from all chromosomes
bim.all = read.table(file="Results_plink/MatrixGeno.bim", sep="\t", header=FALSE, stringsAsFactors=FALSE)

# User parameters
Experimenter = "Toto"
Exper_Name = "GWAS_Essai"
Panel_Name = "Panel"
Marker_Type = "SNP"
Software = "FASTLMM"
Statistical_Model = "K"
RefGenMap = "Ref_v2"
# number of chromosome
Nchr = 10

#### Call the function "gnpasso" to create the template table
# Load the gnpasso function
source("../ScriptR/Function_GNPasso_Template.r")
gnpasso(pheno, geno, bim.all, Experimenter, Exper_Name, Panel_Name, Marker_Type, Software, Statistical_Model, RefGenMap, Nchr)

#Output files, Rdata, table with all SNP info in 'template'

#################################################################
## Vizualisation from E. Corlouer ##
#################################################################
rm(list=ls())
#Exemple avec un fichier GNPasso, une variable phénotypique

load("Results_GNPasso/GNPasso_biogai.2012.drought.dta.RData")
#install.packages("qqman")

data <- as.data.frame(template)
rm(template)

for (k in (4:12)) {
  data[,k]<-as.numeric(as.character(data[,k]))
}
for (k in (15:20)) {
  data[,k]<-as.numeric(as.character(data[,k]))
}

data_Adjusted <- data
write.csv(data_Adjusted, file = "biogai.2012.drought.dta.csv", quote = F)


library(qqman)
qq(data_Adjusted$Pvalue, main="2012.drought.dta")
manhattan(data_Adjusted, chr = "Chromosome", bp="Marker_Position", snp = "Marker_Name", p="Pvalue", main.cex = 0.5, main = "MT_TA TSW GWAS FastLMM with LOO Chromosome procedure",col = c("#0F4C81", "#87A5C0","#0F4C81", "#87A5C0","#0F4C81", "#87A5C0","#0F4C81", "#87A5C0","#0F4C81", "#87A5C0","#C00000", "#E59999","#C00000", "#E59999","#C00000", "#E59999","#C00000", "#E59999","#C00000"), suggestiveline = 5.55344086, genomewideline = 6.561614121)
rm(list=ls())
