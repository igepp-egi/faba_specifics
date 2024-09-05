####################### To utilize GPART, it is essential to operate within the R environment of version 4.0.2, as GPART has not been updated to be compatible with more recent releases of R ################

#loading the libraries

library(igraph)
library(plyr)
library(Homo.sapiens)
library(stringr)
library(gpart)

##################
col<-"AMS"
# Open a README file for writing
sink(paste(col,"_BigLD_LDblocks_Results.txt",sep=""), append = FALSE)

#SNP info dataset
SNPinfo <- read.csv("../10_AMS_v2_LD_kinship_structure/LD/Results_plink/SNP_Info_UniquePosition_Ordered_Plink.csv")
SNPinfo<-SNPinfo[,c(2,1,3)]
colnames(SNPinfo)<-c("chrN", "rsID", "bp")

print("Dimension of the SNP info file")
head(SNPinfo)
dim(SNPinfo) #249247 3

#Genotyping dataset
load("../10_AMS_v2_LD_kinship_structure/Kinship/Matrix_Genotyping_012NA_AMS_FW_.RData")
geno<-test2

print("Dimension of the genotyping file")
geno[1:5,1:5]
dim(geno)  #212 249247
identical(colnames(geno),as.character(SNPinfo$rsID))

#Gene dataset
Psat.2Cameor.genes <- read.delim("./CAMEOR_V2_ANNOTATION_1.gff3", header=FALSE)

# Keep only the genes from the gff3 file
Psat.2Cameor.genes <- Psat.2Cameor.genes[grepl("gene", Psat.2Cameor.genes$V3, ignore.case = TRUE), ]

# Modify the chromosome name, e.g. Lcu.2RBY.Chr1 in 1
for (k in (1:7)) {
  Psat.2Cameor.genes$V1<-gsub(paste("chr",k,sep=""), k, Psat.2Cameor.genes$V1)
  Psat.2Cameor.genes$V1<-as.factor(Psat.2Cameor.genes$V1)
}

# Remove markers located on scaffolds
Psat.2Cameor.genes <- Psat.2Cameor.genes[!grepl("scaffold", Psat.2Cameor.genes$V1, ignore.case = TRUE), ]

# Extract the values after "ID=" and before ";Name" for all rows
Psat.2Cameor.genes$V9 <- str_extract(Psat.2Cameor.genes$V9, "(?<=ID=)[^;]+")

# Remove unused levels from factor columns
Psat.2Cameor.genes <- droplevels(Psat.2Cameor.genes)

# Reset row names
rownames(Psat.2Cameor.genes) <- NULL

# Format the geneinfo table
geneinfo <-Psat.2Cameor.genes[,c(9,1,4,5)]
colnames(geneinfo)<-c("genename", "chrN", "start.bp", "end.bp")

print("Dimension of the gene info file")
dim(geneinfo) #48881     4
summary(geneinfo)
head(geneinfo)

print("To utilize GPART, it is essential to operate within the R environment of version 4.0.2, as GPART has not been updated to be compatible with more recent releases of R")

# Big LD algorithm using a 15 < LD bin size < 500 SNPS

Gresf1_15_LD_500 <- GPART(geno=geno, 
                           SNPinfo=SNPinfo, 
                           geneinfo=geneinfo, 
                           LD = "r2",
                           GPARTmode="LDblockBased",
                           Blockbasedmode="onlyBlocks",
                           minsize=15,
                           maxsize=500,
                           CLQmode = "maximal")

setwd("./Results_BigLD")
write.table(Gresf1_15_LD_500,file="AMS_BigLD_LDblocks_15_LD_500.csv",sep=";",row.names = FALSE)
# Close the README file
sink()

# Plot the LDheatmap
Gresf1_15_LD_500 <- read.csv("./AMS_BigLD_LDblocks_15_LD_500.csv", sep=";")

#AMS_DS_RB84
GWAS_SNP<-c("Chr3_433512360",
            "Chr4_364396365",
            "Chr5_24977625",
            "Chr6_12208861",
            "Chr7_279743119")


Chr_GWAS_SNP<-c("3",
                "4",
                "5",
                "6",
                "7")

Pos_GWAS_SNP<-c("433512360",
                "364396365",
                "24977625",
                "12208861",
                "279743119")



GWAS_Results_AMS_DS_RB84<-as.data.frame(cbind(GWAS_SNP,Chr_GWAS_SNP,Pos_GWAS_SNP))

#AMS_DS_Ae109
GWAS_SNP<-c("Chr1_426493035",
            "Chr2_488687021",
            "Chr3_443681564",
            "Chr4_199560606",
            "Chr4_397074433",
            "Chr4_80357270",
            "Chr4_80361300",
            "Chr5_487703199",
            "Chr5_67453776",
            "Chr6_379370857") 


Chr_GWAS_SNP<-c("1",
                "2",
                "3",
                "4",
                "4",
                "4",
                "4",
                "5",
                "5",
                "6")

Pos_GWAS_SNP<-c("426493035",
                "488687021",
                "443681564",
                "199560606",
                "397074433",
                "80357270",
                "80361300",
                "487703199",
                "67453776",
                "379370857") 

GWAS_Results_AMS_DS_Ae109<-as.data.frame(cbind(GWAS_SNP,Chr_GWAS_SNP,Pos_GWAS_SNP))

Variable<-c("AMS_DS_RB84",
            "AMS_DS_Ae109")

# Open a README file for writing
sink(paste(col,"_BigLD_LDHeatmap_Results.txt",sep=""), append = FALSE)

for (i in 1:length(Variable)) {
  for (k in 1:nrow(get(paste("GWAS_Results_",Variable[i],sep="")))) {
  print("Variable:")
  print(Variable[i])
  print("LD heatmap for the marker:")
  print(get(paste("GWAS_Results_",Variable[i],sep=""))[k,1])
  print("LDblock around the selected marker:")
  filtered_df <- Gresf1_15_LD_500[Gresf1_15_LD_500$chr == as.numeric(get(paste("GWAS_Results_",Variable[i],sep=""))[k,2]) & as.numeric(get(paste("GWAS_Results_",Variable[i],sep=""))[k,3]) >= Gresf1_15_LD_500$start.bp & as.numeric(get(paste("GWAS_Results_",Variable[i],sep=""))[k,3]) <= Gresf1_15_LD_500$end.bp,]
  print(filtered_df)
  geno_test<-geno[,(which(colnames(geno) == filtered_df$start.rsID)):(which(colnames(geno) == filtered_df$end.rsID))]
  SNPinfo_test<-SNPinfo[(which(SNPinfo$rsID == filtered_df$start.rsID)):(which(SNPinfo$rsID == filtered_df$end.rsID)),]
  LDblockHeatmap(geno = geno_test, 
                 SNPinfo = SNPinfo_test,
                 geneinfo = geneinfo,
                 blockresult = Gresf1_15_LD_500,
                 LD="r2",
                 geneshow = TRUE,
                 filename = paste(Variable[i],"_LDHeatmap_",get(paste("GWAS_Results_",Variable[i],sep=""))[k,1],sep=""),
                 type="png")
  }
}
# Close the README file
sink()





