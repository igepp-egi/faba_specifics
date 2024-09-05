####################### Here we can load the last released version of R software ################

#### Haplotype contruction using Haplo_Construction function

##################
col<-"AMS"
# Open a README file for writing
sink(paste(col,"_Haplotype_Analysis_Results.txt",sep=""), append = FALSE)

#Genotyping dataset
load("../9_AMS_v2_filtered_variants_again/Results_QC/geno_QC_212Lines_249247SNPs.Rdata")
geno<-geno.pruned.QC

print("Dimension of the genotyping file")
geno[1:5,1:5]
dim(geno)  #212 249247

#Significant SNP detected in GWAS
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

#List of selected genotypes to color in blue the associated dots in the haploplot
print("Lines sharing the favorable allele at Chr7_279743119 to color in blue in the haploplot")
Blue_lines_AMS_DS_RB84 <- read.csv("./Haplotype_Analysis_blue_lines_AMS_DS_RB84.csv", sep=";")
print(Blue_lines_AMS_DS_RB84)
Blue_lines_AMS_DS_RB84 <- Blue_lines_AMS_DS_RB84$Code

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

#List of selected genotypes to color in blue the associated dots in the haploplot
print("Lines sharing the favorable allele at Chr4_80357270 to color in blue in the haploplot")
Blue_lines_AMS_DS_Ae109 <- read.csv("./Haplotype_Analysis_blue_lines_AMS_DS_Ae109.csv", sep=";")
print(Blue_lines_AMS_DS_Ae109)
Blue_lines_AMS_DS_Ae109 <- Blue_lines_AMS_DS_Ae109$Code

### Haplotype analysis

Variable<-c("AMS_DS_RB84",
            "AMS_DS_Ae109")

source("Function_Haplo_Construction.R")
source("Function_Haplo_Comparison.R")

#Import the BigLD resuls
Gresf1_15_LD_500 <- read.csv("./Results_BigLD/AMS_BigLD_LDblocks_15_LD_500.csv", sep=";")

#If no lines are selected, write:
#Blue_lines<-NULL

for (j in 1:length(Variable)) {
  for (k in 1:nrow(get(paste("GWAS_Results_",Variable[j],sep="")))) {
    print("Studied variable:")
    print(Variable[j])
    print("Haplotype analysis for the marker:")
    print(get(paste("GWAS_Results_",Variable[j],sep=""))[k,1])
    print("LDblock around the selected marker:")
    filtered_df <- Gresf1_15_LD_500[Gresf1_15_LD_500$chr == as.numeric(get(paste("GWAS_Results_",Variable[j],sep=""))[k,2]) & as.numeric(get(paste("GWAS_Results_",Variable[j],sep=""))[k,3]) >= Gresf1_15_LD_500$start.bp & as.numeric(get(paste("GWAS_Results_",Variable[j],sep=""))[k,3]) <= Gresf1_15_LD_500$end.bp,]
    print(filtered_df)
    geno_test<-geno[,(which(colnames(geno) == filtered_df$start.rsID)):(which(colnames(geno) == filtered_df$end.rsID))]
    file<-paste("./Results_Haplo_Construction/",Variable[j],"_",get(paste("GWAS_Results_",Variable[j],sep=""))[k,1],sep = "")
    
    #Haplotype construction
    Haplo_Construction(data=geno_test,file=file)
    
    #Haplotype comparison
    #Load the emmeans associated with each line
    Emmeans_Variable <- read.csv("../10_AMS_v2_LD_kinship_structure/LD/Emmeans_Disease_Severity_212_ordered_Renamed.csv", sep=";")
    Emmeans_Variable<- cbind(Emmeans_Variable$variety,eval(parse(text = (paste0("Emmeans_Variable$",Variable[j],sep="")))))
    colnames(Emmeans_Variable)<-c("Code",Variable[j])
    #Load the haplo detail file associated with each line
    Haplo_Detail <- read.csv(paste("./Results_Haplo_Construction/",Variable[j],"_",get(paste("GWAS_Results_",Variable[j],sep=""))[k,1],"_haplo_detail.csv",sep = ""), sep=";")
    #Match the emmeans computed for each line
    data <- merge(Emmeans_Variable, Haplo_Detail, by.x = "Code", by.y = "Code", all.x = TRUE)
    
    file<-paste("./Results_Haplo_Comparison/",Variable[j],"_",get(paste("GWAS_Results_",Variable[j],sep=""))[k,1],sep = "")
    
    #Select the lines to be colored in blue in the haploplot
    Blue_lines<-get(paste("Blue_lines_",Variable[j],sep=""))
      
    #Haplotype comparison
    Haplo_Comparison(data=data,file=file,colored_geno=Blue_lines,variable=Variable[j])
  }
}
# Close the README file
sink()




