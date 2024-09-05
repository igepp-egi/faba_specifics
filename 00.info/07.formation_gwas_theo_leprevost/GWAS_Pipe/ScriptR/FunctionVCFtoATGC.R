######################################################
## Convert VCF File for quality control using R
#####################################################
#####################################################
## Author: Theo Leprevost
####################################################
# Date: May 2023

### README !!!
## AIM: 
# Convert VCF file for quality control using R

## INPUT:
# (1) "data" = A VCF  matrix
# (2) "Names_Lines" = The names of the individuals

## OUTPUT:
# (1) Creation of a Rdata file including the genotyping matrix (missing values are coded NA), the genotyping code is in letters (ex. AG, CT, etc...), the markers are biallelic; the markers are in columns and the observations (individuals) are in rows.

# (2) Creation of a Rdata file including details on the markers i.e. 1st column "SNP.name", 2nd column "Chromosome" (chromosome number), 3rd column"Position" (physical position).

# Open a log file for writing
sink("log.txt", append = TRUE)

convert <- function(data, Names_Lines){
  
  #data loading
  library("stringr")
  library("dplyr")
  
  print(paste("The genotyping matrix contains:",length(data[,1]),"rows",length(data[1,]),"columns"))

  snp<-cbind(as.character(data[,3]),as.character(data[,1]),as.character(data[,2]))
  colnames(snp)<-c("SNP.name","Chromosome","Position")

  C<-length(data[1,])
  L<-length(data[,1])
  mat<-matrix(nrow=L,ncol=C)
  compt<-0

  for (i in (1:L)){
    o<-10
    for (j in (o:C)) {
      if ((str_sub(data[i,j], 1, 3)=="0/0")==TRUE) {
       mat[i,j]<-paste(as.character(data[i,4]),as.character(data[i,4]),sep="")
      }
      if ((str_sub(data[i,j], 1, 3)=="1/1")==TRUE) {
       mat[i,j]<-paste(as.character(data[i,5]),as.character(data[i,5]),sep="")
      }
      if ((str_sub(data[i,j], 1, 3)=="0/1")==TRUE) {
       mat[i,j]<-paste(as.character(data[i,4]),as.character(data[i,5]),sep="")
      }
      if ((str_sub(data[i,j], 1, 3)=="1/0")==TRUE) {
       mat[i,j]<-paste(as.character(data[i,5]),as.character(data[i,4]),sep="")
      }
      if ((str_sub(data[i,j], 1, 3)!="0/0")&(str_sub(data[i,j], 1, 3)!="1/1")&(str_sub(data[i,j], 1, 3)!="0/1")&(str_sub(data[i,j], 1, 3)!="1/0")==TRUE) {
        mat[i,j]<-NA
      }
    } 
    if ((i==(compt+(round(L/20,0))))==TRUE) {
      compt<-compt+(round(L/20,0))
      print(paste("The line number ",compt," of the matrix has been processed"))
    }
  }
  geno<-mat[,-(1:9)]
  geno<-as.data.frame(geno)
  geno<-t(geno)
  colnames(geno)<-snp[,1]
  rownames(geno)<-Names_Lines[,1]
  
  print("Nine informative columns were removed and the genotyping matrix was transposed")
  
  print(paste("The genotyping converted matrix contains :",length(geno[,1]),"rows",length(geno[1,]),"columns"))
  
  save(geno, file=paste("./AMS_v2_ATGC_",length(geno[,1]),"Lines_",length(geno[1,]),"SNPs.RData",sep=""))
  save(snp, file=paste("./AMS_v2_Info_",length(snp[,1]),"SNPs.RData",sep=""))
}

# Close the log file
sink()
