###########################################################
###########################################################
# Estimation of the kinships without the chromosome tested #
###########################################################
# Author: Renaud Rincent modified as a function and optimized in term of speed calculation by Sandra Negro
###########################################################
## Aim: 
# Estimation of IBD matrices
# Choice between the algorythms from Astle and Balding 2009 or Van Raden et al. 2008

## INPUT:
# (1) "geno" = genotyping file with code 0, 1, 2 (0 and 2 are homozygotes for allele 1 and allele 2 respectively, 1 is the heterozygote). No missing data are allowed. Individuals in lines and markers in columns.
# Rq: for the 50K Illumina chip, PANZEA markers are used for the estimation of the kinship!
#
# (2) "algo" = algorithm to choose: "Astle" or "VR" (Astle and Balding 2009 or Van Raden et al. 2008)
#
# (3) "snp" = marker information file including the column's names: "SNP.name", "Chromosome", "Position". Unplaced markers should be in a fictive chromosome (i.e. chr11 for maize)
#
# (4) "chr" = nb of chromosomes present in the SNP info file (in maize chr=10 or chr=11 if there are unplaced SNP put in a fictive chromosome 11) 

## OUTPUT:
# (1) A general kinship matrix estimated from all chromosomes ("IBD_General_Renaud_Astle_(nb.Individuals).csv")
#
# (2) A list containing one matrix by chromosome tested ("IBD_Renaud_Astle_(nb.Individuals).RData"), if you have unplaced SNP in a fictive chromosome (ex.: chr11 for maize) then you will have 11 kinship matrices. The first matrix of the list was estimated from all chromo except the chr1
#
# (3) A log record file containing information on the analyses i.e. date, mean and range of the general kinship
#
# (4) An histogram of the relatedness coefficient of the general kinship


## REFERENCES:
# Rincent et al. (2014). Recovering Power in Association Mapping Panels with Variable Levels of Linkage Disequilibrium. Genetics, vol. 197 no. 1 375-387.
#
# Astle and Balding (2009). Population Structure and Cryptic Relatedness in Genetic Association Studies. Statistical Sciences, Vol. 24, No. 4, 451–471.
#
# VanRaden (2008). Efficient Methods to Compute Genomic Predictions. Journal of Dairy Science Vol. 91 No. 11 4414–4423.
############################################################
############################################################

kin <- function(geno, algo, snp, chr){
	### Create folder "Kinship"
	dir.create("Results_Kinship_ssChr")
	### Convert the snp in data.frame and geno in matrix
	snp=as.data.frame(snp)
	geno=as.matrix(geno)
	row=rownames(geno)
	geno=apply(geno,2, as.integer)
	rownames(geno)=row
	### Check that the markers are in order (chr number and position)
	snp= snp[order(snp$Chromosome, snp$Position), ]
	### Extract the same SNPs between snp and geno
	snp = snp[snp$SNP.name %in% colnames(geno), ]
	geno = geno[, match(snp$SNP.name, colnames(geno))]	

	### Tranpose the genotype and add the a Chromosome and Position columns
	dataPZE = t(geno)
	dataPZE = cbind(snp$Chromosome, snp$Position, dataPZE)
	### Converting geno from dose code to frequency code (0, 0.5, 1)
	ZdosePanzea = geno/2 

	#Freq=apply(ZdosePanzea,2,function(a) mean(a,na.rm=TRUE))
	Freq=apply(ZdosePanzea,2,mean)

	Denom=Freq*(1-Freq)
	ZdosePanzea3=ZdosePanzea
	ZdosePanzea2=ZdosePanzea
	for (t in 1:ncol(ZdosePanzea3)) {
		ZdosePanzea2[,t]=ZdosePanzea[,t]-Freq[t]
		ZdosePanzea3[,t]=(ZdosePanzea3[,t]-Freq[t])/(Denom[t])
	}


	if(algo=="Astle"){	
		### IBD global matrix (Astle and Balding 2009)
		matA1=tcrossprod(ZdosePanzea2, ZdosePanzea3)/(ncol(ZdosePanzea2))

		### Methode IBD (Astle and Balding 2009) sans Chrom teste
		matA1_sansChrom=list()
		for (ii in 1:length(unique(dataPZE[,1]))) {
			test=which(dataPZE[,1]==ii)
			matA1_sup=tcrossprod(ZdosePanzea2[,test], ZdosePanzea3[,test])
			matA1_sansChrom[[ii]]=(ncol(ZdosePanzea3)*(matA1)-matA1_sup)/(ncol(ZdosePanzea3)-length(test))  
		}
	}else{
		### IBD matrix (Van Raden 2008)
		matA1=tcrossprod(ZdosePanzea2, ZdosePanzea2)/sum(Denom)

		### IBD matrix (Van Raden 2008) sans Chrom teste
		matA1_sansChrom=list()
		for (ii in 1:length(unique(dataPZE[,1]))) {
			test=which(dataPZE[,1]==ii)
			
			matA1_sup=tcrossprod(ZdosePanzea2[,test],ZdosePanzea2[,test])
			matA1_sansChrom[[ii]]=(sum(Denom)*(matA1)-matA1_sup)/(sum(Denom[-test]))
		}

	}

	### Save the general Kinship
	write.table(matA1, file=paste0("./Results_Kinship_ssChr/IBD_General_Renaud_", algo, "_", nrow(geno),".csv"), sep=",", dec=".", col.names=TRUE, row.names=TRUE, quote=FALSE)

	### Save the kinships estimated with all the chr except the chr of the tested snp
	save(matA1_sansChrom, file=paste0("./Results_Kinship_ssChr/IBD_Renaud_", algo, "_", nrow(geno),".RData"))
	
	
	
	### Range and mean
	range.diag.ibd = range(diag(matA1))
	range.ibd = range(matA1[lower.tri(matA1, diag=FALSE)])
	mean.ibd = mean(matA1[lower.tri(matA1, diag=FALSE)])
	

	### Save the log file
	log=rbind(as.character(Sys.Date()), "Analyses: Estimating IBD relatedness matrices", "", "Range of the diagonal of the general IBD: ", range.diag.ibd[1], range.diag.ibd[2], "Range of the relatedness coefficient (IBD): ", range.ibd[1], range.ibd[2], paste0("Mean IBD: ", mean.ibd))

	write.table(log, file="./Results_Kinship_ssChr/Log.txt", col.names=FALSE, row.names=TRUE, quote=FALSE)


	### Histogram
	pdf("./Results_Kinship_ssChr/Hist_Distr_General_kinship.pdf")
	hist(matA1[lower.tri(matA1, diag=FALSE)], main = "Distribution of the IBD", xlab="IBD")
	dev.off()

}



