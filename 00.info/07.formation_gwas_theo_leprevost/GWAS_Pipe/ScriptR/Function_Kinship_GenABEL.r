########################################################################
#########################################################################
## Estimating IBS and IBD matrices using GenABEL ##
#########################################################################
# Author: Sandra Negro
########################################################################
## DATE: 08 October 2015

## AIM:
# Estimation of kinships (IBS or IBD) using GenABEL with plink file format
# and Rincent et al. (2014) approach

## ATTENTION: Package "GenABEL" needed

## INPUT: 
# (1) ped = a genotyping .ped (plink format). You can obtain the plink files using the function "Conversion()" from the file "Function_Transfo_plink.r". The first six columns are: Family ID, Individual ID, first parent ID (use 0, which is considered as missing), Maternal ID (use 0), Sex (use 0), Phenotype. From the seventh and on there are the genotypes. You can have missing genotypes(coded 0 0). 
# Rq: PANZEA markers are used for the estimation of the kinship.

# (2) pedfile = name of the .ped file with the path to the file

# (3) mapfile = name of the .map file with the path to the file

## OUTPUT: 
# (1) A IBD and IBS general kinship matrices
# (2) a log file with details on the analyses and the properties of the kinships
# (2) graph distribution kinship
########################################################################

kin.genable <- function(ped, pedfile, mapfile){

	library(GenABEL)
	dir.create("Results_Kinship_GenABEL")


	### Build a phenophile. You need a phenofile for the GenABEL function "load.gwaa.data" below. There should be column named "id" and "sex".
	id=as.character(ped[,2])
	sex=as.vector(rep(1,length(id)))
	Pheno=as.vector(rep(1,length(id)))
	phenofile=as.data.frame(cbind(id,sex,Pheno))
	head(phenofile)[1:5, 1:3]
	write.table(phenofile, file="./Results_Kinship_GenABEL/phenofile.txt", sep="\t", 	row.names=FALSE, col.names=TRUE, quote=FALSE)


	### Convert the data from plink format into GenABEL format and create an image
	convert.snp.ped(pedfile=pedfile, mapfile=mapfile, outfile="./Results_Kinship_GenABEL/genofile.raw")

	### Combines the phenotypic and genotypic data into an R object
	geno.pheno=load.gwaa.data(phe="./Results_Kinship_GenABEL/phenofile.txt", gen="./Results_Kinship_GenABEL/genofile.raw")

	### This creates a relatedness matrix based on pairwise IBS sharing
	gkin.ibs=ibs(geno.pheno, w="no")
	### This creates a relatedness matrix frequency-weighted (IBD)
	gkin.ibd=ibs(geno.pheno, w="freq")
#Description of the matrix output: below the diagonal you have the average IBS or IBD values (depending on the type of matrix), and above the diagonal you have the number of SNP genotype measured for both members of the pair.


	###################################################################
	#Conversion from GenABEL matrix of similarity to a matrix of 2 identical triangles with a diagonal ~=1 (if not equal to 1, it is because of the heterozygotes)
	#################################################################### 
	ibd=as.matrix(gkin.ibd)
	ibs=as.matrix(gkin.ibs)

	# Create a matrix of 0 
	sous.ibd=matrix(0,nrow(ibd), ncol(ibd))
	sous.ibs=matrix(0,nrow(ibs), ncol(ibs))
	sous.ibd.2=matrix(0,nrow(ibd), ncol(ibd))
	sous.ibs.2=matrix(0,nrow(ibs), ncol(ibs))


	# Fill the new matrix with the lower triangle of the kinship matrix without the diagonal
	sous.ibd[lower.tri(sous.ibd,diag=FALSE)] = ibd[lower.tri(ibd,diag=FALSE)]
	sous.ibs[lower.tri(sous.ibs,diag=FALSE)] = ibs[lower.tri(ibs,diag=FALSE)]
	# Fill the new matrix with the lower triangle of the kinship matrix with the diagonal
	sous.ibd.2[lower.tri(sous.ibd.2,diag=TRUE)] = ibd[lower.tri(ibd,diag=TRUE)]
	sous.ibs.2[lower.tri(sous.ibs.2,diag=TRUE)] = ibs[lower.tri(ibs,diag=TRUE)]

	### Building a symetrical kinship matrix
	new.ibd=sous.ibd + t(sous.ibd.2)
	colnames(new.ibd)=colnames(ibd)
	rownames(new.ibd)=rownames(ibd)

	new.ibs=sous.ibs + t(sous.ibs.2)
	colnames(new.ibs)=colnames(ibs)
	rownames(new.ibs)=rownames(ibs)


	### Range and mean
	range.diag.ibd = range(diag(new.ibd))
	range.ibd = range(new.ibd[lower.tri(new.ibd, diag=FALSE)])
	mean.ibd = mean(new.ibd[lower.tri(new.ibd, diag=FALSE)])

	range.diag.ibs = range(diag(new.ibs))
	range.ibs = range(new.ibs[lower.tri(new.ibs, diag=FALSE)])
	mean.ibs = mean(new.ibs[lower.tri(new.ibs, diag=FALSE)])
	

	### Save the Kinship
	# IBD
	write.table(new.ibd, file="./Results_Kinship_GenABEL/IBD_genABEL.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=FALSE)
	# IBS
	write.table(new.ibs, file="./Results_Kinship_GenABEL/IBS_genABEL.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=FALSE)
	
	### Save the log file
	log=rbind(as.character(Sys.Date()), "Estimating IBS and IBD relatedness matrices", "", "INPUT FILE:", pedfile, mapfile, "", "Range of the diagonal of the IBD: ", range.diag.ibd[1], range.diag.ibd[2], "Range of the diagonal of the IBS: ", range.diag.ibs[1], range.diag.ibs[2], "", "Range of the relatedness coefficient (IBD): ", range.ibd[1], range.ibd[2], "Range of the relatedness coefficient (IBS): ", range.ibs[1], range.ibs[2], paste0("Mean IBD: ", mean.ibd), paste0("Mean IBS: ", mean.ibs))

	write.table(log, file="./Results_Kinship_GenABEL/Log.txt", col.names=FALSE, row.names=TRUE, quote=FALSE)


	### Histogram
	pdf("./Results_Kinship_GenABEL/Hist_Distr_kinships.pdf")
	par(mfrow=c(1,2))
	hist(new.ibd[lower.tri(new.ibd, diag=FALSE)], main = "Distribution of the IBD", xlab="IBD")
	hist(new.ibs[lower.tri(new.ibs, diag=FALSE)], main = "Distribution of the IBS", xlab="IBS")

	dev.off()


}

