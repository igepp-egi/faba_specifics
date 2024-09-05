#####################################################
## BEAGLE v3; inferring missing data
#####################################################
#####################################################
## Author: Sandra Negro
#####################################################
# Date: 30 September 2015

## AIM: 
# Inference of missing data with BEAGLE version 3 !

## BEAGLE is written in Java and runs on most computing platforms (e.g. Windows, Unix, Linux, Solaris, and Mac). If Java is not installed or if it is not version 1.6 or later, the current Java interpreter can be downloaded from the java.sun.com web site. The Java interpreter is called the Java Standard Edition (SE) Runtime Environment (JRE).
## BEAGLE can be downloaded from BEAGLE web site (https://faculty.washington.edu/browning/beagle/b3.html#download)

# Rq: "beagle.jar" HAS TO BE IN THE WORKING DIRECTORY FOR BEAGLE TO WORK or specify the path !!!


## INPUT: 
# (1) "geno" = genotyping file; individuals in rows and markers in columns, genotyping code can be any sequence of characters, it is not restricted to A/C/G/T or to 1/2/3/4 (but each marker have to be biallelic). Missing data (NA) are allowed
#
# (2) "Mb" = number of Megabytes of memory available (ex. Mb=7000)


## OUTPUT: 
# (1) Creation of the "Results_Beagle" folder
# (2) "Matrix_BEAGLE.bgl" = a Beagle genotyping file. First column describes the data on each line: "I" denotes sample identifier data and an "M" denotes marker data. The second column contains the SNP names. From the third column and onwards there are the genotypes: one column per allele
#
# (3) "Matrix_DMimputedBeagle.RData" = a conversion from Beagle output to a genotyping matrix (both alleles of the SNP in one column) with individuals in rows and markers in columns
#
# (4) "Output.log" = a log record file
#
# (5) Four other outputs: "Output.Matrix_BEAGLE.bgl.dose.gz", "Output.Matrix_BEAGLE.bgl.gprobs.gz", "Output.Matrix_BEAGLE.bgl.phased.gz", "Output.Matrix_BEAGLE.bgl.r2". See BEAGLE 3.3.2 manual for a complete description of those files.

## REFERENCES
# S R Browning and B L Browning (2007) Rapid and accurate haplotype phasing and missing data inference for whole genome association studies using localized haplotype clustering. Am J Hum Genet 81:1084-1097.doi:10.1086/521987

#####################################################

beagle <- function(geno, Mb){
	geno=as.matrix(geno)

	### Create the result folder
	dir.create("Results_Beagle")
	

	### Remplace the NA in ?? for beagle format
	rempl=function(mat) {
		mat[which(is.na(mat))]="??"
		return(mat)  
	}
	geno=apply(geno,2,rempl)


	### Create two empty matrices
obj2=matrix(NA,nrow=nrow(geno),ncol=ncol(geno))
obj3=matrix(NA,nrow=nrow(geno),ncol=ncol(geno))

	# Store the first allele in one matrix and the second allele in the second matrix
	for(i in 1:ncol(geno)) {
		obj2[,i]=substr(geno[,i],1,1)
		obj3[,i]=substr(geno[,i],2,2)
	}

	# Add the rownames on the first column
	obj2=cbind(rownames(geno),obj2)
	obj3=cbind(rownames(geno),obj3)
	# Transpose both matrices
	obj2=t(obj2)
	obj3=t(obj3)
	# Merge both matrices and put them in order according to the name of the individuals
	obj4=cbind(obj2, obj3)
	obj4=obj4[,order(obj4[1,])]
	
	# Add two columns: the identifier ("I") / marker data ("M"), and the marker names ("id")
	obj5=cbind(c("I", rep("M",ncol(geno))), c("id",colnames(geno)), obj4)

	write.table(obj5, file="./Results_Beagle/Matrix_BEAGLE.bgl", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


	### Call BEAGLE (version 3) from R and infer the missing data
	command.beagle = paste0("java -Xmx", Mb,"m -jar beagle.jar unphased=./Results_Beagle/Matrix_BEAGLE.bgl missing=? out=./Results_Beagle/Output")
	command.beagle = paste0("java -Xmx", Mb,"m -jar ./PROGRAM/beagle.jar unphased=./Results_Beagle/Matrix_BEAGLE.bgl missing=? out=./Results_Beagle/Output")
	
	# the Xmx option allows Beagle to use a certain amount of memory.
	
	if(Sys.info()[["sysname"]] == "Linux"){system(command.beagle, intern=FALSE)}
	if(Sys.info()[["sysname"]] == "Windows"){shell(command.beagle)}
	if(Sys.info()[["sysname"]] == "Darwin"){system(command.beagle)}

	rm(obj2,obj3,obj4,obj5)


	#################################################################
	## Convert Beagle genotyping output file in matrix with row=ind and col=SNPs ##
	#################################################################
	mat = read.table(file="./Results_Beagle/Output.Matrix_BEAGLE.bgl.phased.gz", header=TRUE, check.names=FALSE, sep="", stringsAsFactors=FALSE)
	mat = as.matrix(mat)

	markers = mat[,2]
	ind = unique(colnames(mat))
	ind = ind[-(1:2)]
	mat = mat[,-(1:2)]

	n  =  nrow(mat)
	m = ncol(mat)/2
	res1  =  matrix(0, n, m)
	res2  =  matrix(0, n, m)
	res3  =  matrix(0, n, m)
  

	conc  <-  function(mat){
		for(i in 1:m){res1[,i] = mat[,(2*i-1)]}
  
  		for(i in 1:m){res2[,i] = mat[,(2*i)]}
  
		for(i in 1:m){res3[,i] = paste( res1[,i], res2[,i], sep = "" )}
		return(res3)
	}

	new.mat=conc(mat)
	rm(mat, res1, res2, res3)
	rownames(new.mat)=markers
	colnames(new.mat)=ind
	new.mat=t(new.mat)

	# Reorder the individual names as it was initially in geno
	geno=new.mat[match(rownames(geno), rownames(new.mat)), ]
	rm(new.mat)
	# Save the genotyping file with the missing data imputed
	save(geno, file="./Results_Beagle/Matrix_DMimputedBeagle.RData")

}




