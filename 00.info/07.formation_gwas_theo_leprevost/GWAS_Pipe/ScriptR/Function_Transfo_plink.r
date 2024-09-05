#############################################################
## Function: conversion to plink binary files
#############################################################
## Author: Sandra Negro
#############################################################
# Date: 30 Septembre 2015

### AIM: 
# Creation of plink binary files for the whole genome and per chromosome i.e. .bed/.bim/.fam
#
### ATTENTION: the plink software should be install and for Windows the executable file (plink.exe) should be in the working directory !!!


### INPUT:
# (1) "geno" = Genotyping matrix (missing values (NA) are allowed), genotyping code is lettres (ex. AG, CT) or any numbers (ex. 11, 12, 22), markers should be biallelic; the markers are in columns and the observations (individuals) are in rows.
#
# (2) "snp" = The SNP information file with details on the markers i.e. 1st column has the column name "SNP.name", 2nd column has the column name "Chromosome" (chromosome number), 3rd column has the column name "Position" (physical position).
#
# ATTENTION: The unplaced SNP should be in a fictive chromosome (i.e. unplaced SNPs in a chr11 for maize).
#
# (3) "chr" = number of chromosome present in the SNP information file (e.i. for maize chr=10 or chr=11 if there are unplaced SNPs)
#

### OUTPUT:
# (1) Creation of a "Results_plink" folder
#
# (2) "SNP_Info_UniquePosition_Ordered_(SNP.number).csv" = Detailed marker file where all SNPs have an unique physical position. SNPs are in (increasing) order by chromosome and position.  
#
# (3) "Matrix_(Ind.number)_(SNP.number).RData" = A new genotyping matrix with deletion of SNPs that have duplicated SNP positions. Markers are in increasing order of the genome.
#
# (4) "Liste_SNP_elimines_(Date).RData" = A list-object including the deleted SNPs that had duplicated SNP positions.
#
# (5) Three type of output binary files (.bed/.bim/.fam) for the whole set of SNPs and for each chromosome. The .ped and .map formats are also outputted.


#############################################################

Conversion <- function(snp, geno, chr){
	geno = as.matrix(geno)
	snp = as.data.frame(snp)
	snp$Chromosome=as.integer(snp$Chromosome)
	snp$Position=as.integer(snp$Position)

	# create the following result folder in the working directory
	dir.create("Results_plink")
	

	####################
	### .map format ###
	####################
	# This format contains four columns: chromosome number, rs# or snp identifier, genetic distance (morgans), base-pair position (bp units). Missing value in the .map format is 0


	### Extract the information of the markers (chr number and position) present in the genotyping file
	snp.info = snp[snp$SNP.name %in% colnames(geno), ]
	### Put the SNP in the order as in snp.info 
	geno = geno[, match(snp.info$SNP.name, colnames(geno))]

	

	### Check for the unique position of each SNP on its corresponding chromosome
	liste.elimines = list()
	elimines = NULL
	for(i in 1:chr){
  		sub = subset(snp.info, snp.info$Chromosome==i)

  		if(dim(sub)[1]-(length(unique(sub$Position)))>0){
    			print(paste("There are", dim(sub)[1]-(length(unique(sub$Position))), "duplicated marker positions on the chromosome", i, sep=" "))
    			unique = sub[which(duplicated(sub$Position) == FALSE), ]
    			elimines = sub[which(duplicated(sub$Position) == TRUE), ]
    
		} else {print(paste("The marker positions are all unique on the chromosome", i, sep=" ")); unique = sub}
  
		liste.elimines[[i]] = elimines
  
		if(dim(sub)[1]!=dim(unique)[1]){ 
			geno = geno[, which(colnames(geno) %in% liste.elimines[[i]]$SNP.name == FALSE)]}
	}


	### Save the files
	save(liste.elimines, file=paste0("./Results_plink/Liste_SNP_elimines_", as.character(Sys.Date()), ".RData"))
	snp.info = snp.info[snp.info$SNP.name %in% colnames(geno), ]
	# Check that the markers are in order in term of chromo and position
	snp.info = snp.info[order(snp.info$Chromosome, snp.info$Position, decreasing=FALSE), ]
	# Put the markers of the genotyping file in order
	geno = geno[, match(snp.info$SNP.name, colnames(geno))]
	
	#write.table(snp.info, file=paste0("./Results_plink/SNP_Info_UniquePosition_Ordered_", dim(geno)[2],".csv"), sep=",", col.names=TRUE, quote=FALSE)
	#save(geno, file=paste0("./Results_plink/Matrix_", dim(geno)[1], "_", dim(geno)[2], ".RData"))

	write.table(snp.info, file=paste0("./Results_plink/SNP_Info_UniquePosition_Ordered.csv"), sep=",", col.names=TRUE, quote=FALSE)
	save(geno, file=paste0("./Results_plink/MatrixGeno.RData"))
	

	### Build the .map format
	map = cbind(snp.info$Chromosome, snp.info$SNP.name, rep(0,nrow(snp.info)), snp.info$Position)
	# Save the map file
	#write.table(map, file=paste0("./Results_plink/Matrix_", dim(geno)[1], "_", dim(geno)[2], ".map"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	write.table(map, file=paste0("./Results_plink/MatrixGeno.map"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	


	#####################
	#### .ped format ###
	#####################

	# The .ped format has six mandatory columns: 1) family ID, 2) Individual ID, 3) Paternal ID, 4) Maternal ID, 5) sex (0 for unknown), 6) phenotype. For the 7th column and on, there are the genotypes (all markers should be biallelic, alleles are white-space delimited, e.g. "A A"). plink runs on Linux and runs on Windows as well. On Windows the executable file ".exe" must be in the same directory than the input files.


	### This function will work if you have missing data in the genotyping file
	if(sum(is.na(geno)>0)){
		rempl = function(dataset){
			pos=which(is.na(dataset))
			dataset[pos]="00"
			return(dataset)
		}
		geno=apply(geno,2,rempl)
	}



	### Build the .ped format
	# If the genotyping code is a lettre/number code (AGCT or 1234 for ex.) WITHOUT SPACE between alleles (e.g. "AA"), there is an option in plink ("--compound-genotypes") that will do the white-space delimited genotype (e.g. "A A")

	ped = cbind(seq(1,nrow(geno)), rownames(geno), rep(0,nrow(geno)), rep(0,nrow(geno)), rep(0,nrow(geno)), rep(0,nrow(geno)), geno)
	# Save the ped file
	#write.table(ped, file=paste0("./Results_plink/Matrix_", dim(geno)[1], "_", dim(geno)[2], ".ped"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	write.table(ped, file=paste0("./Results_plink/MatrixGeno.ped"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	


	### Calling plink from R to generate binary files (.bed, .bim and .fam)

	#file.name = paste0("./Results_plink/Matrix_", dim(geno)[1], "_", dim(geno)[2])
	#out.name = paste0("./Results_plink/Matrix_", dim(geno)[1], "_", dim(geno)[2])
	file.name = paste0("./Results_plink/MatrixGeno")
	out.name = paste0("./Results_plink/MatrixGeno")
	
	if(Sys.info()[["sysname"]] == "Linux"){
	command.plink = paste0("p-link --noweb --file ", file.name, " --compound-genotypes --make-bed --out ", out.name)
	# This will produce directly the .bed, .bim and .fam
	# !!! the file name must be the same for .ped and .map, otherwise specify the type of file (ex.: p-link --ped mydata.ped --map autosomal.map --compound-genotypes --make-bed --out outputmydata)
	system(command.plink, intern=FALSE)}
	
	if(Sys.info()[["sysname"]] == "Windows"){
	command.plink = paste0("plink --noweb --file ", file.name, " --compound-genotypes --make-bed --out ", out.name)
	shell(command.plink)}





	###############################################
	## For performing Rincent et al. 2014 approach
	## Create a genotyping matrix per chromosome
	###############################################

	# create .map
	map.chr=list()
	for(i in 1:chr){
		map.chr[[i]] = subset(map, as.numeric(map[,1])==i)
		# Save the map file
		#write.table(map.chr[[i]], file=paste0("./Results_plink/Matrix", "_", dim(geno)[1], "_", dim(geno)[2], "_",  "Chr", i, ".map"), sep=" ", col.names=F, row.names=F, quote=F)
		write.table(map.chr[[i]], file=paste0("./Results_plink/MatrixGeno_Chr", i,".map"), sep=" ", col.names=F, row.names=F, quote=F)
		
	}

	# create .ped
	geno.chr=list()
	for(i in 1:chr){
		geno.chr[[i]] = geno[, match(map.chr[[i]][,2], colnames(geno))]

		ped = cbind(seq(1,nrow(geno.chr[[i]])), rownames(geno.chr[[i]]), rep(0,nrow(geno.chr[[i]])), rep(0,nrow(geno.chr[[i]])), rep(0,nrow(geno.chr[[i]])), rep(0,nrow(geno.chr[[i]])), geno.chr[[i]])

		# Save the ped file
		#write.table(ped, file=paste0("./Results_plink/Matrix", "_", dim(geno)[1], "_", dim(geno)[2], "_", "Chr", i, ".ped"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
		write.table(ped, file=paste0("./Results_plink/MatrixGeno_Chr", i, ".ped"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	}


	# Call p-link from R
	for(i in 1:chr){
		if(Sys.info()[["sysname"]] == "Linux"){
			command.plink = paste0("p-link --noweb --file ", file.name, "_Chr", i, " --compound-genotypes --make-bed --out ", out.name, "_Chr", i)

			system(command.plink, intern=FALSE)}

		if(Sys.info()[["sysname"]] == "Windows"){
			command.plink = paste0("plink --noweb --file ", file.name, "_Chr", i, " --compound-genotypes --make-bed --out ", out.name, "_Chr", i)

			shell(command.plink)}
	}

	
}














