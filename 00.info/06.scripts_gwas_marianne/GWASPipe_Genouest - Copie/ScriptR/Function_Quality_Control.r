######################################################
## Descriptive analyses and quality control
#####################################################
#####################################################
## Author: Sandra Negro
####################################################
# Date: March 2016

### README !!!
## AIM: 
# Descriptive analyses of the genotyping dataset 
# Applying filters for the quality control and extraction of the selected markers and taxa

## INPUT:
# (1) "geno" = Genotyping matrix (missing values (NA) are allowed), the genotyping code is lettres (ex. AG, CT, etc...), the markers should be biallelic; the markers are in columns and the observations (individuals) are in rows.
#
# (2) "snp" = The SNP information file with details on the markers i.e. 1st column has to have the column name "SNP.name", 2nd column has to have the column name "Chromosome" (chromosome number), 3rd column has to have the column name "Position" (physical position).

# Rq: ATTENTION IF YOU HAVE UNPLACED SNPs NOTED OTHER THAN "NA" (ex. "UNK") PLEASE REPLACE IT BY "NA".

# (3) "input.geno" = Path to the genotyping file to be written in the log record output file
# (4) "input.snp" = Path to the SNP information file to be written in the log record output file
# (5) "Type.analysis" = Type of analysis to be written in the log file (ex: "Quality Control")
#
# (6) Parameters:
# "snp.dm" = threshold (proportion) of SNP missing data ex.: snp.dm = 0.20
# "maf" = threshold for  the frequency of the minor allele (MAF) ex.: maf = 0.01
# "snp.het" = threshold of SNP heterozygosity ex.: snp.het = 0.15
# "ind.dm" = threshold of Ind missing data ex.: ind.dm = 0.10
# "ind.het" = threshold of Ind heterozygosity ex.: ind.het = 0.08

## OUTPUT:
# (1) Creation of a "Results_QC" folder
#
# (2) If there are genotypes other than AA, AC, AT, AG, CC, etc... the "abnormal" genotypes and their indices will be written in a file named "AbnormalGenotypes.txt". ATTENTION: if there are indels, they will be saved in this file.
#
# (3) "Liste_SNP_Duplicated_Position.RData" = If you have duplicated SNP positions, then a list-object including the deleted SNPs that had duplicated SNP positions per chromosome will be outputted.

# (4) "SNP_Info_UniquePosition_Ordered_(SNP.number).RData" = Detailed marker file where the unplaced SNP are put into a fictive chromosome (ex. in maize chr11) with arbitrary physical position, and where all SNPs have an unique physical position. SNPs are in (increasing) order by chromosome and position. 
#
# (5) Tables with counts and rates by marker and by individual: "Effectif_SNP_Diversite.RData", "AfterQC_SNP_Effective_Diversite.RData", "Effectif_Taxa.RData", "AfterQC_Taxa_Effective.RData"
#
# (6) Several graphs on the distribution of missing data, heterozygosity and MAF: "Hist_SNP_Het_MAF_DM.pdf", "Plot_Ind_Het_DM.pdf"
#
# (7) "Matrix_QCed.RData" = genotyping matrix after quality control
#
# (8) "Log.txt" = A log record file
####################################################
QC <- function(geno, snp, input.geno, input.snp, Type.analysis, snp.dm, maf, snp.het, ind.dm, ind.het){
	
	dir.create("Results_QC")
	geno=as.matrix(geno)
	snp=as.data.frame(snp)

	
	### Check for abnormal genotypes
	# Indices of abnormal genotypes
	homozygote <- c('AA','CC','GG','TT')
    heterozygote <- c("AT","AC","AG","CG","TC","TG","TA","CA","GA","GC","CT","GT")	
    pos =which(Reduce('+',lapply(c(homozygote,heterozygote), function(x) geno == x))==0,arr.ind=TRUE)
		
	# Store in a vector the type of abnormal genotypes (ATTENTION: these genotypes might be insertions and deletions)
    if(nrow(pos)>0){
        # Build a table including indices of each abnormal genotypes
        abn.geno = cbind(pos, as.vector(geno[pos]))
        write.table(abn.geno, file="./Results_QC/AbnormalGenotypes.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
        # Replace the abnormal genotypes by "NA" and save the table
        geno[pos] = NA
    }
	# Some markers might have abnormal (NA) genotypes for all individuals, remove these SNPs
	# ATTENTION: these markers might have insertions and deletions
    pos = sapply(1:ncol(geno), function(x) any(!is.na(geno[,x])))    
    geno = geno[, pos]

	

	### Check if the names of the columns and rows are unique
	# How many duplicates in rownames and colnames
	dupl.row=sum(duplicated(rownames(geno)))
	dupl.col=sum(duplicated(colnames(geno)))
	cat("\n")
	print("########################################")
	print(paste("There are", dupl.row, "duplicated row names"))
	print(paste("There are", dupl.col, "duplicated column names"))
	print("########################################")
	cat("\n")
	geno = geno[duplicated(rownames(geno)) == FALSE, duplicated(colnames(geno)) == FALSE]

	

	
	### If there are unplaced SNP, assign a chromo "Nb of chr+1" with arbitrary position
    na.chr.or.pos <- which(is.na(snp$Chromosome)|is.na(snp$Position))
    nb.na <- length(na.chr.or.pos)
    if (nb.na>0){
        snp$Chromosome[na.chr.or.pos] <- length(unique(snp$Chromosome)) # "unique()" counts the "NA", so if there are "NA", "length(unique(snp$Chromosome))=11 (for maize)"
        snp$Position[na.chr.or.pos] <- 1:nb.na
    }
	
	

	### Check for the unique position of each SNP on its corresponding chromosome
	liste.elimines = list()
	elimines = NULL
	snp = snp[snp$SNP.name %in% colnames(geno),]
	for(i in 1:length(unique(snp$Chromosome))){
  		sub = subset(snp, snp$Chromosome==i)

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
	if(length(liste.elimines) > 0){save(liste.elimines, file="./Results_QC/Liste_SNP_Duplicated_Position.RData")}
	
	# Check that the markers are in order in term of chromo and position
	snp = snp[order(snp$Chromosome, snp$Position, decreasing=FALSE), ]

	# Put the markers of the genotyping file in order
	geno.ordered = geno[, match(snp$SNP.name, colnames(geno))]

	#save(snp, file=paste0("./Results_QC/SNP_Info_UniquePosition_Ordered_", dim(geno)[2],".RData"))
	save(snp, file=paste0("./Results_QC/SNP_Info_UniquePosition_Ordered.RData"))
	rm(geno)




	############################## 
	### Count on the columns (SNP)
	##############################
	# List of the type of polymorphism of the SNPs
	list.polymorphism <- sort(unique(as.vector(geno.ordered)))
	# Count of the nb of the different type of polymorphism (AA, AC, AG...)
	Effective_SNP <- sapply(list.polymorphism, function(z) apply(geno.ordered, 2, function(x) length(which(x==z))))
	Effective_SNP = as.data.frame(Effective_SNP)
	# Count of the missing data per marker
	Effective_SNP$DM = apply(geno.ordered, 2, function(x) sum(is.na(x)))

	

	### Frequencies by nucleotide and by SNP
	# List of the type of alleles of the SNPs
	list.alleles = sort(unique(unlist(strsplit(geno.ordered,'')) ))

	
	# Count of the nb of T, A, C, G by column
	Res <-sapply(list.alleles, function(z) apply(geno.ordered,2,function(x) sum(unlist(strsplit(x,''))==z, na.rm=TRUE)))
	Effective_SNP=cbind(Effective_SNP, Res)
	




	### MAF for the biallelic polymorphic loci
	col=paste0("P", colnames(Res))
	Res.P = sapply(1:ncol(Res), function(z) sapply(Res[,z], function(x) {x[which(x==0, arr.ind=TRUE)]=NA; return(x)}))
	colnames(Res.P) = col

	Effective_SNP = cbind(Effective_SNP, Res.P)

	Effective_SNP$MAF <- apply(Res.P, 1, function(x) min(x, na.rm=TRUE)/sum(x, na.rm=TRUE))

	Effective_SNP$MAF[which(Effective_SNP$MAF==1,arr.ind=TRUE)]=NA
	
	

	### Heterozygosity rates
	het = which(list.polymorphism %in% c("AT","TA","AC","CA","AG","GA","CT","TC","CG","GC","GT","TG"))

	Effective_SNP$HET <-apply(Effective_SNP, 1, function(x) sum(x[c(het)])/sum(x[c(1:length(list.polymorphism))]) )

	Effective_SNP$HET[which(is.na(Effective_SNP$MAF)==TRUE)]=NA


	### Diversity estimate (Nei's indice) using all SNPs
	## ATTENTION: usually we estimate Nei's indice on the PANZEA markers only.
	Effective_SNP$Hs <- 2*Effective_SNP$MAF*(1-Effective_SNP$MAF)
	Effective_SNP$Hs[which(is.na(Effective_SNP$MAF)==TRUE)]=NA


	### Missing data rate
	attach(Effective_SNP)
	Effective_SNP$txDM <- DM/nrow(geno.ordered)
	Effective_SNP$txDM[which(is.na(Effective_SNP$MAF)==TRUE)]=NA



	###############################
	### Filter on the SNPs ####
	###############################

	### How many alleles per SNP to detect the triallelic SNPs
	NbAll = apply(Res, 1, function(x) sum(x!=0))

	# Print nb of tri-allelic SNPs
	print(paste("Number of SNPs with more than two alleles:", length(which(NbAll > 2))))
	# Print nb of monomorphic SNPs
	print(paste("Number of monomorphic SNPs:", length(which(NbAll == 1))))

	# add a column with the total nb of alleles per SNP
	Effective_SNP$NbAll=NbAll

	# Save the table on counts
	save(Effective_SNP, file="./Results_QC/Effectif_SNP_Diversite.RData")


	### monomorphic SNP
	LocMono = names(which(NbAll==1))  


	### SNP with a missing data threshold > "snp.dm" 
	LocDM = rownames(Effective_SNP)[Effective_SNP$txDM > snp.dm | (is.na(Effective_SNP$txDM)==TRUE)]  


	### SNP with MAF < "maf"
	LocMAF=rownames(Effective_SNP)[Effective_SNP$MAF < maf | (is.na(Effective_SNP$MAF)==TRUE)]
	  

	### SNP with Het > "snp.het"
	LocHet=rownames(Effective_SNP)[Effective_SNP$HET > snp.het | (is.na(Effective_SNP$HET)==TRUE)]
	


	### Apply the quality control filters
	QCed_Effective_SNP = Effective_SNP[(rownames(Effective_SNP) %in% LocMono)==FALSE & (rownames(Effective_SNP) %in% LocDM)==FALSE & (rownames(Effective_SNP) %in% LocMAF)==FALSE & (rownames(Effective_SNP) %in% LocHet)==FALSE, ]
	# Extract the biallelic markers only
	QCed_Effective_SNP = QCed_Effective_SNP[QCed_Effective_SNP$NbAll %in% 2, ]

	

	### Save QCed table with the counts
	save(QCed_Effective_SNP, file="./Results_QC/AfterQC_SNP_Effective_Diversite.RData")


	### Mean and Range to be written in the log file
	# Het
	mean.snp.het = mean(QCed_Effective_SNP$HET, na.rm=TRUE)
	range.snp.het = range(QCed_Effective_SNP$HET, na.rm=TRUE)
	# Hs
	mean.nei=mean(Effective_SNP$Hs, na.rm=TRUE)
	range.nei=range(Effective_SNP$Hs, na.rm=TRUE)
	# maf
	mean.maf = mean(QCed_Effective_SNP$MAF, na.rm=TRUE)
	range.maf = range(QCed_Effective_SNP$MAF, na.rm=TRUE)
	# DM
	mean.snp.dm = mean(QCed_Effective_SNP$txDM, na.rm=TRUE)
	range.snp.dm = range(QCed_Effective_SNP$txDM, na.rm=TRUE)



	#############################
	### Graphs
	############################
	pdf("./Results_QC/Hist_SNP_Het_MAF_DM.pdf")
	par(mfrow=c(2,2))
	# Heterozygosity
	hist(QCed_Effective_SNP$HET, main="SNP heterozygosity rate distribution", xlab="Heterozygosity Rate", cex.main=1, cex.sub=1.2)
	# MAF
	hist(QCed_Effective_SNP$MAF, main="SNP MAF Distribution", xlab="Minor Allele Frequency", cex.main=1, cex.sub=1.2)
	# Missing data
	hist(QCed_Effective_SNP$txDM, main="SNP missing data rate distribution", xlab="Missing Genotype Rate", cex.main=1, cex.sub=1.2)
	dev.off()


	### In the genotype matrix, delete the SNPs that do not pass the thresholds
	geno.pruned = geno.ordered[, match(rownames(QCed_Effective_SNP), colnames(geno.ordered))]

	rm(geno.ordered)


	######################################
	## Count on the rows (individuals)
	#####################################
	# Count of the nb of the different type of polymorphism (AA, AC, AG...)
	Effective_taxa <- sapply(list.polymorphism, function(z) apply(geno.pruned, 1, function(x) length(which(x==z))))

	Effective_taxa = as.data.frame(Effective_taxa)
	# Count of the missing data per individual
	Effective_taxa$DM = apply(geno.pruned, 1, function(x) sum(is.na(x)))



	### Frequencies by nucleotide
	# Count of the nb of T, A, C, G by row
	Res <-sapply(list.alleles, function(z) apply(geno.pruned,1,function(x) sum(unlist(strsplit(x,''))==z, na.rm=TRUE)))
	Effective_taxa=cbind(Effective_taxa, Res)



	### Heterozygosity rates
	Effective_taxa$HET <-apply(Effective_taxa, 1, function(x) sum(x[c(het)])/sum(x[c(1:length(list.polymorphism))]) )

	

	### Missing data rate
	attach(Effective_taxa)
	Effective_taxa$txDM <- DM/ncol(geno.pruned)


	### Save the table on individual counts
	save(Effective_taxa, file="./Results_QC/Effectif_Taxa.RData")



	### Ind. with a missing data rate < "ind.dm" 
	IndDM = rownames(Effective_taxa)[Effective_taxa$txDM > ind.dm | (is.na(Effective_taxa$txDM)==TRUE)]
	length(IndDM)  


	### Ind with Het > "ind.het"
	IndHet=rownames(Effective_taxa)[Effective_taxa$HET > ind.het | (is.na(Effective_taxa$HET)==TRUE)]
	length(IndHet) 


	### Apply the quality control filters
	QCed_Effective_taxa = Effective_taxa[(rownames(Effective_taxa) %in% IndDM)==FALSE & (rownames(Effective_taxa) %in% IndHet)==FALSE, ]

	

	### Save QCed table with the counts
	save(QCed_Effective_taxa, file="./Results_QC/AfterQC_Taxa_Effective.RData")


	### Save the genotyping file after QC
	geno.pruned.QC = geno.pruned[rownames(geno.pruned) %in% rownames(QCed_Effective_taxa), ]

	save(geno.pruned.QC, file="./Results_QC/Matrix_QCed.RData")


	### Mean and Range to be written in the log file
	# Het
	mean.ind.het = mean(QCed_Effective_taxa$HET, na.rm=TRUE)
	range.ind.het = range(QCed_Effective_taxa$HET, na.rm=TRUE)
	# DM
	mean.ind.dm = mean(QCed_Effective_taxa$txDM, na.rm=TRUE)
	range.ind.dm = range(QCed_Effective_taxa$txDM, na.rm=TRUE)





	############################
	### Graphique
	###########################
	pdf("./Results_QC/Plot_Ind_Het_DM.pdf")
	par(mfrow=c(2,2))
	# DM
	hist(QCed_Effective_taxa$txDM, main="Individuals missing data rate distribution", xlab="Missing Data Rate", cex.main=1, cex.sub=1.2)
	# Heterozygosity
	hist(QCed_Effective_taxa$HET, main="Individuals heterozygosity rate distribution", xlab="Heterozygosity Rate", cex.main=1, cex.sub=1.2)

	barplot(QCed_Effective_taxa$HET,horiz=TRUE, main="Individuals heterozygosity rate", ylab="lines", xlab="Heterozygoty rate", space=0, col="blue",add=FALSE,names.arg=rownames(QCed_Effective_taxa), las=1, cex.names=0.2, cex.main=1, cex.sub=1.2)
	dev.off()



	###########################
	### Write le log file
	###########################	 
	log = rbind(as.character(Sys.Date()), "INPUT FILES:", input.geno, input.snp, Type.analysis, "SNP Missing Data Rate Threshold:", snp.dm, "MAF threshold:", maf, "SNP Heteryzygosity Rate Threshold:", snp.het, "Individual Missing Data Rate Threshold:", ind.dm, "Individual Heterozygosity Rate Thresthold:", ind.het, "Mean SNP heterozygozity:", mean.snp.het, "Range SNP heterozygozity:", range.snp.het[1], range.snp.het[2], "Mean diversity (Nei's indice):", mean.nei, "Range diversity:", range.nei[1], range.nei[2], "Mean Minor Allele Frequency (MAF):", mean.maf, "Range MAF:", range.maf[1], range.maf[2], "Mean SNP missing data:", mean.snp.dm, "Range SNP missing data:", range.snp.dm[1], range.snp.dm[2], "Mean individual heterozygosity:", mean.ind.het, "Range individual heterozygosity:", range.ind.het[1], range.ind.het[2], "Mean individual missing data:", mean.ind.dm, "Range individual missing data:", range.ind.dm[1], range.ind.dm[2])

	write.table(log, file="./Results_QC/Log.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


}


