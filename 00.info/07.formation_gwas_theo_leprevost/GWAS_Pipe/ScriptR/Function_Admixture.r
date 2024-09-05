#############################################################
## Function: Population Structure
#############################################################
## Author: Sandra Negro
#############################################################
# Date: 30 September 2015

## AIM: Identification of sub-groups by Admixture (usually using PANZEA markers)
# Estimating the ancestral factions and the allelic frequencies of the population using Admixture

## ATTENTION: the operating system must be Linux. The program do not run on Windows
# !!! Admxiture software must be install!!! 

## INPUT files:
# (1) "geno" = genotyping matrix as outputted in the folder "Results_plink" (ex. Matrix_247_4491.RData")
#
# (2) "path" = path to the Admixture software. If you do not know, open a terminal and type the command line "which admixture". You will have the path
#
# (3) "Q" = Serie of sub-groups to test (ex.: Q = "1 2 3 4 5 6")
#
# (4) "Qmax" = is the last number of the serie "Q" (ex.: Qmax = 6)
#
# (5) "nbQ" = number of sub-groups we think there is for the stacked bar plot (ex.: nbQ=4)
#
# (6) "Q.threshold" = Ancestral fraction threshold for the stacked bar plot (ex.: Q.threshold=0.3)


## OUTPUT files:
# (1) Creation of the "Results_Admixture" folder

# (2) matrices with the ancestral fractions (probabities of affiliation to a group); each line = an individual (inbred line) and each column = a specific group, the somme of the column equal to one (the somme of the ancestral fractions for all group egal to one)

# (3) matrices with the allelic frequencies; one row = a marker and one column = one sub-group.
# (4) Stacked bar showing the ancestral fractions (Y-axis) of each individual (X-axis)

## REFERENCE
# D.H. Alexander, J. Novembre, and K. Lange.  Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655-1664, 2009

###########################################################

admixture <- function(path, Q, Qmax, geno){

	dir.create("Results_Admixture")
	setwd("./Results_Admixture")

	### Call and run Admixture software. Here you know how many subgroups you have e.g. 5
	#command.adm = paste(path, paste0("../Results_plink/", bed), 5)
	#system(command.adm, intern=FALSE) 


	### To have the standard errors (optional), add the "-B" flag
	#command.adm = paste(path, "-B", paste0("../Results_plink/", bed), 5)
	#system(command.adm, intern=FALSE)


	### !!! If the analyses are slow because of a large number of sub-groups or a high number of markers, then run admixture with a multithreaded mode. If you have four processors
	#command.adm = paste(path, paste0("../Results_plink/", bed), 5, "-j4") 
	#system(command.adm, intern=FALSE) 


	### If you do not have a clue of how many subgroups you might have then use the following command line and add the "--cv" flag for the method of cross-validation of Admixture. !!! use PANZEA markers only if you have the 50K chip
	#binary.name = paste("../Results_plink/Matrix", dim(geno)[1], dim(geno)[2], sep="_")
	binary.name = c("../Results_plink/MatrixGeno")
	command.adm = paste0("for K in ", Q,"; ", "\ do ", path, " --cv ", binary.name, ".bed", " $K | tee log${K}.out; done")
	
	# Here, you will have the values of the CV (from 1 to 10 sub-groups) that you can put on a graph: the K on the x-axis and the CV errors on the y-axis
	system(command.adm, intern=FALSE)

	### Add the row names to the "Q." files
	#adm.name = paste("Matrix", dim(geno)[1], dim(geno)[2], sep="_")
	adm.name = "MatrixGeno"
	
	for(i in 1:Qmax){
	adm=read.table(file=paste(adm.name, i,"Q", sep="."),sep=" ", header=FALSE, stringsAsFactors=FALSE)

	rownames(adm) = rownames(geno)
	write.table(adm, file=paste(adm.name, i,"Q", sep="."),sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
	}
	setwd("..")
}

	
	##################################
	## Stacked Bar
	##################################

StackedBar <- function(geno, nbQ, Q.threshold){

	adm.name = paste("./Results_Admixture/MatrixGeno")

	
	adm.nbQ=read.table(file=paste(adm.name, nbQ,"Q", sep="."),sep="\t", header=FALSE, row.names=1, stringsAsFactors=FALSE)

	### Graph
	subset=list()
	ord=list()
	ordered=NULL
	for(i in 1:nbQ){
		sub=subset(adm.nbQ, adm.nbQ[,i] >= Q.threshold)
		subset[[i]]=sub

		ord[[i]]=subset[[i]][order(subset[[i]][,1], decreasing=T), ]
		ordered=rbind(ordered, ord[[i]])
	}
	ordered=as.matrix(ordered)


	jpeg(paste0("./Results_Admixture/StackedBar_", nbQ, "Q.jpg"))

	par(las=2, cex.axis=0.4, mar=c(7,6,2,1), pin=c(width=5.7, height=1))

	#par(las=2, cex.axis=0.4, mar=c(7,4,2,1), pin=c(width=9.7, height=1.3))
	barplot(t(ordered), col=rainbow(nbQ), names.arg=as.vector(rownames(ordered)), ylab="Ancestry", border=NA)

	dev.off()

}

