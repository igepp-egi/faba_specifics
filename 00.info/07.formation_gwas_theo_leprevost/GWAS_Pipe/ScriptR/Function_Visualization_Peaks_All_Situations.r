##############################################################################
##############################################################################
## Function: Visualization of the association peaks for all the situations by trait ##
##############################################################################
##############################################################################
## Author: Sandra Negro
##############################################################################
## Date: 5 November 2015

## AIM: Visualization of the significant SNPs along the genome for all situations by trait

## INPUT:
# (1) "snp" = SNP information file, which contains the SNP name, chromosome nb and position along the genome. The column names are as followed: "SNP.name", "Chromosome", "Position".
# (2) "ChrStart" = the first chromosome that you would like to visualized on the figure (ex. if you would like to see the significant SNPs from chr1 to chr4 then ChrStart=1)
# (3) "ChrEnd" = the last chromosome that you would like to visualized on the figure (ex. if you would like to see the significant SNPs from chr1 to chr4 then ChrEnd=4)
# (4) "n_effective" =  number of indepedant tests (number of independant markers).
# (5) "alpha" = type I error rate (e.g. 0.05)
# (6) "pheno" = File containing the phenotypes. First column, must have column name "variety" and contains the individual ID and the following columns include the phenotype with header names such as "Site.Year.Treatment.Trait" (ex.: "biogai.2012.drought.plantHT") so that the column names are self informative. Use exactly this format for the column names. If you do not have year or treatments write for ex. "biogai.NoYear.NoTreatment.plantHT".
# (7) "trait" = choose the trait and it must be exactly written like in the columns of pheno ex. "plantHT"
# (8) The folder FASTLMM should be in the working directory. It contains the FASTLMM output results.
# (9) "cex" = numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
# (10) "title" = title of the graph (ex: "Significant SNPs from GWAS using the 50K chip for Plant Height trait")
# (11) "width" = the width of the graph in inches (ex: width=14)

## OUTPUT:
# (1) A graph named "Visualization_SignSNPs_All_Situations_Chr ChrStart-ChrEnd_trait.pdf": visualization of the significant SNPs for all situations (one situation = one site, one year, one treatment) and by trait along the genome.

##############################################################################
peaks.all.situations <- function(snp, ChrStart, ChrEnd, n_effective, alpha, pheno, trait, cex=0.9, title, width=7){
	snp=as.data.frame(snp)
	pheno=as.data.frame(pheno)
	

		### Add a column "Cumul" to the snp object. Cumul will be use to build the x-axis (positin along the genome) of the graph
	
	# Extract a subset of the SNP info from the first chromosome 
	snp.cumul = subset(snp, snp$Chromosome == ChrStart)
	Cumul = snp.cumul$Position
	# Add "Cumul" column to the subset of SNP info
	snp.cumul = cbind(snp.cumul, Cumul)
	# Last "Cumul" position if the first chromosome
	EndChr = snp.cumul$Cumul[nrow(snp.cumul)]
	# Positions on the x-axis of the first labels (e.g. Chr1)
	pos.labels=snp.cumul$Cumul[nrow(snp.cumul)]/2
	# Add the Cumul column to the other chromosomes 
	for(i in (ChrStart+1):ChrEnd){
		sub = subset(snp, snp$Chromosome == i)
		Cumul = EndChr+sub$Position
		sub = cbind(sub, Cumul)
		snp.cumul = rbind(snp.cumul, sub)
		pos = EndChr + (snp.cumul$Cumul[nrow(snp.cumul)] - EndChr)/2
		pos.labels=c(pos.labels, pos)

		EndChr = snp.cumul$Cumul[nrow(snp.cumul)]
		
	}

	
		### Extract the significant SNPs for all the situations by trait
	sign.res = list()
	situation = NULL
	for(i in 2:length(names(pheno))){
	
		Trait = unlist(strsplit(names(pheno)[i], ".", fixed=TRUE))
		if(Trait[4] == trait){
			### Load the fastlmm output results
			res = read.table(paste0("./FASTLMM/Results_FASTLMM/RESULTS_FASTLMM_", names(pheno)[i], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
			#res = read.table(file=paste0("FASTLMM_COMBINES_AvecLigneeIllegitimes_AvecPriseEnCptLigneesPuresDans600K_UniquePosition/RESULTS_FASTLMM_COMBINES_", names(pheno)[i], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
			
			res = res[which(res$Pvalue <= alpha/n_effective), ]
			res = res[which(res$Chromosome >= ChrStart), ]
			res = res[which(res$Chromosome <= ChrEnd), ]
		
			sign.res[[i]] = res
			situation = c(situation, names(pheno)[i])	
		}
	}
	# Extract from the list only the results from the trait of interest	
	sign.res = sign.res[which(names(pheno) %in% situation)]
	

	
		### Create the graph with all the sitations by trait
	pdf(paste0("Visualization_SignSNPs_All_Situations_Chr",ChrStart, "-", ChrEnd, "_",trait , ".pdf"), width=12)	
	par(mar=c(5,12,4,1))
	plot(snp.cumul$Cumul, rep(0,nrow(snp.cumul)),ylim=c(0,length(sign.res)), xlim=c(snp.cumul$Cumul[1], snp.cumul$Cumul[nrow(snp.cumul)]) , col="white", ylab="", xlab=paste("Chromosome:", ChrStart, "to", ChrEnd), cex=cex, xaxt="n", yaxt="n", cex.lab=1.2, main=title)
	axis(2, at=seq(1,length(sign.res)),labels=situation, cex.axis=0.8 ,las=1)
	chr.label = paste0("Chr", ChrStart:ChrEnd)
	axis(1, at=pos.labels , labels=chr.label, tck=0)
	
	EndChr=NULL
	for(i in ChrStart:ChrEnd){
		sub = subset(snp.cumul, snp.cumul$Chromosome == i)
		EndChr=sub$Cumul[nrow(sub)]
		abline(v=EndChr, col="red")
	}
	
	for(i in 1:length(sign.res)){
		xpos = snp.cumul$Cumul[snp.cumul$SNP.name %in% sign.res[[i]]$SNP]
		points(xpos, rep(i, length(xpos)), pch=1)
		abline(h=i)
	}
	dev.off()
}

