######################################################
## Function: built GNPasso outputs from FASTLMM results ##
######################################################
## Author: Sandra Negro
#####################################################
## DATE: February 2016

## AIM: building a GNPasso-like summary table that includes the main GWAS results
#
## INPUT: 
# (1) The folder FASTLMM should be in the working directory. It contains the FASTLMM output results.
# The Results_plink folder containing the binary files ".bim" should also be in the working directory (.bim contains the information on allele 1 (minor) and allele 2 (major)).
#
# (2) bim.all = the .bim file from the "Results_plink" folder that contains the SNPs from all chromosomes (ex.: "./Results_plink/Matrix_247_4778.bim"). 
#
# (3) pheno = the phenotype file: 1st column is the individual ID (the column name is "variety"), the following columns include the phenotype with column names such as "Site.Year.Treatment.Trait" (ex.: biogai.2012.drought.plantHT).
#
# (4) geno = the genotying file that include all the markers tested. The genotyping code can be letters (ex.: AG, CT, etc..) or numbers, but be biallelic (ex: 11, 12, 22). rows = individual ID and columns = marker names
#
# (5) User parameters:
# Experimenter = experimenter name
# Exper_Name = experiment name or project 
# Panel_Name = name of the panel
# Marker_Type = type of the marker (e.g. SNP, SSR, etc...)
# Software = software name (e.g. fastlmm, emmax, asreml)
# Statistical_Model = statistical model used for the GWAS (e.g. K (mixed model with a kinship), K+Q (kinship and genetic structure), etc...)
# RefGenMap = genome reference (e.g. Ref_v2)
# Nchr = number of chromosome (nb of chr present in the SNP info file; in maize, Nchr=10 or Nchr=11 if you have unplaced SNP put in a fictive chromosome)

## OUTPUT: 
# (1) the GNPasso-like GWAS template which contains the main GWAS results
# (2) Log record file containing information on the Experimenter name, date of the analyses, ref genome, nb of individuals used in the analyses, software and statistical model used.

#######################################################

gnpasso <- function(pheno, geno, bim.all, Experimenter, Exper_Name, Panel_Name, Marker_Type, Software, Statistical_Model, RefGenMap, Nchr){

	dir.create("Results_GNPasso")

	for(j in 2:ncol(pheno)){
		### Load the fastlmm output results
		res = read.table(paste0("./FASTLMM/Results_FASTLMM/RESULTS_FASTLMM_", names(pheno)[j], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
		

		## Correction for multiple tests
		# NB: for the FWER (Familly-Wise Error Rate), we seek to guarantee that the probability of having at least one false positive is < alpha
		# NB: for the FDR (False Discovery Rate), we seek to guarantee that the proportion of false positives among the declared positives is at most equal to alpha
		#Bonferroni = p.adjust(res$Pvalue, method="bonferroni") # FWER: Bonferroni
		#BH = p.adjust(res$Pvalue, method="BH") # FDR: Benjamini & Hochberg 1995

	

		# Extract and put the SNPs in the same order than the .bim file
 		bim.all.ordered = bim.all[bim.all[,2] %in% res$SNP, ]
		bim.all.ordered = bim.all.ordered[match(res$SNP, bim.all.ordered[,2]), ]
		# Minor allele
		Allele1 = bim.all.ordered[,5]
		# Major allele
		Allele2 = bim.all.ordered[,6]

		#######################################
		# Rq: Some individuals have missing data for the phenotype. It is possible that some SNPs with very low MAF become monomorphic due to deleting individuals with missing pheno. If we specify the missing phenotype value, fastlmm delete the SNPs that become monomorphic.


		##################################################################
		### Counts on the columns (SNP) of geno
		###########################################
		### Effectif allele 1 (Eff1) and effective allele 2 (Eff2) per situation
		##################################################################
		geno.eff = geno[rownames(geno) %in% pheno$variety[which(is.na(pheno[,j]))]==FALSE , colnames(geno) %in% res$SNP]
		geno.eff = geno.eff[, match(res$SNP, colnames(geno.eff))]

		# List of the type of polymorphism of the SNPs
		list.polymorphism <- sort(unique(as.vector(geno.eff)))
		### Count of the nb of the different type of polymorphism (AA, AC, AG...)
		Effective_SNP <- sapply(list.polymorphism, function(z) apply(geno.eff, 2, function(x) length(which(x==z))))
		Effective_SNP = as.data.frame(Effective_SNP)
		colnames(Effective_SNP)=list.polymorphism


		### Frequencies by nucleotide and by SNP
		# List of the type of alleles of the SNPs
		list.alleles = sort(unique(as.vector(apply(geno.eff, 2,function(x) unlist(strsplit(as.character(x),''))) )))

		# Count of the nb of T, A, C, G by column
		Res <-sapply(list.alleles, function(z) apply(geno.eff,2,function(x) sum(unlist(strsplit(as.character(x),''))==z, na.rm=TRUE)))
		Effective_SNP=cbind(Effective_SNP, Res)

		# Add allele 1 and allele 2
		Effective_SNP$Allele1 = Allele1
		Effective_SNP$Allele2 = Allele2
		# Add the count of allele 1 and allele 2
		Effective_SNP$Eff1 = rep(NA, nrow(Effective_SNP))
		Effective_SNP$Eff2 = rep(NA, nrow(Effective_SNP))

		for(i in 1:nrow(Effective_SNP)){
			for(x in colnames(Res)){
			if(Effective_SNP$Allele1[i] == x){Effective_SNP$Eff1[i] = Res[i,which(colnames(Res)==x)]}
			if(Effective_SNP$Allele2[i] == x){Effective_SNP$Eff2[i] = Res[i,which(colnames(Res)==x)]}
			}
		}


		Marker_Name = res$SNP
		Situation = names(pheno)[j]
		Chr_Nb = res$Chromosome
		Marker_Position = res$Position
		Test = res$WaldStat
		Pvalue = res$Pvalue
		LogPval = -log10(res$Pvalue)
		SNP_Weight = res$SNPWeight
		SNP_Weight_SE = res$SNPWeightSE
		# R2_LR = 1-e(-2*(log Lm-log L0)/n)
		Ind_Nb = res$N[1]
		R2_LR = 1 - exp(-2*(res$AltLogLike - res$NullLogLike)/Ind_Nb)
		trait <- unlist(strsplit(names(pheno)[j], ".", fixed=TRUE))
		Trait=trait[4] 


		### Build the GNPasso table
		template = cbind(Marker_Name, names(pheno)[j], Trait, res$Chromosome, Marker_Position, res$N, Test, Pvalue, LogPval, SNP_Weight, SNP_Weight_SE, R2_LR, Allele1, Allele2, Effective_SNP$Eff1, Effective_SNP$Eff2, res$GeneticVar, res$ResidualVar, res$NullGeneticVar, res$NullResidualVar)

		
		colnames(template) = c("Marker_Name", "Situation", "Trait", "Chromosome", "Marker_Position", "Ind_Nb", "Test", "Pvalue", "LogPval", "SNP_Weight", "SNP_Weight_SE", "R2_LR", "Allele1", "Allele2", "Eff1", "Eff2", "GeneticVar", "ResidualVar", "NullGeneticVar", "NullResidualVar")

		save(template , file=paste0("./Results_GNPasso/GNPasso_", names(pheno)[j], ".RData"))
	}
	
	### Save the log record file
	log = rbind(paste0("Experimenter: ",Experimenter) , as.character(Sys.Date()), "", Exper_Name, Panel_Name, paste0("Marker Type: ", Marker_Type), paste0("Genome Reference: ", RefGenMap), Software, paste0("Statistical model used: " ,Statistical_Model))
	write.table(log, file="./Results_GNPasso/Log.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	
}






