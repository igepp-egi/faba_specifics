######################################################################################
######################################################################################
## GWAS with no corrections in the statistical model
######################################################################################
## Author: Sandra Negro
######################################################################################
## DATE: 13 October 2015


## AIM: GWAS with no correction (no kinship, no genetic structure in the model) 
#
## INPUT:
# (1) pheno = File containing the phenotypes to be considered. Choose the variables to analyse (the columns contains different sites by trait). First column, with header named "variety", contains the individual ID and the following columns include the phenotype with header names such as "Site.Year.Treatment.Trait" (ex.: biogai.2012.drought.plantHT) so that the column names are self informative. Use exactly this format for the column names.
#
# (2) geno = Genotyping matrix without NA. Observations (i.e. individuals) in lines and variables (i.e. SNPs) in columns. Genotypes as dose code; 0, 1, 2 where 0 and 2 are homozygotes for the allele 1 and the allele 2, respectively. Rq: In the function, the code will be divided in 2 to have frequency code; 0, 0.5, 1.
#
# (3) FASTLMM folder should be in the working directory.

## OUTPUT:
# A list "List_Pvalue1.RData" with the pvalue per markers for each situation and each trait analysed.
########################################################################################

model1 <- function(geno, pheno){
	dir.create("Results_GWAS_No_Corrections")

	pheno = as.data.frame(pheno)
	geno = geno[rownames(geno) %in% pheno$variety, ]
	pheno = pheno[match(rownames(geno), pheno$variety), ]
	geno = geno/2

	list.pval1 = list()
	for(j in 2:length(names(pheno))){
		Pvalue1 = NULL
		# ATTENTION: the phenotype might have missing data, hence when we remove the individuals with NA. By doing so, some SNPs in the geno become monomorphic. Hence, we have to remove those SNPs before the analyses.
		# Load the fastlmm results to extract the polymorphic SNPs
		res = read.table(paste0("./FASTLMM/Results_FASTLMM/RESULTS_FASTLMM_", names(pheno)[j], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
		
		phenotype = pheno[-which(is.na(pheno[j])), ]
	
		# Extract the polymorphic SNPs and the individuals having the phenotype
		genotype = geno[rownames(geno) %in% phenotype$variety, colnames(geno) %in% res$SNP]
		
	
		### Model without correction (no kinship, no genetic structure in the model)
		for(i in 1:ncol(genotype)) {
			model=lm(as.numeric(phenotype[,j])~genotype[,i])
			Pvalue1[i]=anova(model)[1,5]
			print(i)
		}
		Pvalue1 = as.data.frame(Pvalue1)
		rownames(Pvalue1) = colnames(genotype)
		colnames(Pvalue1) = names(phenotype)[j]
		list.pval1[[j-1]] = Pvalue1
	}

	save(list.pval1, file="./Results_GWAS_No_Corrections/List_Pvalue1.RData")
}


