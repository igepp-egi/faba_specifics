###################################################################
###################################################################
### QQplots and Manhattan plots using the GWAS results from FASTLMM
###################################################################
## Authors: Stephane Nicolas and Sandra Negro
###################################################################
## DATE: December 2015

## AIM:
# Creation of QQplots and manhattan plots 

## INPUT
# (1) pheno = File containing the phenotypes. First column, with header named "variety", contains the individual ID and the following columns include the phenotype with header names such as "Site.Year.Treatment.Trait" (ex.: biogai.2012.drought.plantHT) so that the column names are self informative. Use exactly this format for the column names.
#
# (2) list.pval1 = a list object including the pvalues obtained from a statistical model with no corrections (no kinship and no genetic structure). In rownames there are the SNP names and in colnames there are the GWAS experiment. Each element of this list contains the pvalues for one site, one year, one treatment and one trait. This list object can be obtained with the function in the file "Function_GWAS_No_Corrections.r"
#
# (3) exp = list of the names of the experiments (=the column names of the pheno file e.g. "biogai.2012.drought.dta") grouped by trait and the site for ex. (e.g. the first element of the list object "exp" contains all experiments done in, for ex. "biogai" site for "DTA" trait: "biogai.2012.drought.dta", "biogai.2013.drought.dta", "biogai.2012.watered.dta", "biogai.2013.watered.dta"). !!! The length of the list "exp" (=the number of the elements of the list) must be equal to the column number of list.pval1 and in the same order of site and trait !!!
#
# (4) The folders "FASTLMM/Results_FASTLMM" with all the GWAS results outputted by FASTLMM should be in the working directory
#
# (5) snp = The SNP information file with details on the markers i.e. 1st column has the column name "SNP.name", 2nd column has the column name "Chromosome" (chromosome number), 3rd column has the column name "Position" (physical position). The class of this file is "data.frame"
# Rq: the markers should be in order in term of physical position on the genome
#
# (6) alpha =  the alpha error (type I)

## OUTPUT
# (1) Quantiles of the pvalues of each experiments (including the experiments of the statistical model without corrections
# (2) Quantiles-Quantiles plots (QQplots)
# (3) Manhattan plots pour toutes les situations par trait
################################################


###################################################################
## Creation of QQplots
###################################################################

#qqplot <- function(pheno, exp, table.pval1, alpha, gao, Bonferonni)
qqplot <- function(pheno, exp, list.pval1, alpha){
	dir.create("Results_QQplot_Manhattan")

	pheno = as.data.frame(pheno)
	
	######################################################################
	### Quantiles of the pvalues. Function from the "TD GWAS" from Stephane
	######################################################################
	CumSumLogLik=function(pvalt) {
		pvalt=pvalt[order(pvalt)]
		logObs=-log10(pvalt)
		logExp=-log10(seq(1,length(pvalt))/(length(pvalt)))
		Coord=cbind(logExp,logObs)
		return(Coord)
	}


	####################################################
	### Load the GWAS results, extract the pval and calculate the quantiles 
	####################################################
	
	for(j in 2:ncol(pheno)){
		
		res = read.table(paste0("FASTLMM/Results_FASTLMM/RESULTS_FASTLMM_", names(pheno)[j], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
		# Extract and combine the snp name and pval
		pval = cbind(res$SNP, res$Pvalue)
		colnames(pval) = c("SNP", "Pvalue")
		# Quantiles
		q_val = CumSumLogLik(as.numeric(pval[,2]))
		
		save(pval, q_val, file=paste0("Results_QQplot_Manhattan/Pval_and_q_val_", names(pheno)[j], ".RData"))
	}


	##################################################
	### Calculate the quantiles for model1 (statistical model without corrections)
	max1.all = NULL
	for(i in 1:length(list.pval1)){
		# Quantiles
		q_val1 = CumSumLogLik(as.numeric(list.pval1[[i]][,1]))
		# vector of the maximum quantile value per experiment
		max1.all = c(max1.all, max(q_val1))
		save(q_val1, file=paste0("Test_Result_julie/Results_GWAS_No_Corrections/q_val1_", colnames(list.pval1[[i]]), ".RData"))
	}

	

	#######################
	### QQplot
	######################
	for(ii in 1:length(exp)){
		# Per plot extract the maximum quantile values
		max.all = NULL
		for(w in exp[[ii]]){
			load(paste0("Results_QQplot_Manhattan/Pval_and_q_val_", w, ".RData"))  		# vector of the maximum quantile value per experiment
			max.all = c(max.all, max(q_val))	
		}

		# Site and trait in the title
		plot.title = unlist(strsplit(w, ".", fixed=TRUE))
 		plot.title = c(plot.title[1], plot.title[4])	

		# Save QQplot
		pdf(paste0("Results_QQplot_Manhattan/QQplot_", paste(plot.title[1], plot.title[2], sep="_"),".pdf"))
		# "ref" is the Null Hypothese 
		ref=cbind(seq(0,max(max.all, max1.all[ii]),0.01),seq(0,max(max.all, max1.all[ii]),0.01))
		
		plot(ref,cex=0.1,type="l",main=plot.title ,xlab="-log10(expected p-value)",ylab="-log10(observed p-value)",xlim=c(0,max(max.all, max1.all[ii])))

		load(paste0("Test_Result_julie/Results_GWAS_No_Corrections/q_val1_", colnames(list.pval1[[ii]]), ".RData"))
		points(q_val1, cex=0.5, col="darkgreen")

		
		for(w in exp[[ii]]){
			
			load(paste0("Results_QQplot_Manhattan/Pval_and_q_val_", w, ".RData")) # one of the R-object in .RData is "q_val"
			points(q_val, cex=0.5, col=rainbow(length(exp[[ii]]))[which(exp[[ii]]==w)])
		}

		abline(h=-log10(alpha/nrow(q_val)), lty=2,col="red")
		#abline(h=gao, lty=2,col="orange")
		legend("bottomright",legend=c("diagonal (H0)", paste("No corrections", colnames(list.pval1[[ii]]), sep=" "), exp[[ii]]), fill=c("black", "darkgreen", rainbow(length(exp[[ii]]))),cex=1)

		dev.off()
	}
}



#################################################
## Manhattan plots ##
########################################
manathan.plot=function(snp, pheno, alpha, cex=0.3){
	snp = as.data.frame(snp)
	snp = snp[order(snp$Chromosome, snp$Position),]
	pheno = as.data.frame(pheno)

	for(j in 2:ncol(pheno)){
		res=read.table(paste0("FASTLMM/Results_FASTLMM/RESULTS_FASTLMM_", names(pheno)[j], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
		#res=read.table(file=paste0("FASTLMM_COMBINES_AvecLigneeIllegitimes_AvecPriseEnCptLigneesPuresDans600K_UniquePosition/RESULTS_FASTLMM_COMBINES_", names(pheno)[j], "_AllsansChrom.csv"), sep=",", header=TRUE, stringsAsFactors=FALSE)
		# Put order (chromosome and position) the markers in "res"
		res=res[order(res$Chromosome, res$Position, decreasing=FALSE),]

		# Extract the markers of the "snp" object that are present in "res"
		snp.res = snp[snp$SNP.name %in% res$SNP, ]
		res = res[res$SNP %in% snp$SNP.name,]

		pos.sp=split(snp.res$Position,snp.res$Chromosome)
		assoc.sp=split(res$Pvalue,snp.res$Chromosome)
		chrs=unique(snp.res$Chromosome)
		nb=length(chrs)
		impairs=seq(1,nb,by=2)
		SeuilBonf=-log10(alpha/length(snp.res$Position))
		maxi=max(c(-log10(res$Pvalue),SeuilBonf))
		# Build and save the plot
		pdf(paste0("Results_QQplot_Manhattan/Manhattan_", names(pheno)[j],".pdf"))
		par(mfrow=c(1,(nb+1)),mar=c(2,4,2,0))
		plot(pos.sp[[1]],-log10(assoc.sp[[1]]),cex=cex,xlab="",ylab="",ylim=c(0,maxi),las=1,axes=FALSE,type="n")
		axis(2,0:maxi,0:maxi,las=1)
		mtext("-log10(Pvalue)",2,line=2)
		par(mar=c(2,0,2,0))		
		for (i in 1:nb){
			pos.tmp=pos.sp[[i]]
    			if (i%in%impairs){
				plot(pos.tmp,-log10(assoc.sp[[i]]),cex=cex,xlab="",ylab="",ylim=c(0,maxi),las=1,axes=FALSE)
    				mtext(paste("chr",i),1,line=0)} else {
				plot(pos.tmp,-log10(assoc.sp[[i]]),col="red",cex=cex,xlab="",ylab="",ylim=c(0,maxi),las=1,axes=FALSE)
    				mtext(paste("chr",i),1,line=0)  
    			}
    			#if (i==8) { # VGT1 et VGT2
				#abline(v=131984851,col="green")
				#abline(v=123506141,col="red")
			#} 
			#abline(h=SeuilBonf,col="red")
		} 
		dev.off()
		rm(snp.res, pos.sp, assoc.sp)
	}
}

