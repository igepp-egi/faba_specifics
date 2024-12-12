############################################################
####### Principal Coordinates Analysis #########
############################################################
# Author: Sandra Negro
############################################################
# Date: 7 October 2015
## AIM: Multidimensional scaling of the kinship matrix (analysis of the similarity between inbred lines) and visualization of the major genetic groups on the Graph.

##  "stats" package should be installed and charged in your R working space

## INPUT:
# (1) "kin" = the kinship matrix (symetrical)
# (2) "adm.name" = the Admixture output name with the path and without the extension (ex.: adm.name = "./Results_Admixture/Matrix_247_4779_2015-06-26")
# (3) "nbQ" = Nb of genetic groups you think you have
# (4) "user.class" = a matrix containing two columns with a header: the first column contains the individuals names and the second column contains the user-classification, which can be characters or numbers.
# (5) "user.classification" = TRUE or FALSE (to specify if the user has a "user.class" file

## OUTPUT:
# (1) Creation of a "Results_PCoA" folder
# (2) PCoA graphs
# (3) "PCoA-Points.csv" = matrix with the individuals in rows and the axes (10 dimension here) in columns
# (4) "ProportionVariance.txt" = proportions of the variance explained by the axes

## REFERENCE: 
# Gower (1967) Multivariate Analysis and Multidimensional Geometry. Journal of the Royal Statistical Society. Series D (The Statistician), Vol. 17, No. 1 (1967), pp. 13-28 

############################################################

PCoA <- function(kin, adm.name, nbQ, user.classification){
	dir.create("./Results_PCoA")
	### PCoA on the distance matrix
	Dist=1-kin
	PCOA=cmdscale(Dist, 10, eig=TRUE) # 10 dimensions of the space
	#str(PCOA)
	# Proportion of the variance explained by the axes
	propvar=PCOA$eig/sum(PCOA$eig)
	
	### Save the proportions of the variance explained by the axes
	write.table(PCOA$points , file="./Results_PCoA/PCoA-Points.csv", sep=",", row.names=TRUE, col.names=FALSE, quote=FALSE)
	propaxes = cbind(paste0("axe", seq(1,length(propvar), by=1)) ,propvar)
	write.table(propaxes, file="./Results_PCoA/ProportionVariance.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	
	### Table with a group allocated to each inbred line
	Color=matrix(NA,nbQ,nrow(kin))
	colnames(Color) = rownames(kin)
	col=list()
	##col.i=NULL
	col.i=c("red","blue", "green", "black", "purple", "orange", "darkgreen", "turquoise", "gray", "darkblue", "maroon", "pink", "tomato", "yellowgreen")

	for(i in 2:nbQ){
		##col.i = rainbow(i)
		colouration = col.i[1:i]
		col[[i]] =colouration	

		admixture=read.table(file=paste(adm.name, i,"Q", sep="."),sep="\t", header=FALSE, row.names=1, stringsAsFactors=FALSE)

		Color[i,]=col[[i]][apply(admixture,1,function(a) which(a==max(a)))]		
	}
	
	### Graph
	#jpeg(paste0("./Results_PCoA/PCoA_",nbQ,"Q.jpg"))
	pdf(paste0("./Results_PCoA/PCoA_",nbQ,"Q.pdf"))
	plot(PCOA$points[,1:2],col="white",main="PCoA on the distance matrix",xlab=paste("1st axis:",signif(propvar[1],2)*100,"% of the variance"),ylab=paste("2nd axis:",signif(propvar[2],2)*100,"% of the variance",sep=" "))
	text(PCOA$points,labels=rownames(kin),cex=0.7,col=Color[nbQ,])
	#legend(title="Groups:","topleft",col=col[[nbQ]], legend=col[[nbQ]],lty=1)
	dev.off()

	
	
	### PCoA Graph with colouration based on the user classification
	if(user.classification == TRUE){
		# Check that the individuals in the kinship and the user classification file are identical and in the same order
		user.class = user.class[user.class[,1] %in% rownames(kin), ]
		user.class = user.class[match(rownames(kin), user.class[,1]), ]
		if(dim(kin)[1] != dim(user.class)[1]){print("ATTENTION: the Nb of individuals are different in the kinship and the user.class file")}
	
		user.col = col.i[1:length(unique(user.class[,2]))]
		user.class = as.data.frame(user.class)
		user.class[,2] = as.factor(user.class[,2])
		levels(user.class[,2]) = user.col
	
		pdf("./Results_PCoA/PCoA_User_Classification.pdf")
		plot(PCOA$points[,1:2],col="white",main="PCoA on the distance matrix",xlab=paste("1st axis:",signif(propvar[1],2)*100,"% of the variance"),ylab=paste("2nd axis:",signif(propvar[2],2)*100,"% of the variance",sep=" "))
		text(PCOA$points,labels=rownames(kin),cex=0.7,col=as.character(user.class[,2]))
		#legend(title="Groups:","topleft",col=col[[nbQ]], legend=col[[nbQ]],lty=1)
		dev.off()
	}
	
}





