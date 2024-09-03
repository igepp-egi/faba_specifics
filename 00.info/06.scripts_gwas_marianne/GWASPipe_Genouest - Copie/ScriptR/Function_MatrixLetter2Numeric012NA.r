##################################################################
#### Fonction: MatrixLetter2Numeric012NA
#### Author: Stephane Nicolas
#### Date: 16/01/13
#### Description: Transformation du genotypage au format 0/1/2/NA        
#### Input: 
#### 1) Une matrice de genotypage IND x SNP code AA, TT, GG, CC, AT et NA en donnees manquantes (rownames=NomInd non redondants), colnames (noms SNPs)
#### 2) Panel = Nom du Panel
#### 3) Origin = Origine du genotypage TOP ou FW
#### Output:
#### 1)Une matrice de genotypage IND x SNP code 0 (homozygote), 1 (Heterozygote), 2 (Homozygote 2d allele), NA (donnees manquantes)
#### 2)Ecrit dans le repertoire de travail le fichier de genotypage
#### 3)Ecrit dans le repertoire de travail le fichier de correspondance
#### To do: Liste avec sortie Corr et Matrice Genotypage 
##################################################################

MatrixLetter2Numeric012NA=function(test,Panel,Origin){
Type="Matrix_Genotyping_012NA"
test2=as.matrix(test)
if(class(test[,2])=="factor"){print("matrice d'entree a des colonnes de classe facteurs, transformation au format caractere");test=apply(test,2,as.character)}
#test2[which(test2=="AG"|test2=="GA"|test2=="CG"|test2=="GC"|test2=="CT"|test2=="TC"|test2=="GT"|test2=="TG"|test2=="AT"|test2=="TA"|test2=="AC"|test2=="CA",arr.ind=TRUE)]
#Transformation des heterozygotes en 1 partir d'une matrice de genotypage
test2[which(test2=="AG"|test2=="GA"|test2=="CG"|test2=="GC"|test2=="CT"|test2=="TC"|test2=="GT"|test2=="TG"|test2=="AT"|test2=="TA"|test2=="AC"|test2=="CA",arr.ind=TRUE)]="1"
#Transformation des donnees manquantes en 4 ? partir d'une matrice de genotypage
test2[which(is.na(test2),arr.ind=TRUE)]="4"
#Stockage des correspondance Lettre, 0,2
#Attention ne marche que s'il n'y a pas de sample name redondant...
Allele0=apply(test2,2,function(x){x[x!="1"&x!="4"][1]})
Allele2=apply(test2,2,function(x){x[x!="1"&x!="4"&x!=x[x!="1"&x!="4"][1]][1]})
CorrAllele=cbind(Allele0,Allele2)
#length(CorrAllele[which(is.na(CorrAllele[,2])),])#Monomorphe
unique(test2[,which(is.na(CorrAllele[,2]))[1]])
#Transformation homozygotes en 0,2
#Attention ne marche que s'il n'y a pas de sample name redondant...
test2=apply(test2,2,function(x){x[x==x[x!="1"&x!="4"][1]]="0";x[x==x[x!="1"&x!="4"&x!="0"][1]]="2";x})
#dim(test2)
#Transformation donnees manquante 4 en NA
test2[which(test2==4,arr.ind=TRUE)]=NA

#Verification que toutes les lettres ont ete transformes
if(length(which(test2=="AG"|test2=="GA"|test2=="CG"|test2=="GC"|test2=="CT"|test2=="TC"|test2=="GT"|test2=="TG"|test2=="AT"|test2=="TA"|test2=="AC"|test2=="CA",arr.ind=TRUE))==0){print("Transformation Heterozygotes Ok")}else{print("Transformation Heterozygotes Pb")}
if(length(which(test2=="AA"|test2=="GG"|test2=="CC"|test2=="TT",arr.ind=TRUE))==0){print("Transformation Homozygote Ok")}else{print("Transformation Homozygote Pb")}

#test2[which(test2=="AA"|test2=="GG"|test2=="CC"|test2=="TT",arr.ind=TRUE)[1],]
#test2[,which(test2=="AA"|test2=="GG"|test2=="CC"|test2=="TT",arr.ind=TRUE)[2]]
#test[,which(test2=="AA"|test2=="GG"|test2=="CC"|test2=="TT",arr.ind=TRUE)[2]]

##Export fichier de genotypage 0/1/2/NA Format###
#Genotyping
#test2=cbind(rownames(test2),test2)
if(dim(test)[1]==dim(test2)[1]){print("Transformation fichier 012NA Ok")}else{print("Transformation fichier 012NA NOk")}

#Ecriture des fichiers genotypages et correspondance
#write.table(test2,paste(as.character(Type),as.character(Panel),as.character(Origin),Sys.Date(),".txt",sep="_"),sep="\t",na="NA",dec=".",row.names=TRUE,col.names=TRUE,quote=FALSE)

test2=as.matrix(test2)
row = rownames(test2)
test2 = apply(test2, 2, as.integer)
rownames(test2)=row


#save(test2,file=paste(as.character(Type),as.character(Panel),as.character(Origin),Sys.Date(),".RData",sep="_"))
save(test2,file=paste(as.character(Type),as.character(Panel),as.character(Origin),".RData",sep="_"))
write.table(CorrAllele,paste("CorrespondanceAllele_LetterATGC_012",as.character(Type),as.character(Panel),as.character(Origin),".txt",sep="_"),sep="\t",na="NA",dec=".",row.names=TRUE,col.names=TRUE,quote=FALSE)

return(test2)

##InfoLocus
#write.table(InfoLoc[,1:3],paste("SNP_info",as.character(Panel),as.character(Origin),Sys.Date(),".txt",sep="_"),sep="\t",na="NA",dec=".",row.names=FALSE,col.names=FALSE,quote=FALSE)
}


