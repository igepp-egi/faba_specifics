# README

Je t'envoie pas mal de choses :

- PEAMUST_vfaba_F1b.vcf.gz contient les SNPs et les INDELS avec - de 15% d'Htz, une MAF à 1% et - de 10% de missing data
- PEAMUST_vfaba_F2a.vcf.gz est la matrice que j'utilise pour les calculs de structure / PCA... C'est F1b mais avec un filtre sur le LD réalisé avec plink. On enlève par fenêtre glissante de 50 SNPs tous les SNPs avec plus de 0.2 de R².
- PCA_new.Rdata est la PCA réalisé par le package "adegenet". Celà permet ensuite de calculer la DAPC qui a été utilisé pour déterminer les groupes. On l'utilise aussi comme covariable GWAS
- PEAMUST_vfaba_F1_SNP_imputed_MAF5_kinship.cXX.txt est une matrice de kinship calculé par GEMMA
- PEAMUST_vfaba_F1_SNP_imputed_MAF5.vcf.gz est la matrice pour la GWAS. C'est F1b uniquement SNP mais avec une imputation par beagle puis une MAF à 5%.

Pour ce qui est des groupes génétiques, il ont été réalisé par DAPC (package adegenet) qui réalise un clustering sur PCA pour déterminer les groupes. Le nombre de groupe a été validé après observation des données sur une carte géographique. L'arbre n'a été produit que dans un second temps car nous avons utilisé une méthode imaginé au labo avec un collègue phylogénéticien. Nous avons réalisé des figures mélangeant structure et phylogénie, je t'en ai mis une en pièces jointes.