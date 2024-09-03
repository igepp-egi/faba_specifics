Tuto GWASpipe

1. Installation de la machine virtuelle

Se connecter sur GenOuest

Genostack = outil pour gérer les machines virtuelles
Pour accéder à genostack: via le lien ou sur genouest: tools>cloud>using genouest cloud.
Lien genostack pour créer et suivre la machine virtuelle: https://genostack.genouest.org/auth/login/?next=/

Dans Genostack Compute > Instances > launch instances 
Remplir

    names: gwas_pipe (nom qu'on veut donner) (create new volume "NO")

    source: rechercher image de Fabrice "IGEPP_GQ" > cliquer sur petite fleche

    flavor (taille): "m2.2xlarge" pour colza 30k SNP mais on peut choisir dans les 3 dernières en fonction de nos données > flèche pour valider

    lauch instance

On a créer la machine virtuelle mais on doit ensuite faire la config pour pouvoir l'utiliser sous MobaXterm

Dans MobaXterm

Si premiere connexion :

    creer une seession qui connecte Mobaxterm à genouest

    mettre les clés ssh (suivre tuto genouest)


Si déjà configuré

    se connecter

    ssh fefora@<adresseIPmachineVirtuelle>

    sudo yum install nano (installer nano - editeur de texte) > y

    modifier un texte pour créer les utilisateurs mais attention un seul utilisateur autorisé à la fois

    sudo nano /etc/ssh/sshd_config

    Changer le paramètre suivant :

    To disable tunneled clear text passwords, change to no here!

    PasswordAuthentification yes > no

    cliquer crtl + x > y

    puis enter pour enregistrer ficher et en sortir

    redemarrer sshd - créer un port pour acceder à RStudio en ligne

    sudo service sshd restart

    créer nos identifiants et mots de passe

    sudo useradd <memeNomQueMachineVirtuelle> -m -g igepp

    créer mot de passe

    sudo passwd <memeNomQueMachineVirtuelle>

    enter nouveau mot de passe

    exit


Aller sur internet (enlever les points dans l'ip et remplacer par des tirets)
https://myapp-[adresseIpMachineVirtuelle].vm.openstack.genouest.org/
Remarque version de R 1.2.50.33 => Il faudrait mettre à jour la version R car tres vielle et ca peut poser probleme si on veut faire du développement  downstream pour plot etc par exemple en sortie de FaST-LMM mais verifier que FaST-LMM tourne sur une version plus récente de R

rentrer identifiant <memeNomQueMachineVirtuelle> ou igepp et mot de passe rentré précedemment
Et on est maintenant dans RStudio en ligne dans la machine virtuelle
Télécharger le script dans RStudio

Les différentes étapes du pipeline :

    filtres des SNP sur les données génotypiques, input 

    snp.dm données manquantes ici tolérées à 10%

    maf à 5% on peut baisser à 1 si pannel tres grand combien d'individus minimum dans une classe allélique

    ATTENTION - il faut prétraiter les données phénotypiques - là on a une valeur par génotype et il faut avoir traitement les données phénotypiques avant - les residus doivent etre distribués normalement (données propre de phénotypaqge) - si pas distribué normalement il faut faire une autre méthode que FaST-LMM

    on fait le QC 

    Remarque la filtre drastique sur l'hétérozygocie ce qu'il faut adapter

    filtre et imput SNP avec BEAGLE



## Pipeline step #1: clean-up and formatting of VCF raw data
    inputs = 1) VCF +  2) Phenotype
1. VCF to genotype: convert VCF to Rdata
1. SNP description = before filtering (in R)
2. Filter and impute SNPs  (heterozygoty, maf, missing data) (in R) BEAGLE
3. Check consistency of genotype and phenotype identifiers (in R).
4. Write ..ped / .bed / .bam / .fam etc. files to be used by PLINK (bash)
5. Import .ped / .bed / .bam / .fam data files in PLINK (bash)

## Pipeline step #2: get population kinship (in R)
inputs = pre-processed PLINK files (.bed etc.) 
ouputs = matrix of similarities (.csv) between individuals (all genome + per chromosome)

2 steps :
  - using GenABLE (on all chromosomes)
  - per chromosome ssCHr (different algorithms possible Astle or Vasraden)

## Pipeline step #3: get population structure

Option 1:   Admixture (bayesian method, in bash)

inputs = pre-processed PLINK files (.bed etc.)  + number of genetic groups.

Outputs = matrix of groups

1. get population structure using a PCoA before the analysis (Python or R)

2. run admixture (bash)

3. create genetic stacked bar plot based on number of groups. 

Option 2: PCoA (inferential method, in Python)

inputs = pre-processed PLINK files (.bed) + number of genetic groups. 

outputs = plot + matrix of dissimilarities .csv

## Pipeline step #4: perform GWAS analysis
inputs = genotype format PLINK + kinship .csv + structure file .csv (PCoA/Admixture) 
FASTLMM (in R)
outputs = one .csv file of SNPs (one line per SNP and per phenotype)

reformat Fast-LMM output using GNPasso
R2 is phenotypic variance explained per SNP
Remark: Pvalue is the one not corrected - should do FDR 5 to 20% -> Rapeseed  (Gao 2008 & 2010 not use Bernouilli) - would be good to look at pvalue districution if good or not before/after correction

## Pipeline step #5: Create manhattan plot
Inputs = .csv of SNPs with p-values
Output = Manhanttan plot














