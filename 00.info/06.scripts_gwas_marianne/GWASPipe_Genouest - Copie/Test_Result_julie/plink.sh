
#! /bin/bash 


eval "$(conda shell.bash hook)"_
conda activate /usr/local/env/plink_1.90b6 



	

plink --noweb --file /home/jferreir/GWAS_Pipe_S_Nicolas/Test_Result_julie/Results_plink/MatrixGeno_Chr10 --compound-genotypes --make-bed --out /home/jferreir/GWAS_Pipe_S_Nicolas/Test_Result_julie/Results_plink/MatrixGeno_Chr10
