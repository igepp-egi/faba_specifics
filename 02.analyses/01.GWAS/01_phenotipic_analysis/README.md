# README

This folder contains the different analyses performed on the ['big screen' of the 251 genotypes](../../01.data/2023_screen_specifics_251genotypes/):  

1) The first analysis called '01.dataset_filtering' and has two subfolders for each aphid species:  
    - 'Testing all filters' test the effect of several filters (IQR, manual filtering and no filter) on the normality of the distribution. 
    - 'final filters' contains the filters used to generate the 3 datasets (N >= 4, 5 or 6 biological replicates).  
2) The second analysis called '02.ajusted_means' has two subfolders, one for each aphid species and contains the extracted adjusted means (emmeans) from the linear models fitted. These values will be used for GWAs as they correspond to the BLUEs (Best Linear Unbiaised Estimators) for each genotype.   
3) The third analysis called '03.n_larvae_per_genetic_cluster' shows the number of offspring from one L1 on the different genetic groups.




