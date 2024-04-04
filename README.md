# Faba bean resistance to aphids

## Summary

Within the SPECIFICS project (PIA3 Programme de Recherche Prioritaire "Cultiver et Protéger Autrement" SPECIFICS), a collection of 259 genotypes of _Vicia faba_ (faba bean) were phenotyped for their resistance to _Acyrthosiphon pisum_ and _Aphis fabae_.

Resistance is estimated by measuring the offspring number obtained after 15 days at 18°C from one single L1 24h larvae of A. pisum or A. fabae.   

These results will be used to conduct a GWAS analysis on the collection of the 259 exome-sequenced faba bean genotypes.

The list of faba genotypes can be found in "01.data/03.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv".

## Results A. pisum 

### Analysis 1 ("quantitative robust"): N>=6 and no zeros

Criteria:
- Only genotypes with more than 6 replicates are kept.  
- When total number of larvae is equal to 0, the row is removed and value is discarded. 

<img src="02.analyses/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot03_apisum_most_contrasted_genotypes.png" alt="number of larvae per genotype" width="80%"> 

### Analysis 02 ("qualitative zeros"): N>=6 and proportions of zeros

Criteria/method:
- Only genotypes with more than 6 replicates are kept.  
- The proportion of zeros from the total number of measurements is calculated.

<img src="02.analyses/01.acyrtosiphon_pisum_on_faba_collection/02_analysis2_proportion_zeros/plot02_number_of_true_zeros.png" alt="number of true 0 per genotype" width="80%"> 


