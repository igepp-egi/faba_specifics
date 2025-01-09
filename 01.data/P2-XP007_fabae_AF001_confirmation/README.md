# P2-XP007: A. fabae AF001 fecundity test on Vicia faba (10/10 genotypes) 

## Objective 

The 10 most sensitive and 10 most resistant genotypes of Vicia faba identified in the 2023-2024 screen (see `01.data/01.aphis_fabae...` and `02.analyses/01.aphis_fabae.../01_analysis1`) were chosen to be confirmed within the same experiment. 

Indeed, a considerable batch effect was observable in the big screen performed on the 251 genotypes in 2023-2024.  

Once the most contrasted genotypes are confirmed, they can be further investigated in follow-up experiments. 

## Protocol

### L1 larvae aphid rearing

P123 is a clone of _Acyrtosiphon pisum_ from the pea biotype (cf Ollivier et al. 2022 TAG).  

- Step 1: G0 adults were transferred on Castell faba bean plants for 24h to obtain several single one-day old L1 larvae (G1-L1).   
- Step 2: these G1-L1 larvae developped into adults (G1-adults) during approx. 15 days at 18-/+ 1.5 degrees.   
- Step 3: these G1-adults were transfered for 24h on new faba Castell plants to obtain G2-L1 larvae.
- Step 4: one L1 larvae was placed on each plant for 15 days.   

## Plants 

The 10 most sensitive and 10 most resistant genotypes of _Vicia faba_ identified in the 2024 screen (see `01.data/01.aphis_fabae...` and `02.analyses/01.aphis_fabae.../01_analysis1`) were chosen to be confirmed within the same experiment. 

|ivip       |phenotype          | median_larvae|
|:----------|:------------------|-------------:|
|IVIP102198 |Ten most resistant |          22.0|
|IVIP102302 |Ten most resistant |          22.0|
|IVIP102737 |Ten most resistant |          22.5|
|IVIP101659 |Ten most resistant |          24.0|
|IVIP102547 |Ten most resistant |          24.0|
|IVIP102361 |Ten most resistant |          25.0|
|IVIP102516 |Ten most resistant |          25.0|
|IVIP100947 |Ten most resistant |          26.0|
|IVIP101293 |Ten most resistant |          26.0|
|IVIP101482 |Ten most resistant |          27.0|
|IVIP103268 |Ten most sensitive |          46.0|
|IVIP101240 |Ten most sensitive |          46.5|
|IVIP101059 |Ten most sensitive |          47.0|
|IVIP101738 |Ten most sensitive |          47.0|
|IVIP101799 |Ten most sensitive |          47.0|
|IVIP100355 |Ten most sensitive |          47.5|
|IVIP100285 |Ten most sensitive |          49.0|
|IVIP100887 |Ten most sensitive |          49.0|
|IVIP100376 |Ten most sensitive |          50.0|
|IVIP100336 |Ten most sensitive |          51.0|

## Growth chamber used 

Growth chamber used: number 132 (after installation up to counting).

## Data dictionary

- genotype: short code used for practical ease of use in experiments (1 to 251).     
- seed_batch: the seed lot used in the experiment (from own greenhoue multiplication).   
- ivip: seed stock number accession identifier.   
- aphid: aphid species and clone used in experiment.   
- phenotype: resistance class as infered from fecundity test results from the main 251 screen.
	- PR: partially resistant  
	- S: sensitive  
- rep: plant individual replicate number (biological replicate).  
- age_plant: age of the faba bean plant in days.  
- start_date: start date of the experiment (installation date).
- end_date: end date of the experiment (total offspring count) 
- adult_alive: number of adults alive seen on the plant.  
- adult_winged: number of winged adult seen on the plant. 
- offspring_count: total number of larvae (all apterous stage) counted visually from one single G2-L1 larvae.  
- comment: optional free-text comments. 


## Reference of the protocol

Ollivier, R., Glory, I., Cloteau, R. et al. A major-effect genetic locus, ApRVII, controlling resistance against both adapted and non-adapted aphid biotypes in pea. Theor Appl Genet 135, 1511–1528 (2022). https://doi.org/10.1007/s00122-022-04050-x