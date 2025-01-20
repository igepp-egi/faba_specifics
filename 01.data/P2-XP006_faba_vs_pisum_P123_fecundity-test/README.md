# P2-XP006: A. pisum P123 fecundity test on Vicia faba (10/10 genotypes)

## Objective 

The 10 most sensitive and 10 most resistant genotypes of __Vicia faba__ identified in the 2022-2023 screen (see `01.data/01.acyrthosiphon...` and `02.analyses/01.acyrthosiphon.../01_analysis1`) were chosen for confirmation. 

These genotypes will be confirmed within the same experiment. Indeed, a considerable batch effect was observable in the big screen performed on the 251 genotypes in 2022-2023.  

Once the most contrasted genotypes are confirmed, they can be further investigated in follow-up experiments. 

## Protocol

### L1 larvae aphid rearing

P123 is a clone of _Acyrtosiphon pisum_ from the pea biotype (cf Ollivier et al. 2022 TAG).  

- Step 1: G0 adults were transferred on Castell faba bean plants for 24h to obtain several single one-day old L1 larvae (G1-L1).   
- Step 2: these G1-L1 larvae developped into adults (G1-adults) during approx. 15 days at 18-/+ 1.5 degrees.   
- Step 3: these G1-adults were transfered for 24h on new faba Castell plants to obtain G2-L1 larvae.
- Step 4: one L1 larvae was placed on each plant for 15 days.   

## Genotypes to confirm 

The 10 most sensitive and 10 most resistant genotypes of _Vicia faba_ identified in the 2024 screen (see `01.data/01.acyrthosiphon...` and `02.analyses/01.acyrthosiphon.../01_analysis1`) were chosen to be confirmed within the same experiment. 

|ivip       |phenotype          | median_larvae|
|:----------|:------------------|-------------:|
|IVIP101957 |Ten most resistant |          15.5|
|IVIP101240 |Ten most resistant |          16.0|
|IVIP101986 |Ten most resistant |          17.5|
|IVIP100204 |Ten most resistant |          18.5|
|IVIP100123 |Ten most resistant |          19.0|
|IVIP100310 |Ten most resistant |          19.5|
|IVIP100388 |Ten most resistant |          19.5|
|IVIP101550 |Ten most resistant |          19.5|
|IVIP100396 |Ten most resistant |          20.0|
|IVIP100382 |Ten most resistant |          21.0|
|IVIP100234 |Ten most sensitive |          35.0|
|IVIP100899 |Ten most sensitive |          35.0|
|IVIP101599 |Ten most sensitive |          35.0|
|IVIP102547 |Ten most sensitive |          35.0|
|IVIP100872 |Ten most sensitive |          36.0|
|IVIP102492 |Ten most sensitive |          36.0|
|IVIP101320 |Ten most sensitive |          38.0|
|IVIP101736 |Ten most sensitive |          38.0|
|IVIP102516 |Ten most sensitive |          38.0|
|IVIP103069 |Ten most sensitive |          40.0|

## Growth chamber 

Growth chamber used: L057 132 (after installation up to counting).

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
- larvae_total_number: total number of larvae (all apterous stage) counted visually from one single G2-L1 larvae.  
- comment: optional free-text comments. 

## Reference of the protocol

Ollivier, R., Glory, I., Cloteau, R. et al. A major-effect genetic locus, ApRVII, controlling resistance against both adapted and non-adapted aphid biotypes in pea. Theor Appl Genet 135, 1511–1528 (2022). https://doi.org/10.1007/s00122-022-04050-x