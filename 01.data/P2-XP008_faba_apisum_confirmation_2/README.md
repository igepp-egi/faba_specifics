# P2-XP008: Second fecondity experiment _A.pisum_

## Objectives

In the first confirmation test for __A. pisum__ 4 genotypes were not available.  
Therefore, a second confirmation test is required to test the missing genotypes as well as 4 others genotypes that were already confirmed (as batch effect controls).

## Growth chamber: 

L057 chamber 132

## Temperature and humidity sensors

3 temperature and humidity sensors Milesight (EM300-TH) were placed next to the trays on Thu 23-01-2025 at 15:50. 

- Sensor 40: 1st shelf (bottom) right at the back of the room.  
- Sensor 43: 2st shelf (bottom) right at the back of the room.
- Sensor 42: 1st shelf

Scheme of the growth chamber: 

Left: left trays  
Right: right trays  

**Back of the growth chamber**

|  Left | Left  | Corridor  | Right                      | Right  |
|-------|-------|-----------|----------------------------|--------|
|       |       |           |  40: 1st shelf bottom      |        |
|       |       |           |                            |        |
|       |       |           |                      | 43: 2st shelf middle |
|       |       |           |                            |        |
|       |       |           |                            |        |
|       |       |           |                            |        |
|       |       |           |                            |        |
|       |       |           |                            |        |
|       |       |           |  42: 2st shelf bottom      |        |
|       |       |  **Door** |                            |        |


Note: sensors moved to same shelf (top shelf) during the experiment.  

Temperature and hygrometry data were recorded from installation to offspring counting (15 days). Data are available in this folder under `temperature_hygrometry_data.csv`.

T040_mean: average temperature of sensor #40.  
T043_mean: average temperature of sensor #43.  
T042_mean: average temperature of sensor #42.  

H043_mean: average hygrometry of sensor #43.  
H042_mean: average hygrometry of sensor #42.  
H040_mean: average hygrometry of sensor #40.  

### L1 larvae aphid rearing

P123 is a clone of _Acyrtosiphon pisum_ from the pea biotype (cf Ollivier et al. 2022 TAG).  

- Step 1: G0 adults were transferred on Castell faba bean plants for 24h to obtain several single one-day old L1 larvae (G1-L1).  
- Step 2: these G1-L1 larvae developped into adults (G1-adults) during approx. 15 days at 18-/+ 1.5 degrees.   
- Step 3: these G1-adults were transfered for 24h on new faba Castell plants to obtain G2-L1 larvae.  
- Step 4: one L1 larvae was placed on each plant for 15 days.   


## Genotypes 

Genotypes to test:    
- 221 / IVIP102492, South America  
- 148 / IVIP101701,  East  
- 61 / IVIP100396,  South  
- 38 / IVIP100310  North  

Plants to check experimental batch effects with P2-XP006:  
- 7 / IVIP100123 / North (Cluster 4)  
- 57 / IVIP100388 / South (Cluster 3)  
- 232 / IVIP103069 / North (Cluster 4)  
- 113 / IVIP101320 / South (Cluster 3)  

The decision to not retest genotypes #221 and #148 was taken due to the fact that they belong to genetic groups for which there are no PR in their groups = no partially resistant in South America and East groups. 

**Final list of genotypes to test:**   

|  short_code |  seed_batch                         | ivip        | cluster_name  | cluster_number  |
|-------------|-------------------------------------|-------------|---------------|-----------------|
| 61          |  P2F0824 (N=15) and JF2129 (N=5) 	| IVIP100396  |  South        | 3               |
| 38   		  |  P2F0821	                       	| IVIP100310  |  North        | 4               |
| 7  		  |  PF0755    							| IVIP100123  |  North        | 4               |
| 57  		  |  KF0170   							| IVIP100388  |  South        | 3               |
| 232 		  |  PF0783     						| IVIP103069  |  North        | 4               |
| 113 		  |  PF0765     						| IVIP101320  |  South        | 3               |


**Remarks:**   
- the reference genotype Hedin/2 used to obtain the first genome assembly of faba bean is also tested (seed batch KF2229) for its its resistance/sensitivity to A. pisum P123.   
- genotype #61 came from two batches: P2F0824 (N=15 seeds) and JF2129 (N=5).   


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
