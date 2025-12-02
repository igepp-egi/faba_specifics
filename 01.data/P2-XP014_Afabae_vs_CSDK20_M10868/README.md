# P2-XP014 CSDK20 and ILB540 vs A. fabae AF001

## Summary and objectives 

Goal: to find differentially expressed genes and differentially accumulated metabolites from _Vicia faba_ plantlets (aged 15 days) infested or not with 100 L3 _Aphis fabae_ larvae. 

mRNA-seq and untargeted metabolomics to be performed on the same samples. 

## Material 

biological material:
- _Aphis fabae_ AF001 (aphid clone sampled in Le Rheu Spring 2024).
- _Vicia faba_: 
	- CS20DK (IVIP102198): identified as the most partially resistant to AF001 in previous screenings.   
	- M10868 (IVIP100355) identified as the most sensitive in previous screenings.  

material:

Plant were grew in small plastic pots (6,5cm * 6,5cm * 8cm), rigide perfored plastic bag cut from both sides were placed around the pot once the plant had germinated and was maintained using 2 tutors (see picture in the experimental_design file)

## Complete randomized block design 

See [./XP_design_P2-XP014.pdf](XP design file)

## Temperatures and hygrometry

18°C D:N and long days 16:8 D:N. 

In chamber 133 (from seed to 15-day old seedlings): probes #38 to 43 (September 8 to September 22, 2025).  
In chamber 131 (after aphid installation), probes #33 to 37 (September 22 to September 27, 2025).  

see temperature_hygrometry_P2-XP009_and_P2-XP010 files for more informations on temperature and hygrometry

## protocol

presentation of the parameters:
2 plant genotypes (cf plant material) * 2 conditions (aphid treatment and control treatment) * 4 time points (3hpi, 24hpi, 72hpi, 120hpi) * 10 repetitions --> 160 datapoints

-plant were planted 14 days before the experiment
-8 days later the plastic bags were placed around the plants
-Before dropping aphids all control plastic bag were closed using 3-4 staples
-For the aphid treatment 100 3 days old larvae were dropped on each plants. After the dropping the plastic bags were closed using staples. 
-At 4 time point (counted in hour post infestation (hpi)) the plants were frozen in liquid nitrogen. Aphids were removed with a brush and the frozen plant were stored in folded aluminium foil. During the sample harvest the aluminium foils were kept in a nitrogen refrigerated polystyrene box before beeing transfered in -80 freezer to be stored.
- Plants extracts were then grinded at a later date in mortar and pesle with liquid nitrogen. During the process the powder was separated into 3 tubes 2 2ml tubes (1 for transcriptomic and 1 for metabolomic) and a backup tube of 15ml. 
-metabolomic tubes were hyophilized during 48h. Pierced 2ml tubes cap were used to allow water evacuation without powder contamination.
-From the lyophilized dried powder about 10mg was placed in a numeroted well plate for each datapoint. The exact weight is indicated in the samples_to_conditions file.

## Comments/observations:
At the last time point we had contamination on some control plants. We found few A. fabae larvae on the plants these are noted on the samples_to_conditions file.


## mRNA-seq protocol

Macherey-Nagel NucleoSpin RNA protocol. 



## Metabolomics protocol

**We came up with the following:**  
1. We will make one extract and use it for the phytohormones and the metabolomics.
2. Please use 10 mg (dry weight) plant material. Please make an excel table with the exact weights. Just in case you have less, just use less, we will try to adapt the procedure to it (but please indicate).
3. I will send you 96 well plates. Please put weight the plant material in the vials of the 96 well plate and close the lid tightly (at the moment the lid is loose on top of the wells).
4. It will take 2 to 3 days to prepare the samples for one experiment. The phytohormone measurements will take about 2 days. The metabolomic measurement will be done more or less at the same time. One mode takes about 2.5 days. If you want to have both modes (positive and negative ionisation mode) it will take a week.
5. And for the last experiment it might be saver to do it end of January. (just in case there is some delay in one of the steps, we come into the Christmas time and we come into trouble, if we do it in December)


**What we need to decide:**
Is it enough for the metabolome to measure just in one mode (negative or positive) or do you want both? One mode would be enough to tell whether they are different. But not all compounds ionize in both modes. So if you are interested to get "all" compounds, which you can detect with that method, we have to run the samples in both modes.

And since there are a lot of samples to analyze it would be good if you could buy two columns for the HPLC which we need to do the analyses. Then we can start with new, clean columns.
The following columns are needed for the planned analyses:

1) phytohormones:
Name: Kinetex 2.6um Biphenyl 100 A
Dimensions: 100 x 4.6mm
Manufacturer: Phenomenex
Product/order-Number: 00D-4622-E0

2) Metabolomics:
Name: Zorbax Eclipse XDB-C18
Dimensions: 21. x 100mm, 1.8um particles
Manufacturer: Agilent
Product/order-Number: 981758-902
 
