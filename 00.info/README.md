# Passeport and germplasm information

This directory list the IVIP codes, the geographic origin, names of landraces/cultivars and the genetic cluster of the different lines. 

# Data dictionary

- **ivip:** INRAE Dijon "Centre de Ressources Génétique" code.  
- **historical_id:** INRAE Dijon former code.   
- **species:** Vicia faba species and subspecies.   
- **full_name:** code, cultivar name or traditional name. 
- **status:** landrace or improved cultivar. When unambiguous status is present, both are conserved.    	
- **geo_origin:** country of origin when known. NA: not available. 
- **cluster_dapc:** genetic cluster as determined by the Discriminant Analysis of Principal Components (DAPC). 
- **cluster_name:** cluster name through mapping of the geographic accessions on a world map.   
- **winterH_springP:** winter ("hiver", H) or spring ("printemps", P) landrace / cultivar. 


# References

## DAPC method (from adegenet package)  
Grünwald NJ., Goss EM. 2011. Evolution and population genetics of exotic and re-emerging pathogens: Novel tools and approaches. Annual Review of Phytopathology 49:249–267. Available at: http://www.annualreviews.org/doi/abs/10.1146/annurev-phyto-072910-095246?journalCode=phyto  

Jombart T., Devillard S., Balloux F. 2010. Discriminant analysis of principal components: A new method for the analysis of genetically structured populations. BMC genetics 11:94. Available at: http://www.biomedcentral.com/1471-2156/11/94  

Pritchard JK., Stephens M., Donnelly P. 2000. Inference of population structure using multilocus genotype data. Genetics 155:945–959. Available at: http://www.genetics.org/content/155/2/945.abstract  
