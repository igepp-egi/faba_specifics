# README

Measurements in the Near Infrared (NIR) using FT-NIR (Fourier-Transform Near Infrared Imaging).  
Explained [here (in French)](https://www.bruker.com/fr/products-and-solutions/infrared-and-raman/ft-nir-spectrometers/What_is_FT-NIR_spectroscopy.html)   


# Instruments 

## Bruker MTA FT-NIR device 

Burker MTA 
Closest reference is for the [MTA II](https://www.bruker.com/en/products-and-solutions/infrared-and-raman/ft-nir-spectrometers/mpa-ii-multi-purpose-analyzer.html)

### Solid samples  
There are two types of recipients:   
- For small quantities, there are "HPLC vials" made of glass that can be cleaned with compressed air (wear safety goggles!) and rinsed with alcohol, then left to dry.  
- For larger quantities, such as seeds, there is a rotating sample cup holder and sample cups. The sample cups have a quartzite surface (which does not absorb any laser). Various sample cup diameters are available with a minimal diameter of 4cm.
- A round metal disk should be placed on top of 
- An adapter for Petri dish also exists.   

### Liquid samples

- It is also possible to use 1ml of liquid to measure absorbance.   
- The sample should not be mainly composed of water but rather be a mixture of organic compounds. For instance, two molecules diluted in water do not work  
IN312-P

## MicroNIR Pro

### Demo videos:    
- [Connect the MicroNIR](https://www.youtube.com/watch?v=vpf0aBHA4Xs&list=PLMNcKV8MjC-KjinkbT7S-5d52IIROskhe&index=13)  

- [Spectral acquisition](https://www.youtube.com/watch?v=HfZtHaUbnUc&list=PLMNcKV8MjC-KjinkbT7S-5d52IIROskhe&index=15)  

### Comments: sur la prise de spectres

- En fonction du stade de développement des plantes, noter l’étage foliaire
- Aplanir la feuille lors de la mesure
- Faire plusieurs acquisitions pour un même échantillon (ex : 3 spectres sur 3/6 endroits de la feuille)

Tester différents protocoles car influence de l’épaisseur :  
- Sans rien derrière, en appuyant juste avec la main pour aplanir la feuille, les doigts sur le plastique noir du support
- Avec le spectra long derrière : cet accessoire peut servir quand il n’y a pas assez de retour de la lumière (ie. feuilles très fines)
- Avec un support noir mat positionné derrière la feuille, ce qui permet de ne capturer que la lumière réfléchie.
 
A tester également : influence de la face (dessus/dessous) de la feuille pour le spectre ?

Réglage de l’appareil :  
- Nombre de scans :
- Par défaut, l’appareil est réglé sur 100 (le scan dure 1 à 2 secondes).

Gamme spectrale : 950 à 1650 nm

Possibilité d’augmenter à 400 scans et de déplacer la sonde sur la feuille le temps de l’analyse. Toutefois, si nos feuilles sont petites, ce n’est pas à privilégier.

Mme Taghzouti précise qu’aucun client ne lui a fait de retour quant à d’éventuelles brûlures des feuilles après plusieurs acquisitions au même endroit (à confirmer sur notre matériel).

 

Sur la nomenclature des data

Faire en sorte que le nom des spectres contienne le type de protocole d’acquisition (sous forme de numéro par exemple) pour pouvoir comparer les spectres et les résultats.
# New XP measure: measure directly from a leaf section from faba bean or pea leaf  infested or not with aphids ? 

## Reflectance from a leaf disk
Ideas to measure reflectance from a leaf round section:  
- Use a 4cm sharp steel round disk (used for cookie pastry for instance) to obtain a leaf circle of _ca_ 4cm diameter.   
- 4cm corresponds to the smallest cup size equipped with quartzite.  
- There is a need for a 4cm round metal disk covered with inox in order to reflect the laser beam and record absorbance.   

It is possible to extract a 2D-matrix "reflectance measurements x wavelengths (cm-1)" in .csv.


