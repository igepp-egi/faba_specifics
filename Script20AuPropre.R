###############################
####Importation des données####
###############################
pisum_zero = read.table("pisum_zéro.txt", header=TRUE, stringsAsFactors = TRUE)
pisum_zero$Code <- factor( pisum_zero$Code, levels=unique( pisum_zero$Code))
pisum_zero$Genotype <- factor( pisum_zero$Genotype, levels=unique( pisum_zero$Genotype))
pisum_zero$Nb_larves <- as.integer (pisum_zero$Nb_larves)
pisum_zero$Ad_vivant <- as.numeric (pisum_zero$Ad_vivant)
pisum_zero$Ad_ailé <- as.numeric (pisum_zero$Ad_ailé)
pisum_zero$Nb_noeuds <- as.numeric(pisum_zero$Nb_noeuds)
summary(pisum_zero)

fabae_zero = read.table("fabae_zéro.txt", header=TRUE, stringsAsFactors = TRUE)
fabae_zero$Code <- factor( fabae_zero$Code, levels=unique( fabae_zero$Code))
fabae_zero$Genotype <- factor( fabae_zero$Genotype, levels=unique( fabae_zero$Genotype))
fabae_zero$Nb_larves <- as.integer (fabae_zero$Nb_larves)
fabae_zero$Ad_vivant <- as.numeric (fabae_zero$Ad_vivant)
fabae_zero$Ad_ailé <- as.numeric (fabae_zero$Ad_ailé)
fabae_zero$Nb_noeuds <- as.numeric(fabae_zero$Nb_noeuds)
summary(fabae_zero)



##############
##Graph####
##############
#Pour pisum####
##sans les zéros 
subset10 <- pisum_zero[pisum_zero$Genotype %in% levels(pisum_zero$Genotype)[1:6], ]
subset11 <- pisum_zero[pisum_zero$Genotype %in% levels(pisum_zero$Genotype)[7:11], ]

boxplot(Nb_larves ~ Genotype, data = subset10,
        xlab = "Génotypes",
        ylab = "Nombre de larves",
        main = "Nombre de larves pour les génotypes (Pisum)",
        col = "lightgreen",  
        notch = FALSE, levels = c(7, 38, 61, 172, 14, 55, 57, 134, 170, 77, 113, 138, 148, 153, 221, 232, 241, 84), ylim = c(0,14))
boxplot(Nb_larves ~ Genotype, data = subset11,
        add = TRUE, 
        col = "lightblue")  

#Pour Fabae####
##sans les zéros 
subset7 <- fabae_zero[fabae_zero$Genotype %in% levels(fabae_zero$Genotype)[1:5], ]
subset8 <- fabae_zero[fabae_zero$Genotype %in% levels(fabae_zero$Genotype)[6:10], ]

boxplot(Nb_larves ~ Genotype, data = subset7,
        xlab = "Génotypes",
        ylab = "Nombre de larves",
        main = "Nombre de larves pour les génotypes (Fabae)",
        col = "lightgreen",  
        notch = FALSE, levels = c(198,208,189,230,45,49,52,155), ylim = c(0,50))
boxplot(Nb_larves ~ Genotype, data = subset8,
        add = TRUE, 
        col = "lightblue")  



###Analyses stats####
#####GLM.nb#### 
#Pour pisum sans zéro####
mod <- glm.nb(Nb_larves ~ Genotype,data = pisum_zero)
summary(mod)

vecteur <- residuals(mod)
vecteur
hist(vecteur)
qqPlot(vecteur)
shapiro.test(mod$residuals)
bartlett.test(Nb_larves ~ Genotype,data = fabae_zero)
plot(mod)
###conditions vérifiées 

library(emmeans)
# Calculer les moyennes estimées des effets marginaux (EMM)
emm_results <- emmeans(mod, "Genotype", na.rm=TRUE)
emm_results

plot(emm_results)

library(multcomp)
library(contrast)
pairs(emm_results)
##pas de différences significatives entre génotypes 

#Pour fabae sans zéro####
mod1 <- glm.nb(Nb_larves ~ Genotype,data = fabae_zero)
summary(mod1)

vecteur1 <- residuals(mod1)
vecteur1
hist(vecteur1)
qqPlot(vecteur1)
shapiro.test(mod1$residuals)
bartlett.test(Nb_larves ~ Genotype,data = fabae_zero)
plot(mod1)
###conditions vérifiées 

library(emmeans)
# Calculer les moyennes estimées des effets marginaux (EMM)
emm_results1 <- emmeans(mod1, "Genotype", na.rm=TRUE)
emm_results1

plot(emm_results1)

library(multcomp)
library(contrast)
pairs(emm_results1)
##pas de différences significatives entre génotypes 


##Modèle GLMTMB avec effet aléatoire rep#### 
library(glmmTMB)
library(MASS)
library(car)
library(multcomp)

mod2 <- glmmTMB(data=fabae_zero, Nb_larves ~ Genotype + (1|Genotype/Rep),family=nbinom1(link="log")) 
summary(mod2) # différence entre certains génotype 
hypotheses1 <- glht(mod2, linfct = mcp(Genotype = "Tukey"))
summary(hypotheses1)

vecteur2 <- residuals(mod2)
vecteur2
hist(vecteur2)

mod3 <- glmmTMB(data=pisum_zero, Nb_larves ~ Genotype + (1|Genotype/Rep),family=nbinom1(link="log")) 
summary(mod3) # différence entre certains génotype 
hypotheses2 <- glht(mod3, linfct = mcp(Genotype = "Tukey"))
summary(hypotheses2)

vecteur3 <- residuals(mod3)
vecteur3
hist(vecteur3)



###################################
###Génotype communs :106 et 224####
###################################

####Graph : Représenter les génotypes communs#####
##Genotype 106####
gen106_pisum <- pisum_zero[pisum_zero$Genotype == 106, ]
gen106_fabae <- fabae_zero[fabae_zero$Genotype == 106, ]

library(ggplot2)
# Tracer le graphique
ggplot() +
  geom_boxplot(data = gen106_pisum, aes(x = "Pisum", y = Nb_larves), fill = "white") +
  geom_boxplot(data = gen106_fabae, aes(x = "Fabae", y = Nb_larves), fill = "white") + 
  geom_jitter(data = gen106_pisum, aes(x = "Pisum", y = Nb_larves), 
              position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_jitter(data = gen106_fabae, aes(x = "Fabae", y = Nb_larves), 
              position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(x = "", y = "Nombre de larves", title = "Boxplot pour le génotype 106") +
  theme_minimal() + 
  scale_x_discrete(labels = c("Fabae", "Pisum")) + 
  ylim(0, 48) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, margin = margin(t = 10))) + 
  theme(axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 10)))

##Genotype 224####
gen224_pisum <- pisum_zero[pisum_zero$Genotype == 224, ]
gen224_fabae <- fabae_zero[fabae_zero$Genotype == 224, ]

ggplot() +
  geom_boxplot(data = gen224_pisum, aes(x = "Pisum", y = Nb_larves), fill = "white") +
  geom_boxplot(data = gen224_fabae, aes(x = "Fabae", y = Nb_larves), fill = "white") + 
  geom_jitter(data = gen224_pisum, aes(x = "Pisum", y = Nb_larves), 
              position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_jitter(data = gen224_fabae, aes(x = "Fabae", y = Nb_larves), 
              position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(x = "", y = "Nombre de larves", title = "Boxplot pour le génotype 224") +
  theme_minimal() + 
  scale_x_discrete(labels = c("Fabae", "Pisum")) + 
  ylim(0, 40) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, margin = margin(t = 10))) + 
  theme(axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12, margin = margin(r = 10)))


##Mise en forme du tableau : suppression des lignes 11 à 15 pour fabae car NA#### 
gen106_fabae1 <- gen106_fabae[-c(11, 12, 13, 14, 15), ]
gen106_fabae2 <- gen106_fabae1[, -which(names(gen106_fabae1) == "Nb_feuilles")]

merged_data <- rbind(gen106_pisum, gen106_fabae2)

gen224_fabae1 <- gen224_fabae[-c(12, 13, 14, 15), ]
gen224_fabae2 <- gen224_fabae1[, -which(names(gen224_fabae1) == "Nb_feuilles")]

merged_data1 <- rbind(gen224_pisum, gen224_fabae2)


##Analyses stats

###Avec effet aléatoire génotype####
#106
mod4 <- glmmTMB(data=merged_data, Nb_larves ~ Aphid + (1|Genotype/Rep),family=nbinom1(link="log")) 
summary(mod4)
hypotheses3 <- glht(mod4, linfct = mcp(Aphid = "Tukey"))
summary(hypotheses3)

#224
mod5 <- glmmTMB(data=merged_data1, Nb_larves ~ Aphid + (1|Genotype/Rep),family=nbinom1(link="log")) 
summary(mod5) #même résultats qu'avant  
hypotheses4 <- glht(mod5, linfct = mcp(Aphid = "Tukey"))
summary(hypotheses4)



###################################
##Régression entre nombre de larve et nombre de noeuds#### 
reg <- lm(Nb_larves ~ Nb_noeuds, data = mean_pisum)
plot(Nb_noeuds ~Nb_larves, data = mean_pisum)
##pas de lien entre l'âge de la féverole et le nombre de larve 

reg1 <- lm(Nb_larves ~ Nb_feuilles, data = mean_fabae)
plot(Nb_feuilles ~Nb_larves, data = mean_fabae)
##pas de lien entre l'âge de la féverole et le nombre de larve 

summary(reg)
#r^2=0,16
summary(reg1)
#r^2=0,28
###################################