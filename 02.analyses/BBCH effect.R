library(lme4)

#################################################################################$
###############################. A. fabae . #####################################
#################################################################################


#data preparation ----
data_BBCH_fabae <- read.csv("01.data/P2-XP009_afabae_third_confirmation/P2-XP009_fecundity_data_afabae.csv", sep = ",") %>%
  mutate(adult_alive = as.factor(adult_alive)) %>%
  mutate(genotype = as.character(genotype)) %>%
  filter(adult_alive != "0")
summary(data_BBCH_fabae)
data_BBCH_fabae$bbch_10days<-as.character(data_BBCH_fabae$bbch_10days)
data_BBCH_fabae$bbch_25days<-as.character(data_BBCH_fabae$bbch_25days)


#effect of BBCH on larvae count (no effect)----
model<- lm(larvae_total_number~bbch_10days, data = data_BBCH_fabae)
anova(model)

model<- lm(larvae_total_number~bbch_25days, data = data_BBCH_fabae)
anova(model)


model<- lm(larvae_total_number~bbch_25days + block+ genotype, data = data_BBCH_fabae)
anova(model)

model<- lm(larvae_total_number~bbch_10days + block+ genotype, data = data_BBCH_fabae)
anova(model)

#boxplot (small differences but no clear correlation btw BBCH and larvae count most prevosky have a bbch of 13 therefore there is a bias----

data_BBCH_fabae %>% 
  filter(bbch_10days %in% c(9, 10, 10.5, 11, 12)) %>%  # Filter out values that are not in the specified list
  mutate(bbch_10days = factor(bbch_10days, levels = c(9, 10, 10.5, 11, 12))) %>%  # Order the 'bbch_10days' factor
  ggplot(aes(x = bbch_10days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 10 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


data_BBCH_fabae %>% 
  filter(bbch_25days %in% c(13, 14, 15, 16, 17)) %>%  # Filter values 13, 14, 15, 16, 17 for 'bbch_25days'
  mutate(bbch_25days = factor(bbch_25days, levels = c(13, 14, 15, 16, 17))) %>%  # Order the 'bbch_25days' factor
  ggplot(aes(x = bbch_25days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 25 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#to remove the bias we can extract the emmeans with the genotype as random factor ----

# model evaluation before emmeans calculation to get the ajusted means for BBCH (best model = genotype as random factor)----
summary(data_BBCH_fabae)
model0<- lm(larvae_total_number~1, data = data_BBCH_fabae)
model1<- lmer(larvae_total_number~bbch_25days + (1 |block), data = data_BBCH_fabae)
model2<- lmer(larvae_total_number~bbch_25days + (1 |genotype), data = data_BBCH_fabae)
model3<- lmer(larvae_total_number~bbch_25days + (1 |block)+(1 |genotype), data = data_BBCH_fabae)
model1b<- lmer(larvae_total_number~bbch_10days + (1 |block), data = data_BBCH_fabae)
model2b<- lmer(larvae_total_number~bbch_10days + (1 |genotype), data = data_BBCH_fabae)
model3b<- lmer(larvae_total_number~bbch_10days + (1 |block)+(1 |genotype), data = data_BBCH_fabae)

compare_performance(model0, model1, model2,model3,model1b, model2b,model3b, rank = FALSE, verbose = TRUE)

#emmeans (no difference when we remove the genotype effect) ----
model2<- lmer(larvae_total_number~bbch_10days + (1 |genotype), data = data_BBCH_fabae)
model2b<- lmer(larvae_total_number~bbch_25days + (1 |genotype), data = data_BBCH_fabae)

emm_genotypeBBCH10 <- emmeans(model2, ~ bbch_10days)
emm_dfBBCH10 <- as.data.frame(emm_genotypeBBCH10)

emm_genotypebbch25 <- emmeans(model2b, ~ bbch_25days)
emm_dfBBCH25 <- as.data.frame(emm_genotypebbch25)




#################################################################################$



###############################. pisum. #####################################
#################################################################################

#data preparation ----
data_BBCH_pisum08 <- read.csv("01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", sep = ",") %>%
  mutate(adult_alive = as.factor(adult_alive)) %>%
  mutate(genotype = as.character(genotype)) %>%
  filter(adult_alive != "0")
summary(data_BBCH_pisum08)
data_BBCH_pisum08$bbch_10days<-as.character(data_BBCH_pisum08$bbch_10days)
data_BBCH_pisum08$bbch_25days<-as.character(data_BBCH_pisum08$bbch_25days)


data_BBCH_pisum10 <- read.csv("01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", sep = ",") %>%
  mutate(adult_alive = as.factor(adult_alive)) %>%
  mutate(genotype = as.character(genotype)) %>%
  filter(adult_alive != "0")
data_BBCH_pisum10$bbch_10days<-as.character(data_BBCH_pisum10$bbch_10days)
data_BBCH_pisum10$bbch_25days<-as.character(data_BBCH_pisum10$bbch_25days)
summary(data_BBCH_pisum10)

data_BBCH_pisum10 <- data_BBCH_pisum10 %>%
  filter(bbch_25days != "10") # 1 genotype removed, I will need to study it later


combined_data <- bind_rows(data_BBCH_pisum08, data_BBCH_pisum10)

#effect of BBCH on larvae count ----
model<- lm(larvae_total_number~bbch_10days, data = data_BBCH_pisum08)
anova(model)

model<- lm(larvae_total_number~bbch_25days, data = data_BBCH_pisum08)
anova(model)

model<- lm(larvae_total_number~bbch_10days, data = data_BBCH_pisum10)
anova(model)

model<- lm(larvae_total_number~bbch_25days, data = data_BBCH_pisum10)
anova(model)

model<- lm(larvae_total_number~bbch_10days, data = combined_data)
anova(model)

model<- lm(larvae_total_number~bbch_25days, data = combined_data)
anova(model)

#boxplot (small differences but no clear correlation btw BBCH and larvae count most prevosky have a bbch of 13 therefore there is a bias----

data_BBCH_pisum08 %>% 
  filter(bbch_10days %in% c(9, 10, 10.5, 11, 12)) %>%  # Filter out values that are not in the specified list
  mutate(bbch_10days = factor(bbch_10days, levels = c(9, 10, 10.5, 11, 12))) %>%  # Order the 'bbch_10days' factor
  ggplot(aes(x = bbch_10days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 10 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


data_BBCH_pisum08 %>% 
  filter(bbch_25days %in% c(13, 14, 15, 16, 17)) %>%  # Filter values 13, 14, 15, 16, 17 for 'bbch_25days'
  mutate(bbch_25days = factor(bbch_25days, levels = c(13, 14, 15, 16, 17))) %>%  # Order the 'bbch_25days' factor
  ggplot(aes(x = bbch_25days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 25 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

data_BBCH_pisum10 %>% 
  filter(bbch_10days %in% c(9, 10, 10.5, 11, 12)) %>%  # Filter out values that are not in the specified list
  mutate(bbch_10days = factor(bbch_10days, levels = c(9, 10, 10.5, 11, 12))) %>%  # Order the 'bbch_10days' factor
  ggplot(aes(x = bbch_10days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 10 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


data_BBCH_pisum10 %>% 
  filter(bbch_25days %in% c(13, 14, 15, 16, 17)) %>%  # Filter values 13, 14, 15, 16, 17 for 'bbch_25days'
  mutate(bbch_25days = factor(bbch_25days, levels = c(13, 14, 15, 16, 17))) %>%  # Order the 'bbch_25days' factor
  ggplot(aes(x = bbch_25days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 25 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability





combined_data %>% 
  filter(bbch_10days %in% c(9, 10, 10.5, 11, 12)) %>%  # Filter out values that are not in the specified list
  mutate(bbch_10days = factor(bbch_10days, levels = c(9, 10, 10.5, 11, 12))) %>%  # Order the 'bbch_10days' factor
  ggplot(aes(x = bbch_10days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 10 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


combined_data %>% 
  filter(bbch_25days %in% c(13, 14, 15, 16, 17)) %>%  # Filter values 13, 14, 15, 16, 17 for 'bbch_25days'
  mutate(bbch_25days = factor(bbch_25days, levels = c(13, 14, 15, 16, 17))) %>%  # Order the 'bbch_25days' factor
  ggplot(aes(x = bbch_25days, y = larvae_total_number)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Total Number by IVIP (A. pisum P123)",
    x = "IVIP (BBCH 25 Days)",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability





#emmeans even when we remove genetic effect it seem like we have an effect of the BBCH 25 on the larvae count ...----
model2<- lmer(larvae_total_number~bbch_10days + (1 |genotype), data = combined_data)
model2b<- lmer(larvae_total_number~bbch_25days + (1 |genotype), data = combined_data)

emm_genotypeBBCH10 <- emmeans(model2, ~ bbch_10days)
emm_dfBBCH10 <- as.data.frame(emm_genotypeBBCH10)

emm_genotypebbch25 <- emmeans(model2b, ~ bbch_25days)
emm_dfBBCH25 <- as.data.frame(emm_genotypebbch25)


