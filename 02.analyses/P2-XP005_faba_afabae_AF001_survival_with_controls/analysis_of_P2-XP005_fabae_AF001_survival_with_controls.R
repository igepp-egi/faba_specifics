suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("multcomp"))    
suppressPackageStartupMessages(library("MASS"))    


################################
## Genetic structure data import
################################

### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
igepp_to_ivip <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv") %>% 
  rename("ivip" = "seedbank_id") %>% 
  rename("genotype" = "igepp_id")

genetic_structure <- read.csv("00.info/02.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv") %>% 
  rename("ivip" = "genotype")

igepp_ivip_genetic_group <- inner_join(igepp_to_ivip, genetic_structure)
rm(igepp_to_ivip, genetic_structure)

################################
## Import of survival data
################################

fabae_counts <- read.csv("01.data/P2-XP005_fabae_afabae_survival_with_controls/P2-XP005_survival_assay.csv", 
                         stringsAsFactors = F) %>% 
  left_join(., igepp_ivip_genetic_group) %>% 
  mutate(cluster = paste("group",cluster,sep = "")) %>% 
  mutate(short_code = as.character(short_code)) %>% 
  mutate(genotype = paste(ivip_code, short_code, sep = "_")) %>%
  # add a column labeled with test and controls if the column nb_zero_afabae contains 0/ string
  mutate(test_control = ifelse(nb_zero_afabae == "0/12", "control", "test")) %>%
  arrange(desc(test_control))

good_genotype_order <- fabae_counts %>% pull(genotype) %>% unique()

fabae_counts$genotype <- factor(fabae_counts$genotype, levels = good_genotype_order)
  
######
# Plot
######

p_adults <- 
  fabae_counts %>% 
  ggplot(., aes(x = genotype, y = nb_adults_observed, fill = test_control)) +
  geom_boxplot() +
  geom_jitter(width = 0.01, alpha=0.5) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Survival test (A. fabae AF001)") +
  labs(x = "Vicia faba genotype", y = "Number of A. fabae adults (from 10 larvae)") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,5,10)) +
  coord_flip() +
  # create a y axis numbered from 0 to 10 with 1 unit intervals
  scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0,10)) +
  geom_hline(yintercept = 10, colour = "blue", linetype = "dashed"); p_adults

ggsave(plot = p_adults, 
       filename = "02.analyses/P2-XP005_faba_afabae_AF001_survival_with_controls/adult_survival_plot.png")
ggsave(plot = p_adults, 
       filename = "02.analyses/P2-XP005_faba_afabae_AF001_survival_with_controls/adult_survival_plot.pdf")


############################################
# Stats: GLM and GLMM models and comparisons
############################################

# GLM models: no random effects
glm_null_model <- glm.nb(formula = nb_adults_observed ~ 1, 
                         data = fabae_counts)
glm_geno_model <- glm.nb(formula = nb_adults_observed ~ genotype, 
                         data = fabae_counts)


# GLMM models: with biological replicate (plant) as random effect
glmer_null_model <- lme4::glmer.nb(formula = nb_adults_observed ~ 1 + (1|biol_replicate), 
                                   data = fabae_counts,  
                                   control = lme4::glmerControl(optimizer = "bobyqa"))
glmer_geno_model <- lme4::glmer.nb(formula = nb_adults_observed ~ genotype + (1|biol_replicate), 
                                   data = fabae_counts,  
                                   control = lme4::glmerControl(optimizer = "bobyqa"))

anova(glm_null_model, glm_geno_model)
anova(glmer_null_model, glmer_geno_model)

# Conclusions: no difference between random model and genotype model.

#############################
# Stats: post-hoc comparisons
#############################

library("emmeans")
pairs(emmeans(glm_geno_model, "genotype")) %>% 
  broom::tidy() %>%
  # filter rows that contains string 'Castell'
  dplyr::filter(., grepl("Castell",contrast))

library("emmeans")
pairs(emmeans(glmer_geno_model, "genotype")) %>% 
  broom::tidy() %>%
  # filter rows that contains string 'Castell'
  dplyr::filter(., grepl("Castell",contrast))

# Conclusions: no difference between genotypes


