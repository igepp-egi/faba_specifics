suppressPackageStartupMessages(library("tidyverse"))
library(multcompView)    
library("ggdist")

################################
## Genetic structure data import
################################

### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
igepp_to_ivip <- read.csv("03.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv") %>% 
  rename("ivip" = "seedbank_id") %>% 
  rename("genotype" = "igepp_id")

genetic_structure <- read.csv("03.info/02.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv") %>% 
  rename("ivip" = "genotype")

igepp_ivip_genetic_group <- inner_join(igepp_to_ivip, genetic_structure)
rm(igepp_to_ivip, genetic_structure)

################################
## Import of survival data
################################

fabae_counts <- read.csv("01.data/P2-XP003_fabae_afabae_survival/P2-XP003_survival_assay_faba_ten_larvae_AF001.csv", 
                         stringsAsFactors = F) %>% 
  inner_join(igepp_ivip_genetic_group) %>% 
  mutate(cluster = paste("group",cluster,sep = "")) %>% 
  mutate(short_code = as.character(short_code)) %>% 
  mutate(genotype = paste(genotype, short_code, sep = "_"))

######
# Plot
######

p_adults <- 
  fabae_counts %>% 
  ggplot(., aes(x = genotype, y = adult_alive, fill = genotype)) +
  geom_boxplot() +
  geom_point(alpha=0.4) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Survival test (A. fabae AF001)") +
  labs(x = "Vicia faba genotype", y = "Number of A. fabae adults (from 10 larvae)") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,5,10)) +
  geom_hline(yintercept = 10, colour = "blue", linetype = "dashed")
p_adults
ggsave(plot = p_fabae, filename = "02.analyses/P2-XP003_faba_afabae_AF001_survival/adult_survival_plot.png")
ggsave(plot = p_fabae, filename = "02.analyses/P2-XP003_faba_afabae_AF001_survival/adult_survival_plot.pdf")


#######
# Stats
#######

# glm_model <- lme4::glmer.nb(formula = adult_alive ~ genotype + (1|genotype:replicate), data = fabae_counts)
# multcomp::glht(model = glm_model, linfct = mcp(genotype = "Tukey"))
# 
# library("emmeans")
# pairs(emmeans(glm_model, "genotype"))

adult_counts =
  fabae_counts %>% 
  group_by(genotype) %>% 
  summarise(adults_average = mean(adult_alive, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(adults_expected = 10)

