suppressPackageStartupMessages(library("tidyverse"))
library(multcompView)    
library("ggdist")
library("glmmTMB")
library("multcomp")

################################
## Genetic structure data import
################################

### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
igepp_to_ivip <- read.csv("01.data/03.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv") %>% 
  rename("ivip" = "seedbank_id") %>% 
  rename("genotype" = "igepp_id")

genetic_structure <- read.csv("01.data/03.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv") %>% 
  rename("ivip" = "genotype")

igepp_ivip_genetic_group <- inner_join(igepp_to_ivip, genetic_structure)
rm(igepp_to_ivip, genetic_structure)

################################
## Import of total larvae counts
################################

# Manon's results
pisum_counts <- read.csv("01.data/06.P2_XP001_faba_pisum_P123/most_sensitive_resistant_vs_apisum_P123.csv", 
                         stringsAsFactors = F) %>% 
  inner_join(igepp_ivip_genetic_group) %>% 
  mutate(cluster = paste("group",cluster,sep = ""))

######
# Plot
######
order_genotypes = 
  pisum_counts %>% 
  group_by(genotype) %>% 
  summarise(sum_counts = sum(larvae_total_number, na.rm=TRUE)) %>% 
  arrange(sum_counts) %>% 
  pull(genotype) %>% 
  unique()

p_p123 <- 
  pisum_counts %>% 
  mutate(genotype = factor(genotype, levels = order_genotypes)) %>% 
  ggplot(., aes(x = genotype, y = larvae_total_number, fill = phenotype_class_p123, label = cluster)) +
  geom_boxplot() +
  geom_point(alpha=0.4) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Fecundity test (A. pisum P123)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Acyrtosiphon fabae larvae (from one L1)") +
  geom_text(y = 80, hjust=0) + 
  theme(axis.text = element_text(size = 10))
p_p123
ggsave(plot = p_p123, filename = "02.analyses/06.P2-XP001_pisumP123/P2-XP001_fecundity_plot.png")
ggsave(plot = p_p123, filename = "02.analyses/06.P2-XP001_pisumP123/P2-XP001_fecundity_plot.pdf")


######
# Stat
######

# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-specification
glmm_model = glmmTMB(data = pisum_counts, 
                 formula = larvae_total_number ~ genotype + (1|genotype/replicate)) # nested random effects

# multcomp ghlt
glht(glmm_model, linfct = mcp("genotype" = "Tukey"))

