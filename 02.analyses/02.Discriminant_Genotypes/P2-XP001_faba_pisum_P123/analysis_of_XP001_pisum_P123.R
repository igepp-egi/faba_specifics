suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("glmmTMB"))
suppressPackageStartupMessages(library("multcomp"))

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
pisum_counts <- read.csv("01.data/06.P2_XP001_faba_pisum_P123/P2-XP001_fecundity_results.csv", 
                         stringsAsFactors = F) %>% 
  inner_join(igepp_ivip_genetic_group) %>% 
  mutate(cluster = paste("group",cluster,sep = ""))

n_rep_per_genotype =
  pisum_counts %>% 
  dplyr::mutate(genotype_code = paste(genotype, short_code, sep = "_")) %>% 
  group_by(genotype_code) %>% 
  tally(name = "larvae_total_count")

######
# Plot
######
order_genotypes = 
  pisum_counts %>% 
  dplyr::mutate(genotype_code = paste(genotype, short_code, sep = "_")) %>% 
  group_by(genotype_code) %>% 
  summarise(sum_counts = sum(larvae_total_count, na.rm=TRUE)) %>% 
  arrange(sum_counts) %>%
  pull(genotype_code) %>% 
  unique()

p_p123 <- 
  pisum_counts %>% 
  dplyr::mutate(genotype_code = paste(genotype, short_code, sep = "_")) %>% 
  mutate(genotype_code = factor(genotype_code, levels = order_genotypes)) %>% 
  ggplot(., aes(x = genotype_code, y = larvae_total_count, fill = phenotype_class, label = cluster)) +
  geom_boxplot() +
  geom_point(alpha=0.4) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Fecundity test (A. pisum P123)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Acyrtosiphon fabae larvae (from one L1)") +
  geom_text(y = 14, hjust=0) + 
  scale_y_continuous(limits=c(0,15)) + 
  theme(axis.text = element_text(size = 10)) +
  scale_fill_brewer(type = 'qual', palette = 2) 
p_p123
ggsave(plot = p_p123, filename = "02.analyses/06.P2-XP001_pisumP123/P2-XP001_fecundity_plot.png")
ggsave(plot = p_p123, filename = "02.analyses/06.P2-XP001_pisumP123/P2-XP001_fecundity_plot.pdf")


######
# Stat
######

# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-specification
glmm_model = glmmTMB(data = pisum_counts, 
                 formula = larvae_total_count ~ genotype + (1|genotype/replicate)) # nested random effects

suppressPackageStartupMessages(library("lme4"))
#glmm_model = glmer.nb(formula=larvae_total_count ~ genotype + (1|replicate), data = pisum_counts)

# multcomp ghlt
glht(glmm_model, linfct = mcp(genotype = "Tukey"))
