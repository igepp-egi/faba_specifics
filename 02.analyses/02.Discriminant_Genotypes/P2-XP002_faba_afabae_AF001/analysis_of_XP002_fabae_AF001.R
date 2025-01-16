suppressPackageStartupMessages(library("tidyverse"))
library(multcompView)    
library("ggdist")

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

fabae_counts <- read.csv("01.data/07.P2_XP002_faba_afabae_AF001/most_sensitive_most_resistant_fabae.csv", 
                         stringsAsFactors = F) %>% 
  inner_join(igepp_ivip_genetic_group) %>% 
  mutate(cluster = paste("group",cluster,sep = "")) %>% 
  filter(larvae_total_number > 0) %>% 
  mutate(short_code = as.character(short_code)) %>% 
  arrange(short_code)

######
# Plot
######
order_genotypes = 
  fabae_counts %>% 
  group_by(genotype) %>% 
  summarise(median_counts = median(larvae_total_number, na.rm=TRUE)) %>% 
  arrange(median_counts) %>% 
  pull(genotype) %>% 
  unique()

p_fabae <- 
  fabae_counts %>% 
  ggplot(., aes(x = short_code, y = larvae_total_number, fill = phenotype_class_af001, label = cluster)) +
  geom_boxplot() +
  geom_point(alpha=0.4) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Fecundity test (A. fabae AF001)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Aphis fabae larvae (from one L1)") +
  geom_text(y = 50, hjust=0) +
  scale_y_continuous(limits=c(0,60)) +
  theme(axis.text = element_text(size = 10))
p_fabae
ggsave(plot = p_fabae, filename = "02.analyses/07.P2-XP002_faba_afabae_AF001/P2-XP002_fabae_fecundity_plot.png")
ggsave(plot = p_fabae, filename = "02.analyses/07.P2-XP002_faba_afabae_AF001/P2-XP002_fabae_fecundity_plot.pdf")


#######
# Stats
#######

