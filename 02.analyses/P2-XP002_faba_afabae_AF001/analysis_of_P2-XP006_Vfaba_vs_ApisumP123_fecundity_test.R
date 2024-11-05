rm(list = ls()) # clean the environment
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("lme4"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("emmeans"))
suppressPackageStartupMessages(library("patchwork"))

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
## Import of total larvae counts
################################

# Remove measures where adult alive = 0
apisum_counts_xp6 <- read.csv("01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", 
                         stringsAsFactors = F) %>% 
  inner_join(igepp_ivip_genetic_group, by = "ivip") %>% 
  mutate(cluster = paste("group",cluster,sep = "")) %>% 
  filter(adult_alive != 0) %>% 
  filter(! is.na(larvae_total_number)) %>%
  mutate(short_code = as.character(short_code))

# Check
apisum_counts_xp6 %>% pull(ivip) %>% unique() %>% length()

######
# Plot
######
# Number of biological replicate per genotype (to add to the graph)
n_per_genotype <- 
  apisum_counts_xp6 %>% 
  dplyr::select(ivip, rep) %>%
  group_by(ivip) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")


order_genotypes = 
  apisum_counts_xp6 %>% 
  group_by(ivip) %>% 
  summarise(median_counts = median(larvae_total_number, na.rm=TRUE)) %>% 
  arrange(median_counts) %>% 
  pull(ivip) %>% 
  unique()

p1 <- 
  apisum_counts_xp6 %>% 
  inner_join(n_per_genotype, by = "ivip") %>%
  mutate(ivip = factor(ivip, levels = order_genotypes)) %>%
  ggplot(., aes(x = ivip, y = larvae_total_number, fill = phenotype_to_confirm)) +
  geom_boxplot() +
  geom_point(alpha=0.4) +
  geom_text(aes(label = paste("n=",n_rep,sep = ""), y = 60), vjust = -0.5, hjust = 0.5, size = 4) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Fecundity test (A. pisum P123)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of A. pisum larvae (from one L1)") +
  scale_y_continuous(limits=c(0,60)) +
  theme(axis.text = element_text(size = 10));p1

ggsave(plot = p1, filename = "02.analyses/P2-XP006_Vfaba_ApisumP123_fecundity_test/plot01_fecundity_test_16_genotypes.pdf")


#############################################
# Compare with big screen results (2022-2023)
#############################################

source("02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/analysis01_n_higher_than_6_apisum_.R")

# Combine A. pisum results (big screen 2023, object 'apisum_fecundity_n6') and P2-XP006 (object 'apisum_counts_xp006')

# Number of biological replicate per genotype (to add to the graph)
# Keep only short_code, ivip, larvae_total_number columns in results of P2-XP006
n_per_genotype_xp006 <- 
  apisum_fecundity_n6 %>% 
  dplyr::select(ivip, replicate) %>%
  group_by(ivip) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

apisum_counts_xp6_ivip <-
  apisum_counts_xp6 %>%
  dplyr::select(short_code, ivip, larvae_total_number) %>%
  mutate(year = "2024") %>% 
  inner_join(., n_per_genotype_xp006, by = "ivip")

# Keep only number of replicates for the 16 genotypes also assessed in P2-XP002 (2024)
# Keep only results for the 16 genotypes also assessed in P2-XP002 (2024)
n_per_genotype_2023 <- 
  n_per_genotype %>% 
  filter(ivip %in% apisum_counts_xp6_ivip$ivip)

apisum_fecundity_n6_ivip <- 
  apisum_fecundity_n6 %>% 
  dplyr::select(short_code, ivip, larvae_total_number) %>% 
  mutate(year = "2023") %>% 
  dplyr::filter(ivip %in% apisum_counts_xp6_ivip$ivip) %>% 
  inner_join(., n_per_genotype_2023, by = "ivip")

### Combine results and compare the two experiment (2023 screen and P2-XP002) 
p_comparison_years <- 
  bind_rows(apisum_fecundity_n6_ivip, apisum_counts_xp6_ivip) %>% 
  ggplot(., aes(x = year, y = larvae_total_number, fill = year)) +
  geom_boxplot() +
  geom_point(alpha=0.4) +
  geom_text(aes(label = paste("n=",n_rep,sep = ""), y = 70), vjust = -0.5, hjust = 0.5, size = 4) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Fecundity test (A. pisum P123)") +
  labs(x = "Year", y = "Total number of A. pisum larvae (from one L1)") +
  theme(axis.text = element_text(size = 10)) +
  facet_wrap(~ivip, scales = "fixed") +
  scale_y_continuous(limits=c(0,80)); p_comparison_years

ggsave(plot = p_comparison_years, 
       filename = "02.analyses/P2-XP006_Vfaba_ApisumP123_fecundity_test/plot02_comparison_2023-2024.pdf")

### Combine results and plot all N larvae per genotype (2023 screen and P2-XP002 in 2024 combined)
n_per_genotype_2023_and_2024 <- 
  bind_rows(n_per_genotype_2023, n_per_genotype_xp006) %>% 
  group_by(ivip) %>%
  summarise(n_rep = sum(n_rep))

order_genotypes_2023_and_2024 = 
  bind_rows(apisum_fecundity_n6_ivip, apisum_counts_xp6_ivip) %>% 
  group_by(ivip) %>% 
  summarise(median_counts = median(larvae_total_number, na.rm=TRUE)) %>% 
  arrange(median_counts) %>% 
  pull(ivip) %>% 
  unique()

p_comparison_genotypes <- 
  bind_rows(apisum_fecundity_n6_ivip, apisum_counts_xp6_ivip) %>% 
  dplyr::select(short_code, ivip, larvae_total_number, year) %>%
  inner_join(., n_per_genotype_2023_and_2024, by = "ivip") %>%
  mutate(ivip = factor(ivip, levels = order_genotypes_2023_and_2024)) %>%
  ggplot(., aes(x = ivip, y = larvae_total_number)) +
  geom_boxplot() +
  geom_jitter(width=0.05, alpha=0.4, aes(shape = year), size=2) +
  geom_text(aes(label = paste("n=",n_rep,sep = ""), y = 55), vjust = -0.5, hjust = 0.5, size = 4) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Fecundity test of A. pisum P123 on Vicia faba genotypes (All data: 2022 to 2024)") +
  labs(x = "Vicia faba genotype", y = "Total number of A. pisum larvae (from one L1)") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_continuous(limits=c(0,60)); p_comparison_genotypes

ggsave(plot = p_comparison_genotypes, 
       filename = "02.analyses/P2-XP006_Vfaba_ApisumP123_fecundity_test/plot03_fecundity_test_16_genotypes_2023_and_2024_facet.pdf")

order_genotypes_2023_and_2024_bis = 
  bind_rows(apisum_fecundity_n6_ivip, apisum_counts_xp6_ivip) %>% 
  group_by(short_code) %>% 
  summarise(median_counts = median(larvae_total_number, na.rm=TRUE)) %>% 
  arrange(median_counts) %>% 
  pull(short_code) %>% 
  unique()

p_comparison_genotypes2 <- 
  bind_rows(apisum_fecundity_n6_ivip, apisum_counts_xp6_ivip) %>% 
  dplyr::select(short_code, ivip, larvae_total_number, year) %>%
  inner_join(., n_per_genotype_2023_and_2024, by = "ivip") %>%
  mutate(short_code = factor(short_code, levels = order_genotypes_2023_and_2024_bis)) %>%
  ggplot(., aes(x = short_code, y = larvae_total_number)) +
  geom_boxplot() +
  geom_jitter(width=0.05, alpha=0.4, aes(shape = year), size=2) +
  geom_text(aes(label = paste("n=",n_rep,sep = ""), y = 55), vjust = -0.5, hjust = 0.5, size = 4) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Fecundity test of A. pisum P123 on Vicia faba genotypes (All data: 2022 to 2024)") +
  labs(x = "Vicia faba genotype", y = "Total number of A. pisum larvae (from one L1)") +
  theme(axis.text = element_text(size = 10)) +
  scale_y_continuous(limits=c(0,60)); p_comparison_genotypes2

ggsave(plot = p_comparison_genotypes2, 
       filename = "02.analyses/P2-XP006_Vfaba_ApisumP123_fecundity_test/plot04_fecundity_test_16_genotypes_2023_and_2024_shortcode.pdf")