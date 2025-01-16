suppressPackageStartupMessages(library('ggtree'))
suppressPackageStartupMessages(library('tidytree'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('treeio'))

#################
## Genetic groups
#################

genetic_groups <- read.csv("01.data/03.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv")

#################################################
## Import Jonathan Newick data
## https://4va.github.io/biodatasci/r-ggtree.html
#################################################

faba_tree <- read.newick("01.data/03.genetic_structure_jonathan_dijon/wastral_gene_ratio01_raxml_2.raxml.bestTree")

faba_tree_tibble <- as_tibble(faba_tree) %>% as.phylo()

########################################################
## Plot phylogenetic tree based on Jonathan Newick data
## Tutorials
## https://4va.github.io/biodatasci/r-ggtree.html
## https://yulab-smu.top/treedata-book/chapter4.html#basic-tree-visualization
#################################################

ggplot(faba_tree) + 
  geom_tree() + 
  theme(text = element_text(size = 3)) +
  theme_tree2()  # The horizontal dimension in this plot shows the amount of genetic change,

ggtree(faba_tree) + theme_tree2() + geom_tiplab(size = 2, color = 'blue') 

ggsave(filename = "02.analyses/04.genetic_structure_faba_collection/genetic_structure.pdf", width = 16, height = 12)

####################################################################
### Highlight selected genotypes: sensitive or resistant to A. pisum
####################################################################

apisum_genotypes <- read.csv("02.analyses/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/apisum_most_contrasted_genotypes.csv")
afabae_genotypes <- read.csv("02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/most_extreme_genotypes_aphis_fabae.csv")

source("02.analyses/04.genetic_structure_faba_collection/get_igepp_to_ivip_correspondence.R")
igepp_ivip_genetic_group2 <- 
  igepp_ivip_genetic_group %>% 
  dplyr::rename("label" = "ivip")

# faba_tree_annotated <- 
#   faba_tree_tibble %>% 
#   full_join(., igepp_ivip_genetic_group2, by = "label") %>% 
#   as.treedata() %>% 
#   as.phylo()  












# # to get the node number (useful for coloring it)
# ggtree(faba_tree) + geom_nodelab(aes(label = node),  hjust = 0) 


