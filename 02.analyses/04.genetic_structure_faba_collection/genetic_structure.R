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

########################################################
## Plot phylogenetic tree based on Jonathan Newick data
## Tutorials
## https://4va.github.io/biodatasci/r-ggtree.html
## https://yulab-smu.top/treedata-book/chapter4.html#basic-tree-visualization
#################################################

faba_tree <- read.tree("01.data/03.genetic_structure_jonathan_dijon/wastral_gene_ratio01_raxml_2.raxml.bestTree")

ggplot(faba_tree) + 
  geom_tree() + 
  theme(text = element_text(size = 3)) +
  theme_tree2()  # The horizontal dimension in this plot shows the amount of genetic change,

ggtree(faba_tree) + theme_tree2() + geom_tiplab(size = 2, color = 'blue') + 
  geom_highlight(node=10, fill="steelblue")

ggsave(filename = "02.analyses/04.genetic_structure_faba_collection/genetic_structure.pdf", width = 16, height = 12)
