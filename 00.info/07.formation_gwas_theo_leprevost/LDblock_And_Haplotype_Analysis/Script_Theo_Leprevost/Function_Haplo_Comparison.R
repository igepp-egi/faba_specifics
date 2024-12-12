########### Author: Theo Leprevost
########### Date : 09/25/2023
### Haplotype construction
#Input data : A Haplo_Detail dataframe produced by the R function Function_Haplo_Construction.R
#Output data : An excel table with the statistical results and a haplo graph

Haplo_Comparison<-function(data,file,colored_geno,variable) {

#Load the libraries
library(car)
library(RVAideMemoire)
library(dplyr)
library(multcomp)
library(multcompView)
library(emmeans)
library(ggplot2)
library(writexl)
library(data.table)
library(viridis)
library(simplecolors)
library(MuMIn)

#Trim the dataset
data<-data[,c(1,2,ncol(data)-3,ncol(data)-2,ncol(data)-1,ncol(data))]

# Keep the haplotypes showing a frequency above or equal to 5%
filtered_haplotypes <- data[data$Percentage_of_Lines >= 5, ]

# Keep the rare haplotypes
rare_haplotypes <- subset(data, Heterozygous_Haplotype == FALSE & Percentage_of_Lines < 5)
# Get unique rows based on Haplotype_ID and cumulate the number and percentage of lines for each rare haplotype 
sum_rare_haplo_percentage <- sum(distinct(rare_haplotypes, Haplotype_ID, .keep_all = TRUE)$Percentage_of_Lines)
sum_rare_haplo_line <- sum(distinct(rare_haplotypes, Haplotype_ID, .keep_all = TRUE)$Number_of_Lines)

# Add the rare haplotype if the percentage of shared rare haplotype is above 5%
if (nrow(rare_haplotypes)>0 & sum_rare_haplo_percentage>=5){
  rare_haplotypes$Haplotype_ID<-"R"
  rare_haplotypes$Percentage_of_Lines<-sum_rare_haplo_percentage
  rare_haplotypes$Number_of_Lines<-sum_rare_haplo_line
  filtered_haplotypes<-rbind(filtered_haplotypes,rare_haplotypes)
}

### Haplotype comparison
LM <- lm(filtered_haplotypes[,2]~Haplotype_ID,data=filtered_haplotypes)
R2<-data.frame(summary(LM)$adj.r.squared,r.squaredGLMM(LM))
SUM <- data.table(coef(summary(LM)), keep.rownames = 'term')
ANOVA <- anova(LM)
EM<-emmeans(LM,~Haplotype_ID,data=filtered_haplotypes)
CLD <- cld(EM,alpha=0.05,sort=TRUE,Letters=c(letters))
CLD<-data.frame(CLD)
CLD<-CLD[order(CLD[,1],decreasing=F), ]
CLD$.group <- gsub(" ", "", CLD$.group)

# Create a new column to specify fill colors based on haplotype
CLD <- CLD %>%
  mutate(Fill_Color = case_when(
    any(grepl("ab", .group)) | all(.group == "a") ~ "white",  # If there's "ab" or all are "a"
    any(.group == "a") & any(grepl("b", .group)) & !any(grepl("ab", .group)) ~ if_else(.group == "a", "green2", "red1"),  # If "a" is present and there are other groups with "b" but not "ab"
    TRUE ~ "red2"  # All other haplotypes
  ))

# Write the tables using xlsx
write_xlsx(list(SUM.LM = SUM,R2.LM=R2, ANOVA = ANOVA, CLD = CLD),
           paste(file,"_Haplo_Comparison.xlsx",sep=""),
           col_names = TRUE)

# Combine "Code" and "Number_of_Lines" columns into a new label column
filtered_haplotypes$XLabel <- paste(filtered_haplotypes$Haplotype_ID, "\n(N=", filtered_haplotypes$Number_of_Lines,")",sep="")

# Using aggregate
max_values <- aggregate(as.numeric(eval(parse(text =variable))) ~ Haplotype_ID, data = filtered_haplotypes, FUN = max)
colnames(max_values)<-c("Haplotype_ID","max_values")
N<-nrow(max_values)

# Specify genotype fill color
filtered_haplotypes$Dot_Color <- ifelse(filtered_haplotypes$Code %in% colored_geno, "#00FFFF", "black")

# Add the number of line in each haplotype
Haplo_plot<- ggplot(filtered_haplotypes, aes(x=XLabel, y=as.numeric(eval(parse(text = (paste0("filtered_haplotypes$",variable,sep=""))))), fill=Dot_Color)) +
  geom_boxplot(color="black", fill=sc(CLD$Fill_Color),width = 0.6,outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=1) +
  geom_jitter(width = 0.1, size=0.8,alpha = 1,color=filtered_haplotypes$Dot_Color) + 
  guides(fill = guide_legend(override.aes = list(color = "black", size = 3))) +
  ylim(min(na.omit(as.numeric(eval(parse(text = (paste0("filtered_haplotypes$",variable,sep="")))))))-2, max(na.omit(as.numeric(eval(parse(text = (paste0("filtered_haplotypes$",variable,sep="")))))))+4) + # Set the y-axis limits# Add jittered point
  annotate("text", x = 1:N, y = round(as.numeric(max_values$max_values),0)+3, label = CLD$.group, vjust = 1.5, size = 4) +
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab(paste(col,"_Haplo_",get(paste("GWAS_Results_",variable,sep=""))[k,1],sep=""))+
  ylab(paste("EMMs_",variable,sep=""))

ggsave(paste(file,"_Haplo_Comparison.png",sep=""),
       Haplo_plot,
       width = 6,
       height = 5,
       dpi = 1200)
}
