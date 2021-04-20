#!/usr/bin/env Rscript
# Author: Jiqiu Wu
# Contact: wujiqiucau@gmail.com
# This is an integrated R script to plot figures of the CRISPR study!
# Enjoy  ~

# Throw everything away and restart
rm(list=ls())
graphics.off()

# load libraries
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ape)
library(ggsci)
library(vegan)
library(Rmisc)


###########################
########## Fig 1 ##########
###########################

fig1_data <- read.csv("../Files/Crispr_meta_quanlity.csv", header = T, sep = ",")
CRIS_data <- subset(fig1_data, month != "M")

CRIS_data$month <- gsub("B", "0",  CRIS_data$month)
CRIS_data$month <- gsub("4M", "4",  CRIS_data$month)
CRIS_data$month <- gsub("12M", "12",  CRIS_data$month)

# Fig 1a
plot_crispr <- ggplot(transform(CRIS_data, Months = factor(month, levels = c("0", "4", "12")))) + 
  geom_violin(aes(x = gender, y = num_CRISPRs, fill = gender))   + facet_wrap( ~ Months) +
  geom_boxplot(aes(x = gender, y = num_CRISPRs), width = 0.1) + 
  ylab("Number of CRISPRs") + xlab("Age (month)") + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig1/a_crispr_num.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_crispr
graphics.off()

# Fig 1b
plot_spacer <- ggplot(transform(CRIS_data, Months = factor(month, levels = c("0", "4", "12")))) + 
  geom_violin(aes(x = gender, y = spacers_sample, fill = gender))   + facet_wrap( ~ Months) +
  geom_boxplot(aes(x = gender, y = spacers_sample), width = 0.1) + 
  ylab("Number of spacers") + xlab("Age (month)") + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig1/b_spacer_num.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_spacer
graphics.off()

# Fig 1c
plot_crispr_density <- ggplot(transform(CRIS_data, Months = factor(month, levels = c("0", "4", "12")))) + 
  geom_violin(aes(x = gender, y = num_CRISPRs/total_bp * 1000000, fill = gender))   + facet_wrap( ~ Months) +
  geom_boxplot(aes(x = gender, y = num_CRISPRs/total_bp * 1000000 ), width = 0.1) + 
  ylab("Density of CRISPRs") + xlab("Age (month)") + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none")

pdf("../Figures/Fig1/c_crispr_density.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_crispr_density
graphics.off()

# Fig 1d
plot_spacer_density <- ggplot(transform(CRIS_data, Months = factor(month, levels = c("0", "4", "12")))) + 
  geom_violin(aes(x = gender, y = spacers_sample/num_CRISPRs, fill = gender))   + facet_wrap( ~ Months) +
  geom_boxplot(aes(x = gender, y = spacers_sample/num_CRISPRs), width = 0.1) + 
  ylab("Number of spacers in one CRISPR") + xlab("Age (month)") + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none")

pdf("../Figures/Fig1/d_spacer_density.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_spacer_density
graphics.off()



###########################
########## Fig 2 ##########
###########################
# fig 2a
richness_data <- read.csv("../Files/bacterial_richness.csv", header = T, sep = ",")
richness2com <- subset(richness_data, Age != "M")
cris_data2com <- select(CRIS_data, baby_id, gender)
colnames(cris_data2com)[1] <- "BabyID"
cris_data2com <- unique(cris_data2com)
richness2plot <- left_join(richness2com, cris_data2com, by = "BabyID")

richness2plot$Age <- gsub("B", "0",  richness2plot$Age)
richness2plot$Age <- gsub("4M", "4",  richness2plot$Age)
richness2plot$Age <- gsub("12M", "12",  richness2plot$Age)

plot_richeness <- ggplot(transform(richness2plot, Months = factor(Age, levels = c("0", "4", "12")))) + 
  geom_violin(aes(x = gender, y = Bacterial_richness, fill = gender))   + facet_wrap( ~ Months) +
  geom_boxplot(aes(x = gender, y = Bacterial_richness), width = 0.1) + 
  ylab("Bacterial richness") + xlab("Age (month)") + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig2/a_richness.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_richeness 
graphics.off()


# fig 2b

community_size <- ggplot(transform(CRIS_data, Months = factor(month, levels = c("0", "4", "12")))) +
  geom_violin(aes(x = gender, y = total_bp / 1000000 , fill = gender)) + facet_wrap( ~ Months) + 
  geom_boxplot(aes(x = gender, y = total_bp / 1000000), width = 0.1) +
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) +
  ylab("Community size (Mbp)") + xlab("Age (month)")  + scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none") 


pdf("../Figures/Fig2/b_community_size.pdf", width = 4.3/2.54, height = 4.6/2.54)
community_size
graphics.off()

# fig 2c
cris_model <- select(CRIS_data, baby_id, month, num_CRISPRs, spacers_sample)
colnames(cris_model)[1:2] <- c("BabyID", "Age")
data_model <- left_join(cris_model, richness2plot, by = c("BabyID", "Age"))

data_model$Age <- gsub("12", "12 months",  data_model$Age)
data_model$Age <- gsub("4", "4 months",  data_model$Age)
data_model$Age <- gsub("0", "At birth",  data_model$Age)

plot_model_cris <- ggplot(data = data_model, aes(x = Bacterial_richness, y = num_CRISPRs, colour = Age)) + geom_point(size = 0.5) + 
  geom_smooth(aes(group = 1) , method = 'lm')  + 
  ylab("The number of CRISPRs") + xlab("Bacterial richness") + 
  theme_grey(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_colour_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none")

pdf("../Figures/Fig2/c_model_cris.pdf", width = 3.2/2.54, height = 4.6/2.54)
plot_model_cris 
graphics.off()

# fig 2d
plot_model_spacer <- ggplot(data = data_model, aes(x = Bacterial_richness, y = spacers_sample, colour = Age)) + geom_point(size = 0.5) + 
  geom_smooth(aes(group = 1) , method = 'lm')  + 
  ylab("The number of spacers") + xlab("Bacterial richness") + 
  theme_grey(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + scale_colour_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  labs(fill = "Age")

pdf("../Figures/Fig2/d_model_spacer.pdf", width = 5.4/2.54, height = 4.6/2.54)
plot_model_spacer 
graphics.off()



###########################
########## Fig 3 ##########
###########################
# Fig 3b

indiv_shared <- read.csv("../Files/individual_shared.csv", header = F, sep = ",")
colnames(indiv_shared) <- c("baby_id", "B", "M4", "M12","B_4M", "B_12M", "b4M_12M",
                            "B_4M_12M")

indi_shared <- melt(indiv_shared, id = "baby_id")
colnames(indi_shared)[2:3] <- c("Month", "spacers")


indi_shared$Month <- gsub("B_4M_12M", "At birth & 4 months & 12 months",  indi_shared$Month)
indi_shared$Month <- gsub("b4M_12M", "4 months & 12 months",  indi_shared$Month)
indi_shared$Month <- gsub("B_12M", "At birth & 12 months",  indi_shared$Month)
indi_shared$Month <- gsub("B_4M", "At birth & 4 months",  indi_shared$Month)
indi_shared$Month <- gsub("M12", "12 months",  indi_shared$Month)
indi_shared$Month <- gsub("M4", "4 months", indi_shared$Month)
indi_shared$Month <- gsub("B", "At birth", indi_shared$Month)

indi_shared2plot <- subset(indi_shared, spacers != 0)

plot_share <- ggplot(data = indi_shared2plot, aes(x = Month, y = spacers, colour = Month)) +
                       geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.1) + ylab("Intersection number") + 
                       scale_x_discrete(limits = c("At birth", "4 months", "12 months",
                                  "At birth & 4 months","At birth & 12 months",
                                  "4 months & 12 months", "At birth & 4 months & 12 months")) + theme_classic(base_size = 8) +
                       theme(axis.text.x =element_text(angle = 60, hjust =1, vjust = 1))  +  
  scale_colour_manual(values = c("#C35C6A", "#88ABDA", "#422256", "#9EBC19","#422256","#422256","#422256"))  +
  xlab("") + ylab("The numbers of spacers")

pdf("../Figures/Fig3/b_indi_shared.pdf", width = 9.2/2.54, height = 7/2.54)
plot_share
graphics.off()


###########################
########## Fig 4 ##########
###########################
# fig 4a
contig_species <- read.csv("../Files/contig_species.csv", header = F, sep = ",")
colnames(contig_species) <- c("baby_id", "month", "contig_id", "acc_num", "phylum", 
                              "genus", "species")
genus_data <- contig_species %>% group_by(genus) %>% summarise(n()) %>% as.data.frame()
colnames(genus_data) <- c("genus","value")

phylum_data <- contig_species %>% group_by(phylum) %>% summarise(n()) %>% as.data.frame()
colnames(phylum_data) <- c("phylum","value")

phylum_data = phylum_data[order(phylum_data$value, decreasing = TRUE),]
myLabel = as.vector(phylum_data$phylum)   
myLabel = paste(myLabel, "(", round(phylum_data$value / sum(phylum_data$value) * 100, 2), "%)", sep = "")   

colnames(phylum_data)[1] <- "Phylum"

phylum_plot <-  ggplot(phylum_data, aes(x = "", y = value, fill = Phylum)) +
  geom_bar(stat = "identity") +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values = c("#89c3eb", "#c97586", "#F5DF3D", "#93CA76", "#9D9D9D")) + 
  theme_classic(8) + theme(panel.border = element_rect( fill = NA))

pdf("../Figures/Fig4/a_phylum_summary.pdf", width = 7/2.54, height = 6/2.54)
phylum_plot
graphics.off()

# fig 4b
genus_data = genus_data[order(genus_data$value, decreasing = TRUE),]

genus2plot = subset(genus_data, value > 5)

genus_plot <-  ggplot(genus2plot, mapping = aes(x = reorder(genus, -value), y = value)) +
  geom_bar(stat = "identity") + 
  theme_classic(8) +
  theme(axis.text.x = element_text(angle=60, vjust = 0.5))   + 
  theme(panel.border = element_rect( fill = NA)) +
  labs(x = "Genus", y = "Number of CRISPRs", title = "") 


pdf("../Figures/Fig4/b_genus_summary.pdf", width = 10.8/2.54, height = 6/2.54)
genus_plot
graphics.off()



###########################
########## Fig 6 ##########
###########################

# fig 6a
acquire_data_1 <- read.csv("../Files/acquired_1.csv", header = F, sep = ",")
colnames(acquire_data_1) <- c("Changes","value")
acquire_data_1$Changes <- c("No changes", "Changed")

acquire_1 <-  ggplot(acquire_data_1, aes(x = "", y = value, fill = Changes)) +
  geom_bar(stat = "identity") +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.position = "top") + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values = c("#89c3eb", "#c97586")) + 
  theme_classic(8) + theme(panel.border = element_rect( fill = NA))


pdf("../Figures/Fig6/a_acquire_1.pdf", width = 6/2.54, height = 6/2.54)
acquire_1
graphics.off()

# fig 6b
acquire_data_2 <- read.csv("../Files/acquired_2.csv", header = F, sep = ",")
colnames(acquire_data_2) <- c("changes","value")

acquire_2 <-  ggplot(acquire_data_2, mapping = aes(x = reorder(changes, -value), y = value)) +
  geom_bar(stat = "identity", width = 0.5) + 
  theme_classic(8) +
  theme(axis.text.x = element_text(angle=60, vjust = 0.5))   + 
  theme(panel.border = element_rect( fill = NA)) +
  labs(x = "Type of changes", y = "Number of CRISPRs", title = "") 

acquire_2

pdf("../Figures/Fig6/b_acquire_2.pdf", width = 4.5/2.54, height = 6/2.54)
acquire_2
graphics.off()

###########################
########## Fig 7 ##########
###########################
# Fig 7a
across_data <- read.csv("../Files/across.csv", header = F, sep = ",")
colnames(across_data) <-  c("spacer_ID", "phylum", "genus")

across_phylum <- ggplot(data = across_data) + 
  geom_bar(mapping = aes(x = phylum), width = 0.5)  + 
  ylab("Number of spacers ") + xlab("Number of phylum") + 
  theme_classic(8) + theme(panel.border = element_rect( fill = NA))

pdf("../Figures/Fig7/a_spacer_phylum_across.pdf", width = 2.7/2.54, height = 5/2.54)
across_phylum
graphics.off()

# Fig 7b
across_genus <- ggplot(data = across_data) + 
  geom_bar(mapping = aes(x = genus), width = 0.5)  + 
  ylab("Number of spacers ") + xlab("Number of genera") + 
  theme_classic(8) + theme(panel.border = element_rect( fill = NA))

pdf("../Figures/Fig7/b_spacer_genus_across.pdf", width = 5.9/2.54, height = 5/2.54)
across_genus
graphics.off()  


###########################
####### statistics ########
###########################
# wilcox.test
# Fig 1
test_fig1 <- select(fig1_data, baby_id, gender, month, num_CRISPRs, spacers_sample)

month_12 <- subset(test_fig1, month == "12M")
boy <- subset(month_12, gender == "boy")
girl <- subset(month_12, gender == "girl")

boy_cris <- as.numeric(boy[["num_CRISPRs"]])
girl_cris <- as.numeric(girl[["num_CRISPRs"]])

boy_spacer <- as.numeric(boy[["spacers_sample"]])
girl_spacer <- as.numeric(girl[["spacers_sample"]])

wilcox.test(boy_cris, girl_cris, alternative = "greater", paired = F, exact = F)
wilcox.test(boy_spacer, girl_spacer, alternative = "greater", paired = F, exact = F)

# Fig 2
gender2com <- unique(select(CRIS_data, baby_id, gender))
colnames(richness_data)[1] <- "baby_id"
rich_12 <- subset(richness_data, Age == "12M")
rich2test <- left_join(rich_12, gender2com, by = "baby_id")


boy_rich <- subset(rich2test, gender == "boy")
girl_rich <- subset(rich2test, gender == "girl")

boy_rich <- as.numeric(boy_rich[["Bacterial_richness"]])
girl_rich <- as.numeric(girl_rich[["Bacterial_richness"]])

wilcox.test(boy_rich, girl_rich, alternative = "greater", paired = F, exact = F)

commu_data <- select(fig1_data, baby_id, gender, month, total_bp)
commu_12 <- subset(commu_data, month == "12M")

boy_commu <- subset(commu_12, gender == "boy")
girl_commu <- subset(commu_12, gender == "girl")

boy_commu <- as.numeric(boy_commu[["total_bp"]])
girl_commu <- as.numeric(girl_commu[["total_bp"]])

wilcox.test(boy_commu, girl_commu, alternative = "greater", paired = F, exact = F)

# linear model
CRIS_data$month <- gsub("B", "0",  CRIS_data$month)
CRIS_data$month <- gsub("4M", "4",  CRIS_data$month)
CRIS_data$month <- gsub("12M", "12",  CRIS_data$month)

richness_data$sample_id <- paste(richness_data$baby_id, richness_data$Age, sep = "_")

model_data <- left_join(CRIS_data, richness_data, by = "sample_id")

cris_model <- lm(num_CRISPRs ~ Bacterial_richness, data = model_data)
summary(cris_model)

spacer_model <- lm(spacers_sample ~ Bacterial_richness, data = model_data)
summary(spacer_model)

# anova test
cris_anov <- aov(num_CRISPRs ~ gender + month + delivery_mode + feeding_first_week + feeding_4m +
              breastfeeding_12m + antibiotic_0_4M + antibiotic_4_12m, data = CRIS_data)
summary(cris_anov)

spacer_anov <- aov(spacers_sample ~ gender + month + delivery_mode + feeding_first_week + feeding_4m +
                     breastfeeding_12m + antibiotic_0_4M + antibiotic_4_12m, data = CRIS_data)
summary(spacer_anov)
