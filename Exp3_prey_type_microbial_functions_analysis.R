#### Prey type influence microbial community structure and function ####
#### Authors: Jessica R. Bernardin, Sarah M. Gray, Leonora S. Bittleston ####
#### last update : 02/29/2024 ####
#### Prey Type ~ Microbial Functions Analysis

#### Load Required Packages ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "tidyverse", "effects", "growthcurver",
  "plyr", "dplyr", "reshape", "reshape2", "ape", "DiagrammeR", "tidyverse",
  "tidybayes", "coefplot", "standardize", "bayesplot", "MCMCvis", "car",
  "patchwork", "ggpubr", "corrr", "ggcorrplot", "factoextra", "MASS",
  "pairwiseAdonis", "plotrix", "gridExtra", "multcompView", "ggeffects", "this.path", "brms", "pracma",
  "phyloseq", "picante", "ape", "qiime2R", "decontam", "janitor", "readr", "ANCOMBC"
)

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

# Define hex codes for the color palette
ant_color <- "#c24000"  
beetle_color <- "#008388" 
fly_color <- "#5d004f"  
lm_color <- "#00446a"  
cb_color <- "#d88a00"  
pos_color <- "#006a50"  
water_color <- "slategray" 

# Combine into a vector for easy reference
color_palette <- c(ant_color, beetle_color, fly_color, lm_color, cb_color, pos_color, water_color)
color_palette_3 <- c(ant_color, beetle_color, fly_color,  cb_color)
color_palette_4 <- c(lm_color,  cb_color)

#### Insect Nutrient Compositions ####
prey_data <- read.csv("data/Insect_Nutrient_metadata.csv", header=TRUE)

## Fig 1, lipid
ggplot(data=prey_data, aes(x=prey, y=prop_bodymass_is_lipid, fill=prey, color=prey)) + 
  geom_boxplot(outlier.shape=NA,alpha=.6)+geom_jitter(size=3, color="black")+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=15))+
  scale_fill_manual(values=c("#c24000", "#008388" ,"#5d004f"))+
  scale_color_manual(values=c("#c24000", "#008388" ,"#5d004f"))+theme(legend.position = "none")+
  ylab("Prop Lipid")

## Fig 1, exoskeleton
ggplot(data=prey_data, aes(x=prey, y=prop_body_mass_is_exoskeleton, fill=prey, color=prey)) + 
  geom_boxplot(outlier.shape=NA,alpha=.6)+geom_jitter(size=3, color="black")+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=15))+
  scale_fill_manual(values=c("#c24000", "#008388" ,"#5d004f"))+
  scale_color_manual(values=c("#c24000", "#008388" ,"#5d004f"))+theme(legend.position = "none")+
  ylab("Prop Exo")

## Fig 1, protein
ggplot(data=prey_data, aes(x=prey, y=prop_protein, fill=prey, color=prey)) + 
  geom_boxplot(outlier.shape=NA,alpha=.6)+geom_jitter(size=3, color="black")+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=15))+
  scale_fill_manual(values=c("#c24000", "#008388" ,"#5d004f"))+
  scale_color_manual(values=c("#c24000", "#008388" ,"#5d004f"))+theme(legend.position = "none")+
  ylab("Proportion")

#### Nutrient GLMs
prey_data$prey <- as.factor(prey_data$prey)
prey_data$prey <- relevel(prey_data$prey, ref = "fly")

## Exoskeleton
#m1_exo <- brm(prop_body_mass_is_exoskeleton ~ prey, 
#             family = Beta(), 
#             data = prey_data, iter=10000)
#saveRDS(m1_exo, "RDS_files/m1_exoskeleton.RDS")
m1_exo <- readRDS("RDS_files/m1_exoskeleton.RDS")
summary(m1_exo)
exo_pred <- ggpredict(m1_exo)
plot(exo_pred)
pp_check(m1_exo)

mcmc_areas(m1_exo, regex_pars="b_") + geom_vline(xintercept = 0, linetype=3)+ggtitle("%exoskeleton~preytype")
posterior_m1_exo<- as.data.frame(m1_exo)
posterior_m1_exo_melt <- posterior_m1_exo[,1:3]
posterior_m1_exo_melt <- reshape2::melt(posterior_m1_exo_melt)
level_order <- c("b_preybeetle", "b_preyant","b_Intercept") 

##Fig.S1
ggplot(posterior_m1_exo_melt, aes(x = value,y=factor(variable, level=level_order), fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") + xlab("Estimate effect on proportion exoskeleton")+
  theme_classic() + guides(fill="none")+ 
  scale_fill_manual(values=c("#5d004f","#c24000","#008388" ))

## Lidid
#m1_lipid <- brm(prop_bodymass_is_lipid ~ prey, 
#              family = Beta(), 
#              data = prey_data)
#saveRDS(m1_lipid, "RDS_files/m1_lipid.RDS")
m1_lipid <- readRDS("RDS_files/m1_lipid.RDS")
lip_pred <- ggpredict(m1_lipid)
pp_check(m1_lipid)
summary(m1_lipid)
mcmc_areas(m1_lipid, regex_pars="b_")+ geom_vline(xintercept = 0, linetype=3)+ggtitle("%lipid~preytype")
posterior_m2_lipid<- as.data.frame(m1_lipid)
posterior_m2_lipid_melt <- posterior_m2_lipid[,1:3]
posterior_m2_lipid_melt <- reshape2::melt(posterior_m2_lipid_melt)

##Fig.S1
ggplot(posterior_m2_lipid_melt, aes(x = value,y=factor(variable, level=level_order), fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") + xlab("Estimate effect on proportion lipid")+
  theme_classic() + guides(fill="none")+ 
  scale_fill_manual(values=c("#5d004f","#c24000","#008388"))

##Protein
#m3_low <- brm(prop_protein ~ prey, 
#                family = Beta(), 
#                data = prey_data, iter=10000)

#saveRDS(m3_low, "RDS_files/m1_protein.RDS")
m3_low <- readRDS("RDS_files/m1_protein.RDS")
mcmc_areas(m3_low, regex_pars="b_")+ geom_vline(xintercept = 0, linetype=3)+ggtitle("Prop Protein~preytype")
pp_check(m3_low)
low_pred <- ggpredict(m3_low)
plot(low_pred)

summary(m3_low)
posterior_m3_low<- as.data.frame(m3_low)
posterior_m3_low_melt <- posterior_m3_low[,1:3]
posterior_m3_low_melt <- reshape2::melt(posterior_m3_low_melt)

##Fig.S1
ggplot(posterior_m3_low_melt, aes(x = value,y=factor(variable, level=level_order), fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") + xlab("Estimate effect on proportion protein")+
  theme_classic() + guides(fill="none")+ 
  scale_fill_manual(values=c("#5d004f","#c24000","#008388" ))

#### Microbial Physiological Functions
#bacterial and plant metadata
data <- read.csv("data/Prey_Exp4_meta.csv", header = TRUE)

#filter out plants that had damage
data_filt <- data[!(data$plant %in% c("8", "24", "NEG")), ]
data_filt$treatment <- as.factor(data_filt$treatment)
data_filt$plant <- as.factor(data_filt$plant)
data_filt$day <- as.numeric(data_filt$day)
data_filt$chitinase_uM_min <- as.numeric(data_filt$chitinase_uM_min)
data_filt$protease_nM_min <- as.numeric(data_filt$protease_nM_min)
data_filt$flow_cyt_livingevents_uL <- as.numeric(data_filt$flow_cyt_livingevents_uL)
data_filt$flow_cyt_totalevents_uL <- as.numeric(data_filt$flow_cyt_totalevents_uL)
data_filt <- data_filt %>% mutate(prop_living = flow_cyt_livingevents_uL/flow_cyt_totalevents_uL)

prey_order <- c("ant", "beetle", "fly", "lm", "cb", "pos", "water")
prey_order2 <- c("ant", "beetle", "fly", "lm", "cb")
prey_order3 <- c("ant", "beetle", "fly", "lm", "cb", "water")

# Use the factor function to set the order of levels
data_filt$treatment <- factor(data_filt$treatment, levels = prey_order)

#### Chitinase Activity ####
ggplot(data=data_filt, aes(x=day, y=chitinase_uM_min, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Chitinase Activity") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")


#### Protease Activity ####
ggplot(data=data_filt, aes(x=day, y=protease_nM_min, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Protease Activity") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

#### Bacterial Abundance ####
data_filt$dna_con_ug_ml <- as.numeric(data_filt$dna_con_ug_ml)
ggplot(data=data_filt, aes(x=day, y=dna_con_ug_ml, color = treatment, group=treatment))+ 
  geom_jitter() +facet_wrap(~treatment, nrow=1) +
  ylab("DNA concentration (ng/mL)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

data_filt %>% filter(!day %in% c(0, 21, 28, 42, 49)) %>% 
  ggplot(aes(x=treatment, y=dna_con_ug_ml, color = treatment, group=treatment))+ 
  geom_boxplot(outlier.shape=NA) +geom_jitter() +facet_wrap(~day, nrow=1) +
  ylab("DNA concentration (ng/mL)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

ggplot(data=data_filt, aes(x=day, y=flow_cyt_livingevents_uL, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Total number of living cells (per uL)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

#### pH ####
data_filt_ph <- data_filt[complete.cases(data_filt$pH), ]

ggplot(data=data_filt_ph, aes(x=day, y=pH, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE) +geom_jitter(alpha=.5)+ facet_wrap(~treatment, nrow=1) +
  ylab("Pitcher Fluid pH") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

#### Decomposition ####
decom <- data_filt
decom <- decom[complete.cases(decom$prey_mass_g), ]
decom$prey_mass_g <- as.numeric(decom$prey_mass_g)
decom$treatment <- factor(decom$treatment, levels = prey_order)
decom_delta <- data_filt[complete.cases(data_filt$delta_prey_percent), ]
decom_delta$treatment <- factor(decom_delta$treatment, levels = prey_order)

ggplot(data=decom_delta, aes(x=treatment, y=delta_prey_percent , fill=treatment, group=treatment, color=treatment))+ 
  geom_boxplot(outlier.shape=NA, alpha=.7)+ geom_jitter(size=4)+
  ylab("Proportion Mass Lost") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(text = element_text(size=20))+scale_fill_manual(values=c(color_palette))+
  theme(legend.position = "none")

#### ECOPLATES ####
ecodat <- read.csv("data/Exp4_Ecoplate_all.csv", header = TRUE)
eco_filt <- filter(ecodat, treatment %in% c("ant", "beetle", "fly"))
prey_order2 <- c("ant", "beetle", "fly")

# Use the factor function to set the order of levels
eco_filt$treatment <- factor(eco_filt$treatment, levels = prey_order2)
com <- eco_filt
com_meta <- subset(com, select= c(sample, plant_block, treatment, week))
com <- subset(com, select= -c(sample, plant_block, treatment, week))
com[com < 0] <- 0

#run NMDS on ecoplate data, no controls
set.seed(123)
nmds <- metaMDS(com, k=2, trymax=500) #stress=0.13948

#get the data scores from the nmds
all.scores <- as.data.frame(scores(nmds, "sites"))

#add metadata to the scores info
all.scores <- cbind(all.scores, com_meta)
all.scores$week <- as.factor(all.scores$week)
# make a new column that has unique names for week and treatment
all.scores <- all.scores %>% 
  unite(hull.id, treatment,remove = FALSE)

hullsall <- all.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

ggplot(all.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment, shape = week)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "week") + aes(fill = factor(treatment)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+
  scale_color_manual(values=color_palette)+scale_fill_manual(values=color_palette)

ecoall.bd <- betadisper(vegdist(com), com_meta$treatment)
ecoall.bd
boxplot(ecoall.bd)
anova(ecoall.bd)

ecoall.bdweek <- betadisper(vegdist(com), com_meta$week)
ecoall.bdweek
boxplot(ecoall.bdweek)
anova(ecoall.bdweek)

# PERMANOVA for categorical variables (factors)
# set number of permutations
perm <- how(nperm = 999)
adonis <- adonis2(com ~ treatment + week, data=all.scores, permutations = perm, by  = "margin")
adonis
                 

# pairwise adonis
#Table S2
ecoall.ad.pw <- pairwise.adonis2(com ~ treatment, data = com_meta, strata = 'week')
ecoall.ad.pw

ecoall.ad.pw2 <- pairwise.adonis2(com ~ week, data = com_meta)
ecoall.ad.pw2

#### local vs non-local prey EcoPlates ####
eco_filt2 <- filter(ecodat, treatment %in% c("lm", "cb"))
prey_order3 <- c("lm", "cb")

# Use the factor function to set the order of levels
eco_filt2$treatment <- factor(eco_filt2$treatment, levels = prey_order3)
com2 <- eco_filt2
com2_meta <- subset(com2, select= c(sample, plant_block, treatment, week))
com2 <- subset(com2, select= -c(sample, plant_block, treatment, week))
com2[com2 < 0] <- 0

#run NMDS on ecoplate data, no controls
set.seed(123)
nmds2 <- metaMDS(com2, k=2, trymax=500) #stress=0.09028237

#get the data scores from the nmds
all.scores2 <- as.data.frame(scores(nmds2, "sites"))
#add metadata to the scores info

all.scores2 <- cbind(all.scores2, com2_meta)
all.scores2$week <- as.factor(all.scores2$week)
# make a new column that has unique names for week and treatment
all.scores2 <- all.scores2 %>% 
  unite(hull.id, treatment,remove = FALSE)

hullsall <- all.scores2 %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

ggplot(all.scores2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment, shape = week)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "week") + aes(fill = factor(treatment)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+
  scale_color_manual(values=c("#00446a", "#d88a00"))+scale_fill_manual(values=c("#00446a", "#d88a00"))

# beta disper to check dispersion
ecoall.bd <- betadisper(vegdist(com2), com2_meta$treatment)
ecoall.bd
boxplot(ecoall.bd)
anova(ecoall.bd)

ecoall.bdweek <- betadisper(vegdist(com2), com2_meta$week)
ecoall.bdweek
boxplot(ecoall.bdweek)
anova(ecoall.bdweek)

# PERMANOVA for categorical variables (factors)
# set number of permutations
perm <- how(nperm = 999)
adonis <- adonis2(com2 ~ treatment + week, data=all.scores2, permutations = perm, by  = "margin")
adonis

#### EcoPlate Substrate Analysis ####
eco_melt <- melt(eco_filt, id.vars = c("sample", "plant_block", "treatment", "week"))

eco_melt$treatment <- relevel(eco_melt$treatment, ref = "fly")

eco_melt_all <- eco_melt %>%
  mutate(category = case_when(
    grepl("L.Arginine", variable) | grepl("L.Asparagine", variable) | grepl("L.Phenylalanine", variable) | grepl("L.Serine", variable) | grepl("Glycyl.L.glutamic.acid", variable)| grepl("L.Threonine", variable) ~ "amino_acid",
    grepl("Phenylethylamine", variable) | grepl("Putrescine", variable)  ~ "amine",
    grepl("D.Mannitol", variable) | grepl("Glucose.1.phosphate", variable)  | grepl("Beta.Methyl.D.Glucoside", variable) | grepl("D.L.alpha.glycerol.phosphate", variable) | grepl("D.Galactonic.acid.gamma.Lactone", variable) | grepl("I.Erythritol", variable) | grepl("D.Xylose", variable) | grepl("N.Acetyl.D.Glucosamine", variable) | grepl("D.Cellobiose", variable) | grepl("alpha.D.Lactose", variable) ~ "carbohydrate",
    grepl("D.Glucosaminic.acid", variable) | grepl("D.Malic.acid", variable) | grepl("Itaconic.acid", variable) | grepl("Pyruvic.acid.methyl.ester", variable) | grepl("D.Galacturonic.acid", variable) | grepl("alpha.Ketobutyric.acid", variable) | grepl("gamma.Hydroxybutyric.acid", variable) ~ "carboxylic_acid",
    grepl("X4.Hydroxy.benzoic.acid", variable) | grepl("X2.Hydroxy.benzoic.acid", variable)  ~ "phenolic_compound",
    grepl("Glycogen", variable) | grepl("Alpha.cyclodextrin", variable) |  grepl("Tween.40", variable) | grepl("Tween.80", variable) ~ "polymer",
    FALSE ~ "other_category"
  ))
#eco_m2 <-brm(value~treatment*variable, data=eco_melt, chains=4, iter=5000)
#saveRDS(eco_m2, "RDS_files/Exp4_ecosubstrates.RDS")
eco_m2 <- readRDS("RDS_files/Exp4_ecosubstrates.RDS")

pp_check(eco_m2)
pred1 <- ggpredict(eco_m2, terms = c("variable", "treatment"))
prey_order4 <- c("ant", "beetle", "fly")
pred1$group <- factor(pred1$group, levels = prey_order4)

colnames(pred1)[1] <- "variable"
category <- eco_melt_all[,c(5,7)]
category_unique <- category %>% distinct(variable, .keep_all = TRUE)
merged_df <- merge(pred1, category_unique, by = "variable", all.x = TRUE)

merged_df$variable <- factor(merged_df$variable, levels = unique(merged_df[order(merged_df$category), "variable"]))
prey_order3 <- c("ant", "beetle", "fly")

ggplot(merged_df, aes(x = variable, y = predicted, ymin = conf.low, ymax = conf.high, color = factor(group, level=prey_order3))) +
  geom_pointrange(position = position_dodge(width = 0.7), alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values=c("#c24000","#008388" ,"#5d004f" ))

posterioreco_m1 <- mcmc_intervals_data(eco_m2, 
                                       prob_outer=0.95,
                                       prob=0.5)

posterioreco_m1$nonzero <- NA
posterioreco_m1$nonzero[posterioreco_m1$ll>0 & posterioreco_m1$hh>0] <- "nonzero"
posterioreco_m1$nonzero[posterioreco_m1$ll<0 & posterioreco_m1$hh<0] <- "nonzero"
posterioreco_m1$nonzero[is.na(posterioreco_m1$nonzero)] <- "zero"
posterioreco_m1<- posterioreco_m1[34:93,]

ggplot(posterioreco_m1, aes(x = parameter,
                            shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="", values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), labels=c("95% CI does\nnot contain zero",  "95% CI\ncontains zero"))+
  coord_flip() +theme_classic() + 
  theme(axis.text.y = element_text( size=7), axis.text.x=element_text(size=7), axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) + xlab(NULL) + guides(linetype=FALSE) 


#FILTER TO THOSE THAT DON'T CROSS 0
merged_df_filt <- filter(merged_df, variable %in% c( "alpha.D.Lactose",
                                                     "Glycyl.L.glutamic.acid", "L.Threonine",
                                                     "D.Cellobiose"))

ggplot(merged_df_filt, aes(x = variable, y = predicted, ymin = conf.low, ymax = conf.high, color = factor(group, level=prey_order3))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values=c("#c24000","#008388" ,"#5d004f" ))+theme(legend.position = "none")

#### GLMMs ####
#relevel the predictor variables
data_filt2 <- filter(data_filt, treatment %in% c("ant", "beetle", "fly"))
data_filt2$pH <- as.numeric(data_filt2$pH)
data_filt2 <- data_filt2 %>% arrange(desc(treatment))
data_filt2$treatment <- relevel(data_filt2$treatment, ref = "fly")
decom_delta2 <- filter(decom_delta, treatment %in% c("ant", "beetle", "fly"))
decom_delta2$treatment <- relevel(decom_delta2$treatment, ref = "fly")

#chitinase across time
data_filt2$chitinase_uM_min <- ifelse(data_filt2$chitinase_uM_min <= 0, 0.0001, data_filt2$chitinase_uM_min)
m1a <- brm(chitinase_uM_min ~ treatment+day+(1|plant_block),data=data_filt2, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))


#saveRDS(m1a, "RDS_files/Exp4_chitinase_glm.RDS")
m1a <- readRDS("RDS_files/Exp4_chitinase_glm.RDS")
summary(m1a)
pp_check(m1a)

posterior_m1a <- as.data.frame(m1a)
posterior_m1a_melt <- posterior_m1a[,1:3]
posterior_m1a_melt <- reshape2::melt(posterior_m1a_melt)
ggplot(posterior_m1a_melt, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f","#c24000","#008388"))+
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect on Chitinase Acitivity")
##Fig3
level_order <- c("b_treatmentfly","b_treatmentbeetle","b_treatmentant","b_Intercept") 

posterior_m1a_melt_intc <- filter(posterior_m1a_melt, variable =="b_Intercept")
posterior_m1a_melt_other <- filter(posterior_m1a_melt, variable %in% c("b_treatmentant", "b_treatmentbeetle"))

ggplot(posterior_m1a_melt_intc, aes(x = value,y=factor(variable, level=level_order), fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f" ))

ggplot(posterior_m1a_melt_other, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#c24000","#008388"))

pm1aa <- ggeffects::ggpredict(m1a, "treatment")
pm1aaa <- ggeffects::ggpredict(m1a, "day")

pm1a <- ggeffects::ggpredict(m1a, c("day", "treatment"))
colnames(pm1a)[1] <- "day"
colnames(pm1a)[2] <- "chitinase_uM_min"
colnames(pm1a)[5] <- "treatment"

level_order <- c("ant","beetle","fly") 

ggplot() +
  geom_jitter(data = data_filt2, aes(x = day, y = chitinase_uM_min, color = factor(treatment, level=level_order), fill = factor(treatment, level=level_order))) +
  geom_line(data = pm1a, aes(x = day, y = chitinase_uM_min, color = treatment, group = treatment), lwd=1.5) +
  geom_ribbon(data = pm1a, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.1) + # Shaded area for confidence interval
  scale_color_manual(values = c("#c24000", "#008388","#5d004f")) +
  scale_fill_manual(values = c("#c24000", "#008388","#5d004f"))+theme_classic()+
 theme(legend.position ="none")+ xlim(c(7,52))+theme(text = element_text(color = "black", size = 16 ))+
  xlab("Day")+ ylab("Chitinase Activity")

##protease
data_filt2$protease_nM_min <- ifelse(data_filt2$protease_nM_min <= 0, 0.0001, data_filt2$protease_nM_min)

m2b <- brm(protease_nM_min ~ treatment+day+(1|plant_block),data=data_filt2, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4,
         control = list(adapt_delta = 0.999, max_treedepth = 20))
#saveRDS(m2a, "RDS_files/Exp4_protease_glm.RDS")
m2a <- readRDS("RDS_files/Exp4_protease_glm.RDS")
summary(m2a)
summary(m2b)
pp_check(m2a)

posterior_m2a <- as.data.frame(m2a)
posterior_m2a_melt <- posterior_m2a[,1:3]
posterior_m2a_melt <- reshape2::melt(posterior_m2a_melt)
ggplot(posterior_m2a_melt, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f","#c24000","#008388"))+
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect on Protease Acitivity")
##Fig3
posterior_m2a_melt_intc <- filter(posterior_m2a_melt, variable =="b_Intercept")
posterior_m2a_melt_other <- filter(posterior_m2a_melt, variable %in% c("b_treatmentant", "b_treatmentbeetle"))

ggplot(posterior_m2a_melt_intc, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f"  ))

ggplot(posterior_m2a_melt_other, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#c24000","#008388"))

pm2aa <- ggeffects::ggpredict(m2a, "treatment")

pm2a <- ggeffects::ggpredict(m2a, c("day", "treatment"))
colnames(pm2a)[1] <- "day"
colnames(pm2a)[2] <- "protease_nM_min"
colnames(pm2a)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt2, aes(x = day, y = protease_nM_min, color = factor(treatment, level=level_order), fill = factor(treatment, level=level_order))) +
  geom_line(data = pm2a, aes(x = day, y = protease_nM_min, color = treatment, group = treatment), lwd=1.5) +
  geom_ribbon(data = pm2a, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.1) + # Shaded area for confidence interval
  scale_color_manual(values = c("#c24000", "#008388","#5d004f")) +
  scale_fill_manual(values = c("#c24000", "#008388","#5d004f"))+theme_classic()+
  theme(legend.position ="none") + xlim(c(7,52))+theme(text = element_text(color = "black", size = 16 ))+
  xlab("Day")+ ylab("Protease Activity")


#m3a <- brm(flow_cyt_livingevents_uL ~ treatment+day+(1|plant_block),data=data_filt2, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4,
#          control = list(adapt_delta = 0.999, max_treedepth = 20))
#saveRDS(m3a, "RDS_files/Exp4_cyto_glm.RDS")

summary(m3a)
pp_check(m3a)

posterior_m3a <- as.data.frame(m3a)
posterior_m3a_melt <- posterior_m3a[,1:4]
posterior_m3a_melt <- reshape2::melt(posterior_m3a_melt)

posterior_m3a_melt_intc <- filter(posterior_m3a_melt, variable =="b_Intercept")
posterior_m3a_melt_other <- filter(posterior_m3a_melt, variable %in% c("b_treatmentant", "b_treatmentbeetle"))

ggplot(posterior_m3a_melt_intc, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f" ))

ggplot(posterior_m3a_melt_other, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#c24000","#008388"))

pm3aa <- ggeffects::ggpredict(m3a,  "treatment")

pm3a <- ggeffects::ggpredict(m3a, c("day", "treatment"))
colnames(pm3a)[1] <- "day"
colnames(pm3a)[2] <- "flow_cyt_livingevents_uL"
colnames(pm3a)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt2, aes(x = day, y = flow_cyt_livingevents_uL, color = factor(treatment, level=level_order), fill = factor(treatment, level=level_order))) +
  geom_line(data = pm3a, aes(x = day, y = flow_cyt_livingevents_uL, color = treatment, group = treatment), lwd=1.5) +
  geom_ribbon(data = pm3a, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.1) + # Shaded area for confidence interval
  scale_color_manual(values = c("#c24000", "#008388","#5d004f")) +
  scale_fill_manual(values = c("#c24000", "#008388","#5d004f"))+theme_classic()+
  theme(legend.position ="none")+ xlim(c(7,52))


##ph
#m4a <- brm(pH ~ treatment+day+(1|plant_block), data=data_filt2, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)
#saveRDS(m4a, "RDS_files/Exp4_pH_glm.RDS")
m4a <-readRDS("RDS_files/Exp4_pH_glm.RDS")
summary(m4a)
pp_check(m4a)
posterior_m4a <- as.data.frame(m4a)
posterior_m4a_melt <- posterior_m4a[,1:3]
posterior_m4a_melt <- reshape2::melt(posterior_m4a_melt)

ggplot(posterior_m4a_melt, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f","#c24000","#008388"))+
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect on pH")


posterior_m4a_melt_intc <- filter(posterior_m4a_melt, variable =="b_Intercept")
posterior_m4a_melt_other <- filter(posterior_m4a_melt, variable %in% c("b_treatmentant", "b_treatmentbeetle"))

ggplot(posterior_m4a_melt_intc, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f" ))

ggplot(posterior_m4a_melt_other, aes(x = value,y=variable, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#c24000","#008388"))

pm4aa <- ggeffects::ggpredict(m4a, "treatment")

pm4a <- ggeffects::ggpredict(m4a, c("day", "treatment"))
colnames(pm4a)[1] <- "day"
colnames(pm4a)[2] <- "pH"
colnames(pm4a)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt2, aes(x = day, y = pH, color = factor(treatment, level=level_order), fill = factor(treatment, level=level_order))) +
  geom_line(data = pm4a, aes(x = day, y = pH, color = treatment, group = treatment), lwd=1.5) +
  geom_ribbon(data = pm4a, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.1) + # Shaded area for confidence interval
  scale_color_manual(values = c("#c24000", "#008388","#5d004f")) +
  scale_fill_manual(values = c("#c24000", "#008388","#5d004f"))+theme_classic()+
  theme(legend.position ="none")+ xlim(c(7,52))+theme(text = element_text(color = "black", size = 16 ))+
  xlab("Day") + ylab("pH")

##change in biomass

#m5a <- brm(delta_prey_percent ~ treatment, data=decom_delta2, family=Beta(), iter = 10000, chains = 4, cores = 4)
#m5a
#saveRDS(m5a, "RDS_files/Exp4_decomp_glm.RDS")

posterior_m5a <- as.data.frame(m5a)
posterior_m5a_melt <- posterior_m5a[,1:3]
posterior_m5a_melt <- reshape2::melt(posterior_m5a_melt)

posterior_m5a_melt_intc <- filter(posterior_m5a_melt, variable =="b_Intercept")
posterior_m5a_melt_other <- filter(posterior_m5a_melt, variable !="b_Intercept")

#Figure 2B
ggplot(posterior_m5a_melt_intc, aes(x = value, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#5d004f" ))

#Figure 2C
ggplot(posterior_m5a_melt_other, aes(x = value, fill = factor(variable))) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ scale_fill_manual(values=c( "#c24000","#008388" ))

pm5a <- ggeffects::ggpredict(m5a, "treatment")
colnames(pm5a)[1] <- "treatment"
colnames(pm5a)[2] <- "delta_prey_percent"
level_order <- c("ant","beetle","fly") 

#Figure 2A
ggplot(decom_delta2, 
       mapping = aes(x = factor(treatment, level=level_order), 
                     y = delta_prey_percent, 
                     fill = treatment, color=treatment))  +
  geom_violin(trim=FALSE,alpha=.6) + geom_jitter( size=4)+
  scale_color_manual(values=c("#5d004f","#c24000","#008388" )) +
  scale_fill_manual(values=c( "#5d004f","#c24000","#008388" )) +
  geom_point(data=pm5a, aes(x=treatment, y=delta_prey_percent), color="black", size=5) +
  geom_linerange(data=pm5a,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Proportion Change Prey Biomass") + theme_classic()+
  theme(legend.position="none")

#### LM vs CB GLMS ####
##filter data and relevel
data_filt3 <- filter(data_filt, treatment %in% c("lm", "cb"))
data_filt3$pH <- as.numeric(data_filt3$pH)
data_filt3 <- data_filt3 %>% arrange(desc(treatment))
data_filt3$treatment <- relevel(data_filt3$treatment, ref = "lm")

decom_delta3 <- filter(decom_delta, treatment %in% c("lm", "cb"))
decom_delta3$treatment <- relevel(decom_delta3$treatment, ref = "lm")

##chitinase across time
data_filt3$chitinase_uM_min <- ifelse(data_filt3$chitinase_uM_min <= 0, 0.0001, data_filt3$chitinase_uM_min)
#m1b <- brm(chitinase_uM_min ~ treatment+day,data=data_filt3,family = Gamma(link="log"),iter = 10000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
#saveRDS(m1b, file = "RDS_files/brms_lmcb_chit.RDS")

summary(m1b)
pp_check(m1b)

posterior_m1b <- as.data.frame(m1b)
posterior_m1b_melt <- posterior_m1b[,1:3]
posterior_m1b_melt <- reshape2::melt(posterior_m1b_melt)
level_order <- c("b_treatmentcb","b_Intercept") 

ggplot(posterior_m1b_melt, aes(x = value,y=variable, fill = variable)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+
  scale_fill_manual(values = c("#00446a", "#d88a00", "black"))


pm1b <- ggeffects::ggpredict(m1b, c("day", "treatment"))
colnames(pm1b)[1] <- "day"
colnames(pm1b)[2] <- "chitinase_uM_min"
colnames(pm1b)[5] <- "treatment"

level_order <- c("lm","cb") 

ggplot() +
  geom_jitter(data = data_filt3, aes(x = day, y = chitinase_uM_min, color = factor(treatment, level=level_order), fill = factor(treatment, level=level_order))) +
  geom_line(data = pm1b, aes(x = day, y = chitinase_uM_min, color = treatment, group = treatment)) +
  geom_ribbon(data = pm1b, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.2) + # Shaded area for confidence interval
  scale_color_manual(values = c("#00446a", "#d88a00")) +
  scale_fill_manual(values = c("#00446a", "#d88a00"))+theme_classic()+
  facet_wrap(~factor(treatment, level=level_order), nrow=1)+theme(legend.position ="none")


#protease
#m2b <- brm(protease_nM_min ~ treatment+day,data=data_filt3,family = Gamma(link="log"), iter = 10000, chains = 4, cores = 4,
#           control = list(adapt_delta = 0.999, max_treedepth = 20))
#saveRDS(m2b, file = "RDS_files/brms_lmcb_prot.RDS")

summary(m2b)
pp_check(m2b)
conditional_effects(m2b)

posterior_m2b<- as.data.frame(m2b)
posterior_m2b_melt <- posterior_m2b[,1:3]
posterior_m2b_melt <- reshape2::melt(posterior_m2b_melt)

ggplot(posterior_m2b_melt, aes(x = value,y=variable, fill = variable)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+
  scale_fill_manual(values = c("#00446a", "#d88a00", "black"))

pm2b <- ggeffects::ggpredict(m2b, c("day", "treatment"))
colnames(pm2b)[1] <- "day"
colnames(pm2b)[2] <- "protease_nM_min"
colnames(pm2b)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt3, aes(x = day, y = protease_nM_min, color = treatment, fill = treatment)) +
  geom_line(data = pm2b, aes(x = day, y = protease_nM_min, color = treatment)) +
  geom_ribbon(data = pm2b, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.2) + # Shaded area for confidence interval
  scale_color_manual(values = c("#00446a", "#d88a00")) +
  scale_fill_manual(values = c("#00446a", "#d88a00"))+theme_classic()+
  facet_wrap(~treatment, nrow=1)+theme(legend.position ="none")


##bacterial abundance (no difference between block so removed)
#m3b <- brm(flow_cyt_livingevents_uL ~ treatment+day,data=data_filt3,family = Gamma(link="log"), iter = 10000, chains = 4, cores = 4,
#           control = list(adapt_delta = 0.999, max_treedepth = 20))
#saveRDS(m3b, file = "RDS_files/brms_lmcb_propliving.RDS")

summary(m3b)
m3b

pp_check(m3b)

posterior_m3b <- as.data.frame(m3b)
posterior_m3b_melt <- posterior_m3b[,1:3]
posterior_m3b_melt <- reshape2::melt(posterior_m3b_melt)

ggplot(posterior_m3b_melt, aes(x = value,y=variable, fill = variable)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+
  scale_fill_manual(values = c("#00446a", "#d88a00", "black"))


pm3b <- ggeffects::ggpredict(m3b, c("day", "treatment"))
colnames(pm3b)[1] <- "day"
colnames(pm3b)[2] <- "flow_cyt_livingevents_uL"
colnames(pm3b)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt3, aes(x = day, y = flow_cyt_livingevents_uL, color = treatment, fill = treatment)) +
  geom_line(data = pm3b, aes(x = day, y = flow_cyt_livingevents_uL, color = treatment)) +
  geom_ribbon(data = pm3b, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.2) + # Shaded area for confidence interval
  scale_color_manual(values = c("#00446a", "#d88a00")) +
  scale_fill_manual(values = c("#00446a", "#d88a00"))+theme_classic()+
  facet_wrap(~treatment, nrow=1)+theme(legend.position ="none")

#ph
#m4b <- brm(pH ~ treatment+day, data=data_filt3, iter = 10000,family = Gamma(link="log"), chains = 4, cores = 4)
#saveRDS(m4b, file = "RDS_files/brms_lmcb_ph.RDS")

summary(m4b)
pp_check(m4b)

posterior_m4b <- as.data.frame(m4b)
posterior_m4b_melt <- posterior_m4b[,1:3]
posterior_m4b_melt <- reshape2::melt(posterior_m4b_melt)

ggplot(posterior_m4b_melt, aes(x = value,y=variable, fill = variable)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+
  scale_fill_manual(values = c("#00446a", "#d88a00", "black"))

pm4b <- ggeffects::ggpredict(m4b, c("day", "treatment"))
colnames(pm4b)[1] <- "day"
colnames(pm4b)[2] <- "pH"
colnames(pm4b)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt3, aes(x = day, y= pH, color = treatment, fill = treatment)) +
  geom_line(data = pm4b, aes(x = day, y = pH, color = treatment)) +
  geom_ribbon(data = pm4b, aes(x = day, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.2) + # Shaded area for confidence interval
  scale_color_manual(values = c("#00446a", "#d88a00")) +
  scale_fill_manual(values = c("#00446a", "#d88a00"))+theme_classic()+
  facet_wrap(~treatment, nrow=1)+theme(legend.position ="none")

#change in biomass

#m5b <- brm(delta_prey_percent ~ treatment, data=decom_delta3, family = Beta(), iter = 10000, chains = 4, cores = 4)
#saveRDS(m5b, file = "RDS_files/brms_lmcb_biomass.RDS")
pp_check(m5b)

posterior_m5b <- as.data.frame(m5b)
posterior_m5b_melt <- posterior_m5b[,1:2]
posterior_m5b_melt <- reshape2::melt(posterior_m5b_melt)

ggplot(posterior_m5b_melt, aes(x = value,y=variable, fill = variable)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+
  scale_fill_manual(values = c("#00446a", "#d88a00", "black"))

pm5b <- ggeffects::ggpredict(m5b, "treatment")
colnames(pm5b)[1] <- "treatment"
colnames(pm5b)[2] <- "delta_prey_percent"
level_order <- c("lm","cb") 
ggplot(decom_delta3, 
       mapping = aes(x = factor(treatment, level=level_order), 
                     y = delta_prey_percent, 
                     fill = treatment, color=treatment))  +
  geom_violin(trim=FALSE,alpha=.6) + geom_jitter( size=4)+
  scale_color_manual(values=c("#00446a", "#d88a00")) +
  scale_fill_manual(values=c( "#00446a", "#d88a00" )) +
  geom_point(data=pm5b, aes(x=treatment, y=delta_prey_percent), color="black", size=5) +
  geom_linerange(data=pm5b,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Proportion Change Prey Biomass") + theme_classic()+
  theme(legend.position="none")

posteriorm5b <- mcmc_intervals_data(m5b, 
                                     prob_outer=0.95,
                                     prob=0.95)
posteriorm5b$nonzero <- NA
posteriorm5b$nonzero[posteriorm5b$ll>0 & posteriorm5b$hh>0] <- "nonzero"
posteriorm5b$nonzero[posteriorm5b$ll<0 & posteriorm5b$hh<0] <- "nonzero"
posteriorm5b$nonzero[is.na(posteriorm5b$nonzero)] <- "zero"
posteriorm5b<- posteriorm5b[1:2,]

#FIG6
ggplot(posteriorm5b, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_classic() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  guides(linetype=FALSE) 

#### metabarcoding ####
#read in 16S
asv16s <- read_tsv("metabarcoding/exported-files/asv-table-dada2.txt")
asv16s <- asv16s %>% column_to_rownames(var="#OTU ID")
asv16s <- asv16s[order(rownames(asv16s)), ]

meta <- data
meta <- meta %>% remove_rownames %>% column_to_rownames(var="metabarcoding_id")

tax.16s <- read_tsv("metabarcoding/exported-files/taxonomy.tsv")
tax.16s <- tax.16s %>% column_to_rownames(var="Feature ID")
tax.16s <- tax.16s[order(rownames(tax.16s)), ]

row.names(asv16s) == row.names(tax.16s) # sanity check
nrow(asv16s) #=3342 ASVs
nrow(tax.16s) #=3342 ASVs

asv16s.tree <- read.tree("metabarcoding/exported-files/tree.nwk")
asv16s.tree <-root(asv16s.tree, "63ffc596781727668a90d6da135a33b4")## root with an archaeon

#remove samples from experiment 2
asv16s.exp4<- asv16s %>%
  dplyr::select(-matches("exp2"))
#3342 asvs, 102 samples including neg controls

asv16s.exp4 <- asv16s.exp4 %>%
  filter(rowSums(.) != 0)
#1321 asvs, 102 samples

#rename the negative control
colnames(asv16s.exp4)[99] = "neg1"
colnames(asv16s.exp4)[100] = "neg2"
colnames(asv16s.exp4)[101] = "neg3"
colnames(asv16s.exp4)[102] = "neg4"

tax.16s.exp4 <- subset(tax.16s, row.names(tax.16s) %in% rownames(asv16s.exp4)) 

tax.16s.exp4 <- tax.16s.exp4[order(rownames(tax.16s.exp4)), ]
asv16s.exp4 <- asv16s.exp4[order(rownames(asv16s.exp4)), ]

row.names(asv16s.exp4) == row.names(tax.16s.exp4) # sanity check
nrow(asv16s.exp4) #=1321 ASVs
nrow(tax.16s.exp4) #=1321 ASVs

#### exp4_phys1 = USING all data (raw) ####
asv16s.tax1raw <- data.frame(Feature.ID=row.names(tax.16s.exp4),Taxon=tax.16s.exp4[,1])
asv16s.tax2raw <- parse_taxonomy(asv16s.tax1raw)
asv16s.tax3raw <- cbind(ASVs=row.names(asv16s.tax2raw),asv16s.tax2raw)

tree.16s.root.raw <-root(asv16s.tree, "63ffc596781727668a90d6da135a33b4")## root with an archaeon
new_treeraw <- ape::multi2di(tree.16s.root.raw)

TAX.16sraw <- tax_table(as.matrix(asv16s.tax3raw))
OTU.16sraw <- otu_table(as.matrix(asv16s.exp4), taxa_are_rows = TRUE)
SAM.16sraw <- sample_data(meta)
Exp4.physeq1_raw <-merge_phyloseq(phyloseq(OTU.16sraw),SAM.16sraw,TAX.16sraw,new_treeraw)
Exp4.physeq1_raw
#saveRDS(Exp4.physeq1_raw, "RDS_files/Exp4.physeq1_raw.RDS")

#### DECONTAM PACKAGE for identifying contaminants####
decon <- as.data.frame(sample_data(Exp4.physeq1_raw))
decon$LibrarySize <- sample_sums(Exp4.physeq1_raw)
decon <- decon[order(decon$LibrarySize),]
decon$Index <- seq(nrow(decon))
ggplot(data=decon, aes(x=Index, y=LibrarySize, color=treatment)) + geom_point()

# Step 1: Calculate the number of samples in which each ASV is present
present_in_samples <- as.data.frame(apply(asv16s.exp4 > 0, 1, sum))
present_in_samples$sample_number <- present_in_samples[[1]]
present_in_samples <- present_in_samples %>% dplyr::select(sample_number) %>% arrange(sample_number)
hist(present_in_samples)
dim(present_in_samples)#1321 ASV

counts_by_samples <- table(present_in_samples)
counts_by_samples
sorted_counts <- sort(counts_by_samples)
print(sorted_counts)#958 ASVs are only present in one sample
hist(sorted_counts)

#Prevalance Method (this is raw data, no filtering)
sample_data(Exp4.physeq1_raw)$is.neg <- sample_data(Exp4.physeq1_raw)$treatment == "NEG"
contamdf.prev <- isContaminant(Exp4.physeq1_raw, method="prevalence", neg="is.neg")
contamdf.prev
table(contamdf.prev$contaminant) ### 37 ASVs are identified as contaminants, 1284 as not
head(which(contamdf.prev$contaminant), n=37)

ggplot(data=contamdf.prev, aes(x=p)) + 
  labs(x = 'decontam-prevalence Score', y='Number of ASVs') + 
  geom_histogram(binwidth=0.02)

###filter the df to samples id'd as contaminants and find out home many samples they are present in and the read abundance for these
cont.prev.true <-subset(contamdf.prev, contaminant=="TRUE")
abund <- asv16s.exp4
abund$reads <- rowSums(asv16s.exp4)
cont.prev.true.reads <- subset(abund, row.names(abund) %in% rownames(cont.prev.true)) 
cont.true.reads.only <- cont.prev.true.reads %>% dplyr::select(reads)  
cont.prev.true.reads <- cont.prev.true.reads[,-103]

#count up how many samples each ASV is present in
count_non_zero <- function(row) {
  sum(row != 0)
}
cont.prev.true.reads$NonZeroCount <- apply(cont.prev.true.reads[, -1], MARGIN = 1, count_non_zero)

#### Clean up data ####
#these 37 ASV have low frequency and low abundance and were closely associated with the negative controls, they will be removed from the dataset
cont.tax <- subset(tax.16s.exp4, row.names(tax.16s.exp4) %in% rownames(cont.prev.true.reads))

# remove chloroplasts
tax.16s.exp4.2 <- tax.16s.exp4[grep("Chloroplast",tax.16s.exp4$Taxon, invert = T),]
nrow(tax.16s.exp4.2) #1294, removed 27 ASVs

# remove mitochondria
tax.16s.exp4.2 <- tax.16s.exp4.2[grep("Mitochondria",tax.16s.exp4.2$Taxon, invert = T),]
nrow(tax.16s.exp4.2) #1265, removed 29 ASVs

# remove unassigned
tax.16s.exp4.2 <- tax.16s.exp4.2[grep("Unassigned",tax.16s.exp4.2$Taxon, invert = T),]
nrow(tax.16s.exp4.2) #1253, removed 12 ASVs

#remove contaminants identified above
cont.tax <- cont.tax %>% rownames_to_column(var = "ASV")
tax.16s.exp4.2 <- tax.16s.exp4.2 %>% rownames_to_column(var = "ASV")
tax.16s.exp4.2 <- anti_join(tax.16s.exp4.2, cont.tax, by="ASV") 
rownames(tax.16s.exp4.2) <- tax.16s.exp4.2$ASV
tax.16s.exp4.2$ASV<- NULL

#remove negative controls
asv16s.exp4_filt <- asv16s.exp4[,-c(99:102)]
#remove the filtered taxa from the ASV table
asv16s.exp4_filt <- subset(asv16s.exp4_filt, row.names(asv16s.exp4_filt) %in% row.names(tax.16s.exp4.2))
summary(colSums(asv16s.exp4_filt))


#### Exp4.physeq2 = decontaminated data, filtered, negative controls removed ####
row.names(asv16s.exp4_filt) == row.names(tax.16s.exp4.2)
asv16s.tax1raw <- data.frame(Feature.ID=row.names(tax.16s.exp4.2),Taxon=tax.16s.exp4.2[,1])
asv16s.tax2raw <- parse_taxonomy(asv16s.tax1raw)
asv16s.tax3raw <- cbind(ASVs=row.names(asv16s.tax2raw),asv16s.tax2raw)

pruned.tree<- ape::drop.tip(asv16s.tree,asv16s.tree$tip.label[-match(rownames(asv16s.exp4_filt), asv16s.tree$tip.label)])
tree.16s.root.raw <-root(pruned.tree, "63ffc596781727668a90d6da135a33b4")## root with an archaeon
new_treeraw <- ape::multi2di(tree.16s.root.raw)
TAX.16sraw <- tax_table(as.matrix(asv16s.tax3raw))
OTU.16sraw <- otu_table(as.matrix(asv16s.exp4_filt), taxa_are_rows = TRUE)
SAM.16sraw <- sample_data(meta)
Exp4.physeq2 <-merge_phyloseq(phyloseq(OTU.16sraw),SAM.16sraw,TAX.16sraw,new_treeraw)
saveRDS(Exp4.physeq2, "RDS_files/Exp4.physeq2.RDS")
#98 samples, 1216 ASVs

#### Rarify ####
#take out asv less than 10
asv16s.exp4_filt[asv16s.exp4_filt < 10] <- 0# controls for barcode swapping
asv16s.exp4_filt <- asv16s.exp4_filt[as.logical(rowSums(asv16s.exp4_filt != 0)), ] #456 ASV left
tax.16s.exp4_filt <- subset(tax.16s.exp4.2, row.names(tax.16s.exp4.2) %in% row.names(asv16s.exp4_filt))

asv16s.exp4_filt.t <- t(asv16s.exp4_filt)
min(rowSums(asv16s.exp4_filt.t)) #min number of reads per sample = 4431
max(rowSums(asv16s.exp4_filt.t)) #max number of reads per sample = 48534
mean(rowSums(asv16s.exp4_filt.t)) #mean number of reads per sample = 26392.9

set.seed(1115)
asv16s.rt <- rrarefy(asv16s.exp4_filt.t, 4431) ## rarefy at 4431, samples are rows
asv16s.rt <- asv16s.rt[,colSums(asv16s.rt) > 0] ## remove ASVs no longer present, 27 ASV removed
asv16s.rt <- asv16s.rt[order(row.names(asv16s.rt)),] # order samples alphabetically
asv16s.rt <- asv16s.rt[,order(colnames(asv16s.rt))] # order asvs alphabetically
dim(asv16s.rt) #98 samples, 446 ASVs
summary(rowSums(asv16s.rt)) #=4431
summary(colSums(asv16s.rt)) #1

asv16s.r <- t(asv16s.rt) #rows are asvs and columns are samples

tax.16s2.r <- subset(tax.16s.exp4_filt, row.names(tax.16s.exp4_filt) %in% row.names(asv16s.r)) #filter tax table to match asvs
nrow(asv16s.r)#446
nrow(tax.16s2.r)#446

meta.filt <- meta[order(row.names(meta)),] # order samples alphabetically
meta.r <- subset(meta.filt, row.names(meta.filt) %in% colnames(asv16s.r)) #filter meta to relevant samples

#now asv, tax, meta all the same size
row.names(asv16s.rt) == row.names(meta.r)
tax.16s2.r <- tax.16s2.r[order(row.names(tax.16s2.r)),] # order ASVs alphabetically
asv16s.r <- asv16s.r[order(row.names(asv16s.r)),] # order ASVs alphabetically
row.names(asv16s.r) == row.names(tax.16s2.r)

#### 16S alpha diversity ####
shannon.16s <- diversity(asv16s.rt, index="shannon")
ef.16s <- exp(shannon.16s)
summary(ef.16s)
ef.16s <- round(ef.16s)
richness <- colSums(asv16s.r !=0)
md16s <- cbind(meta.r, richness)
md16s <- cbind(md16s, shannon.16s)
md16s <- cbind(md16s, ef.16s)

write.csv(md16s, "data/exp4_md16s.csv", row.names=TRUE)
md16s <- read.csv("data/exp4_md16s.csv", header=TRUE)
rownames(md16s) <- md16s[,1]
md16s <- md16s[order(row.names(md16s)),] # order samples alphabetically
md16s$ID <- rownames(md16s)
row.names(asv16s.rt) == row.names(md16s)

#visualize alpha diversity metrics
prey_order <- c("ant", "beetle", "fly", "lm", "cb", "pos", "water")

# Use the factor function to set the order of levels
md16s$treatment <- factor(md16s$treatment, levels = prey_order)

ggplot(data=md16s, aes(x=day, y=ef.16s, group=treatment, color=treatment)) +
  geom_jitter() + ylab("Effective Number of ASVs (Hill Number 1") + theme_classic()+
  geom_smooth(se=FALSE)+ scale_color_manual(values=c(color_palette))

ggplot(data=md16s, aes(x=day, y=richness, group=treatment, color=treatment)) +
  geom_jitter() + ylab("ASV Richness") + theme_classic()+
  geom_smooth(se=FALSE)+ scale_color_manual(values=c(color_palette))


ggplot(data=md16s, aes(x=treatment, y=richness, fill=treatment, color=treatment)) +
  geom_boxplot(outlier.shape=NA,alpha=.4) +geom_jitter() + ylab("ASV Richness") + theme_classic()+
  scale_color_manual(values=c(color_palette))+scale_fill_manual(values=c(color_palette))+
  facet_wrap(~day, nrow=1)+ theme(legend.position="none")

md16s.filt <- filter(md16s, treatment %in% c("ant", "beetle", "fly"))
md16s.filt2 <- filter(md16s, treatment %in% c("lm", "cb"))

ggplot(data=md16s.filt, aes(x=day, y=ef.16s, group=treatment, color=treatment)) +
  geom_jitter() + ylab("Effective Number of ASVs (Hill Number 1)") + theme_classic()+
  geom_smooth(se=FALSE)+ scale_color_manual(values=c(color_palette))

ggplot(data=md16s.filt, aes(x=day, y=richness, group=treatment, color=treatment)) +
  geom_jitter() + ylab("ASV Richness") + theme_classic()+
  geom_smooth(se=FALSE)+ scale_color_manual(values=c(color_palette))

ggplot(data=md16s.filt2, aes(x=treatment, y=richness, group=treatment, color=treatment, fill=treatment)) +
   ylab("ASV Richness") + theme_classic()+
  geom_boxplot(outlier.shape=FALSE, alpha=.7)+ geom_jitter(size=4)+
  scale_color_manual(values=c("#00446a", "#d88a00")) +
  facet_wrap(~day)+scale_fill_manual(values=c("#00446a", "#d88a00"))+theme(legend.position="none")


#### alpha diveristy GLMM
md16s.filt$treatment <- factor(md16s.filt$treatment, levels=c("ant", "beetle", "fly"))

ggplot(data=md16s.filt, aes(x=treatment, y=richness, fill=treatment, color=treatment)) +
 geom_boxplot(outlier.shape=NA,alpha=.4) +geom_jitter() + ylab("Richness (Total ASVs)") + theme_classic()+
  scale_color_manual(values=c(color_palette))+scale_fill_manual(values=c(color_palette))+
  facet_wrap(~day, nrow=1)+ theme(legend.position="none")+
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("")

ggplot(data=md16s.filt, aes(x=treatment, y=ef.16s, fill=treatment, color=treatment)) +
  geom_boxplot(outlier.shape=NA,alpha=.4) +geom_jitter() + ylab("Effective Number of ASVs (Hill Number 1)") + theme_classic()+
  scale_color_manual(values=c(color_palette))+scale_fill_manual(values=c(color_palette))+
  facet_wrap(~day, nrow=1)+ theme(legend.position="none")

md16s.filt$treatment <- relevel(md16s.filt$treatment, ref = "fly")
md16s.filt$day <- as.factor(md16s.filt$day)

glm2 <- brm(richness ~ treatment + day, family = negbinomial(), data=md16s.filt, iter=10000)

pglm2a <- ggeffects::ggpredict(glm2, "treatment")
pglm2b <- ggeffects::ggpredict(glm2, "day")

posteriorglm2 <- mcmc_intervals_data(glm2, 
                                  prob_outer=0.95,
                                  prob=0.95)

posteriorglm2$nonzero <- NA
posteriorglm2$nonzero[posteriorglm2$ll>0 & posteriorglm2$hh>0] <- "nonzero"
posteriorglm2$nonzero[posteriorglm2$ll<0 & posteriorglm2$hh<0] <- "nonzero"
posteriorglm2$nonzero[is.na(posteriorglm2$nonzero)] <- "zero"
posteriorglm2<- posteriorglm2[1:5,]

#FIGS2
ggplot(posteriorglm2, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_classic() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  guides(linetype=FALSE) 

### alpha lm vs cb
glm3 <- brm(richness ~ treatment + day, family = negbinomial(), data=md16s.filt2, iter=10000)

pglm3a <- ggeffects::ggpredict(glm3, "treatment")
pglm3b <- ggeffects::ggpredict(glm3, "day")

posteriorglm3 <- mcmc_intervals_data(glm3, 
                                     prob_outer=0.95,
                                     prob=0.95)
mcmc_plot(glm3, regex_pars="b_",prob_outer=0.95,
                              prob=0.95)
posteriorglm3$nonzero <- NA
posteriorglm3$nonzero[posteriorglm3$ll>0 & posteriorglm3$hh>0] <- "nonzero"
posteriorglm3$nonzero[posteriorglm3$ll<0 & posteriorglm3$hh<0] <- "nonzero"
posteriorglm3$nonzero[is.na(posteriorglm3$nonzero)] <- "zero"
posteriorglm3<- posteriorglm3[1:3,]

#FIG6
ggplot(posteriorglm3, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_classic() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  guides(linetype=FALSE) 

#### MAKE PHYLOSEQ OBJECTS ####
#Exp4.physeq3.RDS = RARIFIED
md16s<- md16s[order( md16s[,3], md16s[,6] ),]
asv16s.rt <- asv16s.rt[match(rownames(md16s), rownames(asv16s.rt)), ]
asv16s.r <- t(asv16s.rt)
tax.16s2.r <- tax.16s2.r[match(rownames(asv16s.r), rownames(tax.16s2.r)), ]
row.names(asv16s.r) == row.names(tax.16s2.r)
r.asv16s.tax1 <- data.frame(Feature.ID=row.names(tax.16s2.r),Taxon=tax.16s2.r[,1])
r.asv16s.tax2 <- parse_taxonomy(r.asv16s.tax1)
r.asv16s.tax3 <- cbind(ASVs=row.names(r.asv16s.tax2),r.asv16s.tax2)

r.tree.16s <- picante::prune.sample(asv16s.rt, asv16s.tree)# Subset tree to relevant samples
r.tree.16s.root <-root(r.tree.16s, "9a69e8b7438b1b1c7a04d6493e338d54")## root with an myc bac
r.new_tree <- ape::multi2di(r.tree.16s.root)

r.TAX.16s <- tax_table(as.matrix(r.asv16s.tax3))
r.OTU.16s <- otu_table(as.matrix(asv16s.r), taxa_are_rows = TRUE)

r.SAM.16s <- sample_data(md16s)
Exp4.physeq3 <-merge_phyloseq(phyloseq(r.OTU.16s),r.SAM.16s,r.TAX.16s,r.new_tree)
Exp4.physeq3
#446 taxa and 98 samples
saveRDS(Exp4.physeq3, "RDS_files/Exp4.physeq3.RDS")

#### BETA DIVERSITY ####
#no controls
Exp4.physeq3 <- readRDS("RDS_files/Exp4.physeq3.RDS")
Exp4.set1.physeq3 = subset_samples(Exp4.physeq3, treatment %in% c("ant", "beetle", "fly") )
row_sums <- rowSums(otu_table(Exp4.set1.physeq3))
nonzero_rows <- row_sums != 0
Exp4.set1.physeq3 <- prune_taxa(nonzero_rows, Exp4.set1.physeq3)

## Calculate weighted Unifrac distance and run NMDS
wu.dist.16s<- distance(Exp4.set1.physeq3,"uUniFrac")

set.seed(123)
wu.nmds.16s <- metaMDS(wu.dist.16s,
                             k = 2, 
                             trymax = 1000,
                             wascores = TRUE)

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta_phys <- sample_data(Exp4.set1.physeq3)
meta_phys$ID <- rownames(meta_phys)
meta_phys <- data.frame(meta_phys)
data_merge2 <- merge(data.scores2, meta_phys, by = c("ID"))
data_merge2$day <- as.factor(data_merge2$day)

plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette_3)[data_merge2$treatment],
     pch=19, main="uUniFrac all time points")
ordiellipse(wu.nmds.16s, groups = as.factor(data_merge2$treatment),col = c(color_palette_3))
legend("topright", 
       legend=c("ant", "beetle","fly"),
       col= c(color_palette_3),
       pch=19,
       cex=1,
       bty = "n")


plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette_3)[data_merge2$treatment],
     pch=19)
legend("topright", 
       legend=c("ant", "beetle","fly"),
       col= c(color_palette_3),
       pch=19,
       cex=1,
       bty = "n")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "ant", col = "#c24000")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "beetle", col = "#008388")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "fly", col = "#5d004f")

#get the data scores from the nmds
all.scores.beta <- as.data.frame(scores(wu.nmds.16s, "sites"))
all.scores.beta <- cbind(all.scores.beta, meta_phys)
all.scores.beta$day <- as.factor(all.scores.beta$day)
all.scores.beta <- all.scores.beta %>% 
  unite(hull.id, treatment,remove = FALSE)

hullsall <- all.scores.beta %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

ggplot(all.scores.beta, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2") + aes(fill = factor(treatment)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+
  scale_color_manual(values=color_palette)+scale_fill_manual(values=color_palette)

#get the data scores from the nmds
all.scores.beta <- as.data.frame(scores(wu.nmds.16s, "sites"))
all.scores.beta <- cbind(all.scores.beta, meta_phys)
all.scores.beta$day <- as.factor(all.scores.beta$day)
all.scores.beta <- all.scores.beta %>% 
  unite(hull.id, day,remove = FALSE)

hullsall <- all.scores.beta %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

ggplot(all.scores.beta, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = day)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "day", y = "NMDS2") + aes(fill = factor(day)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+
  scale_color_manual(values=c("navyblue", "skyblue", "slategray"))+scale_fill_manual(values=c("navyblue", "skyblue", "slategray"))



#### Beta Diversity Analysis
uu.bd <- betadisper(wu.dist.16s, data_merge2$treatment)
uu.bd
boxplot(uu.bd)
anova(uu.bd)


uu.bd2 <- betadisper(wu.dist.16s, data_merge2$day)
uu.bd2
boxplot(uu.bd2)
anova(uu.bd2)

# PERMANOVA for categorical variables (factors)
# set number of permutations
perm <- how(nperm = 999)
ad.16s.treat <- adonis2(wu.dist.16s ~ treatment + day, data=data_merge2, permutations = perm, by="margin")
ad.16s.treat


wu.ad.pw <- pairwise.adonis2(wu.dist.16s ~ treatment, data = data_merge2, strata = 'day')
wu.ad.pw
write.csv(wu.ad.pw, "data/wu.ad.pw.csv")

#### Beta diversity LM vs CB ####
Exp4.set2.physeq3 = subset_samples(Exp4.physeq3, treatment %in% c("lm", "cb") )
row_sums <- rowSums(otu_table(Exp4.set2.physeq3))
nonzero_rows <- row_sums != 0
Exp4.set2.physeq3 <- prune_taxa(nonzero_rows, Exp4.set2.physeq3)

## Calculate weighted Unifrac distance and run NMDS
wu.dist.16s2<- distance(Exp4.set2.physeq3,"uUniFrac")

set.seed(123)
wu.nmds.16s2 <- metaMDS(wu.dist.16s2,
                       k = 2, 
                       trymax = 1000,
                       wascores = TRUE)#0.115731

## Plot NMDS
meta_phys2 <- sample_data(Exp4.set2.physeq3)
meta_phys2$ID <- rownames(meta_phys2)
meta_phys2 <- data.frame(meta_phys2)
all.scores.beta2 <- as.data.frame(scores(wu.nmds.16s2, "sites"))
all.scores.beta2 <- cbind(all.scores.beta2, meta_phys2)
all.scores.beta2$day <- as.factor(all.scores.beta2$day)
all.scores.beta2 <- all.scores.beta2 %>% 
  unite(hull.id, treatment,remove = FALSE)

hullsall <- all.scores.beta2 %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

ggplot(all.scores.beta2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2") + aes(fill = factor(treatment)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+
  scale_color_manual(values=c("#00446a", "#d88a00"))+scale_fill_manual(values=c("#00446a", "#d88a00"))+theme(legend.positio="none")


ggplot(all.scores.beta2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = day)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "day", y = "NMDS2")+
  scale_color_manual(values=c("navyblue", "skyblue", "slategray"))

#### Beta Diversity Analysis
uu.bd <- betadisper(wu.dist.16s2, all.scores.beta2$treatment)
uu.bd
boxplot(uu.bd)
anova(uu.bd)

# PERMANOVA
# set number of permutations
perm <- how(nperm = 999)
ad.16s.treat2 <- adonis2(wu.dist.16s2 ~ treatment + day, data=all.scores.beta2, permutations = perm, by="margin")
ad.16s.treat2

#### Top 16S ASVs, with taxonomy #####
asv16s.top <- data.frame(sort(colSums(asv16s.rt),decreasing = T))
names(asv16s.top) <- "Sequences"
asv16s.top.tax <- subset(tax.16s2.r, row.names(tax.16s2.r) %in% row.names(asv16s.top))
asv16s.top.tax <- asv16s.top.tax[match(rownames(asv16s.top), rownames(asv16s.top.tax)), ]
asv16s.top <- data.frame(asv16s.top,asv16s.top.tax)

####look at top 20 asv abundance taxonomy in 16S
sum(asv16s.top$Sequences)
asv16s.top.20 <- asv16s.top[1:20,]
sum(asv16s.top.20$Sequences)

Feature.ID <- rownames(asv16s.top.20)
asv16s.top.20 <- cbind(asv16s.top.20, Feature.ID)
asv16s.top.20.tax <- parse_taxonomy(asv16s.top.20)
rownames(asv16s.top.20.tax) == rownames(asv16s.top.20)
asv16s.top.20.tax.abund <- cbind(asv16s.top.20.tax, asv16s.top.20)

#look at top 16S genera
asv.16S.mv <- subset(asv16s.top.20.tax.abund, Genus == "Pseudomonas")
((sum(asv.16S.mv$Sequences)/287729)*100)
((sum(asv.16S.mv$Sequences)/434238)*100)

asv.16S.st <- subset(asv16s.top.20.tax.abund, Genus == "Pedobacter")
((sum(asv.16S.st$Sequences)/287729)*100)
((sum(asv.16S.st$Sequences)/434238)*100)

asv.16S.pd <- subset(asv16s.top.20.tax.abund, Genus == "Undibacterium")
((sum(asv.16S.pd$Sequences)/287729)*100)
((sum(asv.16S.pd$Sequences)/434238)*100)


#### ant, beetle, fly only####
asv16s.top <- asv_glmr
asv16s.top <- data.frame(sort(colSums(asv16s.top),decreasing = T))
names(asv16s.top) <- "Sequences"
asv16s.top.tax <- subset(tax.16s2.r, row.names(tax.16s2.r) %in% row.names(asv16s.top))
asv16s.top.tax <- asv16s.top.tax[match(rownames(asv16s.top), rownames(asv16s.top.tax)), ]
row.names(asv16s.top.tax) == row.names(asv16s.top.tax)
asv16s.top <- data.frame(asv16s.top,asv16s.top.tax)

####look at top 20 asv abundance taxonomy in 16S
sum(asv16s.top$Sequences) #186102 16S reads post clean up
asv16s.top.20 <- asv16s.top[1:20,]
sum(asv16s.top.20$Sequences) #135186 reads in top 20

#top 20 asv reads divided by total reads
(135186/186102)*100 #72.64081

Feature.ID <- rownames(asv16s.top.20)
asv16s.top.20 <- cbind(asv16s.top.20, Feature.ID)
asv16s.top.20.tax <- parse_taxonomy(asv16s.top.20)
rownames(asv16s.top.20.tax) == rownames(asv16s.top.20)
asv16s.top.20.tax.abund <- cbind(asv16s.top.20.tax, asv16s.top.20)

#look at top 16S families
asv.16S.mv <- subset(asv16s.top.20.tax.abund, Genus == "Pseudomonas")
((sum(asv.16S.mv$Sequences)/135186)*100) #19.70248%
((sum(asv.16S.mv$Sequences)/186102)*100) #14.31204%

asv.16S.st <- subset(asv16s.top.20.tax.abund, Genus == "Pedobacter")
((sum(asv.16S.st$Sequences)/135186)*100) #22.97649%
((sum(asv.16S.st$Sequences)/186102)*100) #16.69031%

asv.16S.pd <- subset(asv16s.top.20.tax.abund, Genus == "uncultured")
((sum(asv.16S.pd$Sequences)/135186)*100) #8.792331%
((sum(asv.16S.pd$Sequences)/186102)*100) #6.38682%

 #### ANCOMBC DIFFERENTIAL ABUNDANCE ####
#subset to treatment
Exp4.physeq3_filt <- Exp4.physeq3 %>% subset_samples(treatment %in% c("ant", "beetle", "fly"))
row_sums <- rowSums(otu_table(Exp4.physeq3_filt))
nonzero_rows <- row_sums != 0
Exp4.physeq3_filt <- prune_taxa(nonzero_rows, Exp4.physeq3_filt)

tse = mia::makeTreeSummarizedExperimentFromPhyloseq(Exp4.physeq3_filt)
tse$treatment = factor(tse$treatment, levels = c("fly","ant", "beetle"))

output = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                  fix_formula = "treatment + day",
                  #rand_formula = "(1 | day)",
                  p_adj_method = "fdr", pseudo_sens = TRUE,
                  prv_cut = 0.07,
                  group = "treatment",
                  alpha = 0.05, n_cl = 3, verbose = TRUE,
                  global = TRUE, pairwise = TRUE)

res_prim = output$res
res_pair <- output$res_pair

output_res <- output$res
output_res_sig <- subset(output_res, `diff_(Intercept)`==T | diff_treatmentbeetle==T |
                           diff_treatmentant==T )
output_res_pair_sig <- subset(res_pair, `diff_treatmentant`==T | diff_treatmentbeetle==T |
                                diff_treatmentbeetle_treatmentant==T )
dim(output_res)
dim(output_res_sig)
dim(output_res_pair_sig)

res_prim_filt <- output_res_pair_sig[,c(1:7)]
colnames(res_prim_filt)[2] <- "lfc_antvsfly"
colnames(res_prim_filt)[3] <- "lfc_beetlevsfly"
colnames(res_prim_filt)[4] <- "lfc_beetlevsant"
colnames(res_prim_filt)[5] <- "se_antvsfly"
colnames(res_prim_filt)[6] <- "se_beetlevsfly"
colnames(res_prim_filt)[7] <- "se_beetlevsant"

res_prim_filt_melt <- res_prim_filt %>%
  gather(key, value, -taxon) %>%
  separate(key, into = c("measure", "treatment"), sep = "_") %>%
  spread(measure, value)

ggplot(data=res_prim_filt_melt, aes(y=lfc, x=taxon,  color=treatment)) + geom_boxplot()+
  theme_classic() + facet_wrap(~treatment,nrow=3)+geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2)

#taxonomy
tax <- as.data.frame(tax.16s2.r)
tax <- cbind(Feature.ID=rownames(tax),tax)
tax.p <- parse_taxonomy(tax)
tax.p$ASV <- row.names(tax.p)

#asvs
asv16s.anc <- asv16s.r
asv16s.anc <- as.data.frame(asv16s.anc)
asv16s.anc$ASV <- row.names(asv16s.anc)

#ancombc
names(output_res_pair_sig)[1] <- "ASV"
ancom_tax <- merge(output_res_pair_sig, tax.p, by="ASV")

asv16s.rt.127 <- asv16s.anc[, (grep(pattern="*\\_DAY7|*\\_DAY14|*\\_DAY35", colnames(asv16s.anc)))]
asv16s.rt.127$ASV <- row.names(asv16s.rt.127)
wkB_ancom_tax <- merge(ancom_tax, asv16s.rt.127, by="ASV")

wkB_ancom_tax.filt <- wkB_ancom_tax[,c(1,27:28, 30:127)]
wkB_ancom_tax.filt.melt <- wkB_ancom_tax.filt %>% melt(id=c("ASV", "Family", "Genus"))
AN_meta_127<- md16s.filt[,c(1,6,7)]
AN_meta_127["variable"] <- AN_meta_127[["X"]]
AN_meta_127 <- AN_meta_127[,-1]

wkB_ancom_tax.filt.melt <- left_join(wkB_ancom_tax.filt.melt, AN_meta_127, by="variable")
wkB_ancom_tax.filt.melt$treatment <- factor(wkB_ancom_tax.filt.melt$treatment, levels=prey_order)
wkB_ancom_tax.filt.melt$ASV <- as.factor(wkB_ancom_tax.filt.melt$ASV)
wkB_ancom_tax.filt.melt <- wkB_ancom_tax.filt.melt[complete.cases(wkB_ancom_tax.filt.melt$treatment), ]

ggplot(wkB_ancom_tax.filt.melt, aes(x=ASV, y=log2(0.5+value), fill=treatment)) +
  geom_boxplot(outlier.shape=NA, alpha = .7)+geom_jitter() +theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))+
  scale_fill_manual(values=color_palette)+facet_wrap(~day, ncol=1)

wkB_ancom_tax.filt.melt$day <- as.numeric(wkB_ancom_tax.filt.melt$day)
ggplot(wkB_ancom_tax.filt.melt, aes(x = day, y = log2(0.5 + value), fill = treatment, color = treatment)) +
  geom_jitter(position = position_jitterdodge()) +  # Adjust position here
  theme_classic() + geom_smooth(se=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ASV, nrow = 2) +
  scale_color_manual(values = color_palette)

##filter to a few
wkB_ancom_tax.filt.melt.sub <- filter(wkB_ancom_tax.filt.melt, ASV %in% c("9dcd05e2f900431a5c267a4b853bc2ae", "42a90f85f75e17da09d9ac1106b132f0", "525d33a7c4702329607bab16327e2cb7", "20b3df08cd5fac2be26a2928155bab69", "38c27ceaed634984c1225a82648cf571", "9154281b5cef5092ba03398458099017", "c6fa5fd77a5483314afc1d62c03864eb"))
wkB_ancom_tax.filt.melt.sub$treatment <- as.factor(wkB_ancom_tax.filt.melt.sub$treatment)

ggplot(wkB_ancom_tax.filt.melt.sub, aes(x = ASV, y = log2(0.5 + value), fill = treatment, color = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + # Boxplots
  geom_jitter(position = position_jitterdodge()) + # Jittered points centered around boxplots
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = color_palette_nopos) + # Fill color scale
  facet_wrap(~day, ncol = 1) + # Faceting
  scale_color_manual(values = color_palette_nopos)

#### ancombc lm vs cb ####
Exp4.physeq3_filt2 <- Exp4.physeq3 %>% subset_samples(treatment %in% c("lm", "cb"))
row_sums2 <- rowSums(otu_table(Exp4.physeq3_filt2))
nonzero_rows2 <- row_sums2 != 0
Exp4.physeq3_filt2 <- prune_taxa(nonzero_rows2, Exp4.physeq3_filt2)

output2 <- ancombc2(data = Exp4.physeq3_filt2, assay_name = "counts", tax_level = NULL,
                   fix_formula = "treatment",
                   p_adj_method = "fdr",
                   prv_cut = 0.07,
                   group = "treatment",
                   alpha = 0.05, verbose = TRUE,
                   global = TRUE)

res_prim2 <- output2$res
output_res_sig2 <- subset(res_prim2, `diff_(Intercept)`==T | diff_treatmentcb==T)
dim(res_prim2)
dim(output_res_sig2)

#### RELATIVE ABUNDANCE PLOTS ####
tax <- data.frame(tax_table(Exp4.physeq3))
samp <- data.frame(sample_data(Exp4.physeq3))
otu <- data.frame(otu_table(Exp4.physeq3))
rownames(otu)==rownames(tax)

jewel_tone_palette <- c("#000000","#a9a19c","#493829","#117a65","#bdd09f", "#0c2461", "#d35400",
                        "#1f8a13","#8f3b1b","#855723","#e67e22","#740058",
                        "palegoldenrod","#3498db","#f39c12", "#28b463","#af7ac5",
                        "#1f618d","#d35400","#154360","darkmagenta")
# Plotting the relative abundance at Genus level
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])

otu.gaop <- otu.gao[c(1:20),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
other <- colSums(otu.gao[c(21:nrow(otu.gao)),])#put the all the other genera in a "other" category
otu.gaopo <- rbind(otu.gaop, other)
otu.gaopo_gg <- as.data.frame(otu.gaopo)
otu.gaopo_ggp <- adorn_percentages(otu.gaopo_gg,, 1:ncol(otu.gaopo_gg), denominator="col")#number of columns (samples), calculate RA
other <- otu.gaopo_ggp[21,]
otu.gaopo_ggp <- otu.gaopo_ggp[-21,]
otu.gaopo_ggp <- otu.gaopo_ggp[order(rowSums(otu.gaopo_ggp),decreasing = F),]
otu.gaopo_ggp <- rbind(other, otu.gaopo_ggp)
otu.gaopo_ggp[ "family" ] <- rownames(otu.gaopo_ggp)
otu.gaopo_ggp$family <- factor(otu.gaopo_ggp$family, levels = otu.gaopo_ggp$family)
otu.gaopo_ggp.rmelt<- reshape2::melt(otu.gaopo_ggp, id.vars="family", value.name="Relative_Abundance", variable.name="Sample")
df3 <- left_join(otu.gaopo_ggp.rmelt, samp, by=c("Sample"="X"))
df3$day <- as.factor(df3$day)
df3_filt <- df3 %>% filter(., treatment %in% c("ant", "beetle", "fly", "lm", "cb", "pos"))
ggplot(df3_filt, aes(x=day, y=Relative_Abundance, fill = family)) + 
  geom_bar( position = "fill", stat = "identity", width=1) + 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values = jewel_tone_palette) + facet_wrap(~treatment, nrow=1)+
  guides(fill = guide_legend(ncol = 7)) + theme(legend.position="bottom")

#### Protozoa Community Analysis ####
protist <- read.csv("data/protist_matrix.csv", header=TRUE)
protist$sample_id <- protist$X
protist <- protist[,-1]
prot_meta <- data_filt
prot_meta <- prot_meta[,1:7]
prot_meta <- left_join(protist, prot_meta, by=c("sample_id" = "sample_id"))
rownames(prot_meta) <- prot_meta$sample_id

#get just the community matrix
rownames(prot_meta) <- prot_meta$sample_id
prot_com <- prot_meta[,1:15]

shannon.prot <- diversity(prot_com, index="shannon")
ef.prot <- exp(shannon.prot)
summary(ef.prot)
ef.prot <- round(ef.prot)
richness <- rowSums(prot_com !=0)
prot_meta <- cbind(prot_meta,shannon.prot, ef.prot, richness)
prot_meta <- as.data.frame(prot_meta)

#visualize alpha diversity metrics for protists
ggplot(data=prot_meta, aes(x=day, y=richness, group=treatment, color=treatment)) +
  geom_jitter() + geom_smooth(se=FALSE)+ylab("Richness") + theme_classic()+
  scale_color_manual(values=color_palette)+ facet_wrap(~treatment, nrow=1)+scale_color_manual(values=color_palette)

ggplot(data=prot_meta, aes(x=day, y=ef.prot, group=treatment, color=treatment)) +
  geom_jitter() + geom_smooth(se=FALSE)+ylab("Effective Number of Species") + theme_classic()+
  scale_color_manual(values=color_palette)+ facet_wrap(~treatment, nrow=1)+scale_color_manual(values=color_palette)

prot_meta$treatment <- relevel(prot_meta$treatment, ref = "water")
prot_meta$treatment <- as.factor(prot_meta$treatment)
prot_meta$day <- as.factor(prot_meta$day)

prot_glm1 <- brm(richness ~ treatment +(1|day), family = negbinomial(), data=prot_meta, iter=10000)
mcmc_plot(prot_glm1, regex_pars="b_")

#### Protist BETA DIVERSITY ####
prot_com.filt <- prot_com[rowSums(prot_com) != 0, ]
prot_meta.filt <- prot_meta[match(rownames(prot_com.filt), rownames(prot_meta)), ]
prot.dist<- vegdist(prot_com.filt, method="jaccard", binary=TRUE)

# Perform PCoA
pcoa_result <- cmdscale(prot.dist, k = 2, eig = TRUE)
coordinates <- pcoa_result$points  # Coordinates of samples in PCoA space
groups <- factor(c("ant", "beetle", "fly", "lm", "cb", "pos", "water"))
points(coordinates, col = groups)

coordinates1 <- as.data.frame(coordinates)
coordinates1$sample_id <- rownames(coordinates1)
merged_data <- merge(as.data.frame(coordinates1),prot_meta.filt, by.y = "sample_id", by.x = "sample_id")
ggplot(merged_data, aes(x = V1, y = V2, color = treatment)) +
  geom_jitter(size=3, alpha=.6) +
  labs(title = "PCoA Protist", x = "PC1", y = "PC2") +
  theme_classic() + scale_color_manual(values=color_palette)

## Protist Statistical Analysis
protall.bd <- betadisper(vegdist(prot.dist), prot_meta.filt$treatment)
protall.bd
boxplot(protall.bd)
anova(protall.bd)

proall.bdweek <- betadisper(vegdist(prot.dist), prot_meta.filt$day)
proall.bdweek
boxplot(proall.bdweek)
anova(proall.bdweek)

# PERMANOVA for categorical variables (factors)
# set number of permutations
perm <- how(nperm = 999)
#specify a random variable (plant).
setBlocks(perm) <- with(prot_meta.filt, plant)
prot_adonis <- adonis2(prot.dist ~ treatment*day, data=prot_meta.filt, permutations = perm)
prot_adonis

# pairwise adonis
protall.ad.pw <- pairwise.adonis2(prot.dist ~ treatment * day, data = prot_meta.filt, strata = 'plant')
protall.ad.pw

#### PLANT MORPHOLOGY ####
morph <- read.csv("data/Prey_Exp4_morph.csv", header=TRUE)

ggplot(morph, aes(x=treatment, y=drymass_g, color=treatment)) + geom_violin(trim=FALSE)+
  geom_jitter()+theme_classic()+
  scale_color_manual(values=c(color_palette))

ggplot(morph, aes(x=treatment, y=length_cm, color=treatment)) + geom_violin(trim=FALSE)+
  geom_jitter()+theme_classic()+
  scale_color_manual(values=c(color_palette))
