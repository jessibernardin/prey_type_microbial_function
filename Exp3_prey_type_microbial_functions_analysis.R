#### Prey type influence microbial community structure and function ####
#### Authors: Jessica R. Bernardin, Sarah M. Gray, Leonora S. Bittleston ####
#### last update : December 5, 2023 ####
#### Prey Type ~ Microbial Functions Analysis

#### Load Required Packages ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "tidyverse", "effects", "growthcurver",
  "plyr", "dplyr", "reshape", "reshape2", "ape", "DiagrammeR", "tidyverse",
  "tidybayes", "coefplot", "standardize", "bayesplot", "MCMCvis", "car",
  "patchwork", "ggpubr", "corrr", "ggcorrplot", "factoextra", "MASS",
  "pairwiseAdonis", "plotrix", "gridExtra", "multcompView", "ggeffects", "this.path", "brms", "pracma",
  "phyloseq", "picante", "ape", "qiime2R", "decontam", "janitor", "readr"
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
color_palette <- c(ant_color, beetle_color, fly_color, 
                   lm_color, cb_color,
                   pos_color, water_color)


#### Insect Nutrient Compositions ####
#read in the data
prey_data <- read.csv("Insect_Nutrient_metadata.csv", header=TRUE)

#lipid
ggplot(data=prey_data, aes(x=prey, y=percent_bodymass_is_lipid, color=prey, prey)) + 
  geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=15))+
  scale_color_manual(values=c("#c24000", "#008388" ,"#5d004f"))

#exoskeleton
ggplot(data=prey_data, aes(x=prey, y=percent_body_mass_is_exoskeleton, color=prey, prey)) + 
  geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=15))+
  scale_color_manual(values=c("#c24000", "#008388" ,"#5d004f"))

#graph Lowry
ggplot(data=prey_data, aes(x=prey, y=(LOWRY_protein_per_unit_insect_mass_mg*100), color=prey, prey)) + 
  geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=15))+
  scale_color_manual(values=c("#c24000", "#008388" ,"#5d004f"))


m1_exo <- brm(percent_body_mass_is_exoskeleton~prey, data=prey_data)
summary(m1_exo)
mcmc_areas(m1_exo, regex_pars="b_") + geom_vline(xintercept = 0, linetype=3)+ggtitle("%exoskeleton~preytype")

m2_lipid <- brm(percent_bodymass_is_lipid~prey, data=prey_data)
summary(m2_lipid)
mcmc_areas(m2_lipid, regex_pars="b_")+ geom_vline(xintercept = 0, linetype=3)+ggtitle("%lipid~preytype")

m3_low <- brm((LOWRY_protein_per_unit_insect_mass_mg*100)~prey, data=prey_data)
summary(m3_low)
mcmc_areas(m3_low, regex_pars="b_")+ geom_vline(xintercept = 0, linetype=3)+ggtitle("%protein~preytype")

#relative abundance of prey nutrients
prey <- read.csv("prey_nutrient_relative_abundance.csv", header=TRUE)
prey_melt<- reshape2::melt( prey)

ggplot(prey_melt, aes( x = variable, y = value, fill = X)) + geom_bar( position = "fill", stat = "identity", width=1) + 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values=c("#f4f4f4", "#097969", "#e6bc00", "#8b295a"))













#### Read in Data
#bacterial and plant metadata
data <- read.csv("Prey_Exp4_meta.csv", header = TRUE)

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

# Use the factor function to set the order of levels
data_filt$treatment <- factor(data_filt$treatment, levels = prey_order)


#ecoplate data
ecodat <- read.csv("Exp4_Ecoplate_all.csv", header = TRUE)

#### Chitinase Activity ####
ggplot(data=data_filt, aes(x=day, y=chitinase_uM_min, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Chitinase Activity") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

data_filt$sqrt_chit <- sqrt(data_filt$chitinase_uM_min)

ggplot(data=data_filt, aes(x=day, y=sqrt(chitinase_uM_min), color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("sqrt Chitinase Activity") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

ggplot(data=data_filt, aes(x=treatment, y=cumulative_chit, color = treatment, group=treatment))+ 
  geom_boxplot(outlier.shape=NA) + geom_jitter()+scale_fill_manual(values=c(color_palette))+
  ylab("Cumulative Chitinase Rate") + theme_classic()+scale_color_manual(values=c(color_palette))

#### Protease Activity ####
ggplot(data=data_filt, aes(x=day, y=protease_nM_min, color = treatment, group=treatment))+ 
  geom_smooth() + geom_jitter()+facet_wrap(~treatment) +
  ylab("Protease Rate") + theme_classic()

ggplot(data=data_filt, aes(x=day, y=protease_nM_min, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE) + geom_jitter() +
  ylab("Protease Rate") + theme_classic()

ggplot(data=data_filt, aes(x=treatment, y=cumulative_prot, fill = treatment, group=treatment))+ 
  geom_boxplot() + geom_jitter()+scale_fill_manual(values=c(color_palette))+
  ylab("Cumulative Chitinase Rate") + theme_classic()

ggplot(data=data_filt, aes(x=day, y=protease_nM_min, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Protease Activity") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

ggplot(data=data_filt, aes(x=treatment, y=cumulative_prot, color = treatment, group=treatment))+ 
  geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=.5) +
  ylab("Cumulative Protease Activity") + theme_classic()+scale_color_manual(values=c(color_palette))

#### Bacterial Abundance ####
data_filt$dna_con_ug_ml <- as.numeric(data_filt$dna_con_ug_ml)
ggplot(data=data_filt, aes(x=day, y=dna_con_ug_ml, color = treatment, group=treatment))+ 
  geom_jitter() +facet_wrap(~treatment, nrow=1) +
  ylab("DNA concentration (ng/mL)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

ggplot(data=data_filt, aes(x=day, y=flow_cyt_livingevents_uL, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Total number of living cells (per uL)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")+ ylim(c(0, 20000))

ggplot(data=data_filt, aes(x=day, y=flow_cyt_totalevents_uL, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE)+geom_jitter(alpha=.5) +facet_wrap(~treatment, nrow=1) +
  ylab("Total number of cells (live + dead per uL)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")+ ylim(c(0, 20000))

#### pH ####
ggplot(data=data_filt, aes(x=day, y=pH, color = treatment, group=treatment))+ 
  geom_smooth(se=FALSE) +geom_jitter(alpha=.5)+ facet_wrap(~treatment, nrow=1) +
  ylab("Pitcher Fluid pH") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")

#### Decomposition ####
decom <- data_filt
decom <- decom[complete.cases(decom$prey_mass_g), ]
decom$prey_mass_g <- as.numeric(decom$prey_mass_g)
decom$treatment <- factor(decom$treatment, levels = prey_order)


ggplot(data=decom, aes(x=day, y=prey_mass_g, color = treatment, group=treatment))+ 
  geom_jitter() +facet_wrap(~treatment, nrow=1) +
  ylab("Change in Prey Mass (g)") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")


ggplot(data=decom, aes(x=day, y=prey_mass_g , color = treatment, group=treatment))+ 
  geom_smooth() + geom_jitter()+
  ylab("Prey mass (g)") +theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(text = element_text(size=20))

decom_delta <- data_filt[complete.cases(data_filt$delta_prey_percent), ]
decom_delta$treatment <- factor(decom_delta$treatment, levels = prey_order)

ggplot(data=decom_delta, aes(x=treatment, y=delta_prey_percent , color = treatment, group=treatment))+ 
  geom_boxplot()+ geom_jitter()+
  ylab("Proportion Mass Lost") + theme_classic()+scale_color_manual(values=c(color_palette))+
  theme(text = element_text(size=20))

#### ECOPLATES ####
# no blanks, experimental controls, or damaged plants
eco_filt <- subset(ecodat, treatment != "water" & treatment != "pos")
com <- eco_filt
com_meta <- subset(com, select= c(sample, plant_block, treatment, week))
com <- subset(com, select= -c(sample, plant_block, treatment, week))
com[com < 0] <- 0

#run NMDS on ecoplate data, no controls
set.seed(123)
nmds <- metaMDS(com, k=2, trymax=500) #stress=0.1407289

#get the data scores from the nmds
all.scores <- as.data.frame(scores(nmds, "sites"))

#add metadata to the scores info
all.scores <- cbind(all.scores, com_meta)
all.scores$week <- as.factor(all.scores$week)

# make a new column that has unique names for week and treatment
all.scores <- all.scores %>% 
  unite(hull.id, treatment, week, remove = FALSE)

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
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "week") + aes(fill = factor(hull.id)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)



#ecoplate data just at each individual time point
#week1
eco_filt1 <- subset(ecodat, week == "1")
com2 <- eco_filt1
com2_meta <- subset(com2, select= c(sample, plant_block, treatment, week))
com2 <- subset(com2, select= -c(sample, plant_block, treatment, week))
com2[com2 < 0] <- 0
set.seed(123)
nmds2 <- metaMDS(com2,distance = "bray",
                          k = 2, 
                          trymax = 1000,
                          wascores = TRUE)

ordiplot(nmds2, display = "sites", type = "t")

sc <- as.data.frame(scores(nmds2$points[,1:2]))
sc$ID <- rownames(sc)
com2_meta$ID <- rownames(com2_meta)
com2_meta <- as.matrix(com2_meta)
sc <- merge(sc, com2_meta, by = c("ID"))
sc <- as.data.frame(sc)
sc$week <- as.factor(sc$week)
sc$treatment <- as.factor(sc$treatment)
plot(nmds2$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[sc$treatment],
     pch=19, main="Day 7")

ordiellipse(nmds2, groups = as.factor(sc$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)
legend("topleft", 
       legend=c("ant","beetle","fly", "lm", "cb", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=.5,
       bty = "n")

#week3
eco_filt1 <- subset(ecodat, week == "3")
com2 <- eco_filt1
com2_meta <- subset(com2, select= c(sample, plant_block, treatment, week))
com2 <- subset(com2, select= -c(sample, plant_block, treatment, week))
com2[com2 < 0] <- 0
set.seed(123)
nmds2 <- metaMDS(com2,distance = "bray",
                 k = 2, 
                 trymax = 1000,
                 wascores = TRUE)
sc <- as.data.frame(scores(nmds2$points[,1:2]))
sc$ID <- rownames(sc)
com2_meta$ID <- rownames(com2_meta)
com2_meta <- as.matrix(com2_meta)
sc <- merge(sc, com2_meta, by = c("ID"))
sc <- as.data.frame(sc)
sc$week <- as.factor(sc$week)
sc$treatment <- as.factor(sc$treatment)
plot(nmds2$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[sc$treatment],
     pch=19, main="Day 21")
ordiellipse(nmds2, groups = as.factor(sc$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)

legend("topleft", 
       legend=c("ant","beetle","fly", "lm", "cb", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=.5,
       bty = "n")

#week5
eco_filt1 <- subset(ecodat, week == "5")
com2 <- eco_filt1
com2_meta <- subset(com2, select= c(sample, plant_block, treatment, week))
com2 <- subset(com2, select= -c(sample, plant_block, treatment, week))
com2[com2 < 0] <- 0
set.seed(123)
nmds2 <- metaMDS(com2,distance = "bray",
                 k = 2, 
                 trymax = 1000,
                 wascores = TRUE)
sc <- as.data.frame(scores(nmds2$points[,1:2]))
sc$ID <- rownames(sc)
com2_meta$ID <- rownames(com2_meta)
com2_meta <- as.matrix(com2_meta)
sc <- merge(sc, com2_meta, by = c("ID"))
sc <- as.data.frame(sc)
sc$week <- as.factor(sc$week)
sc$treatment <- as.factor(sc$treatment)
plot(nmds2$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[sc$treatment],
     pch=19, main="Day 35")
ordiellipse(nmds2, groups = as.factor(sc$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)

legend("topleft", 
       legend=c("ant","beetle","fly", "lm", "cb", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=.5,
       bty = "n")

#week7
eco_filt1 <- subset(ecodat, week == "7")
com2 <- eco_filt1
com2_meta <- subset(com2, select= c(sample, plant_block, treatment, week))
com2 <- subset(com2, select= -c(sample, plant_block, treatment, week))
com2[com2 < 0] <- 0
set.seed(123)
nmds2 <- metaMDS(com2,distance = "bray",
                 k = 2, 
                 trymax = 1000,
                 wascores = TRUE)
sc <- as.data.frame(scores(nmds2$points[,1:2]))
sc$ID <- rownames(sc)
com2_meta$ID <- rownames(com2_meta)
com2_meta <- as.matrix(com2_meta)
sc <- merge(sc, com2_meta, by = c("ID"))
sc <- as.data.frame(sc)
sc$week <- as.factor(sc$week)
sc$treatment <- as.factor(sc$treatment)
plot(nmds2$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[sc$treatment],
     pch=19, main="Day 49")
ordiellipse(nmds2, groups = as.factor(sc$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)
legend("topleft", 
       legend=c("ant","beetle","fly", "lm", "cb", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=.5,
       bty = "n")


## Ecoplate Statistical Analysis
# beta disper to check dispersion
ecoall.bd <- betadisper(vegdist(com), com_meta$treatment)
ecoall.bd
boxplot(ecoall.bd)
anova(ecoall.bd)# p=0.0.3743

ecoall.bdweek <- betadisper(vegdist(com), com_meta$week)
ecoall.bdweek
boxplot(ecoall.bdweek)
anova(ecoall.bdweek)# p=4.422e-05 ***

# PERMANOVA for categorical variables (factors)
# set number of permutations
perm <- how(nperm = 999)
#specify a random variable (plant).
setBlocks(perm) <- with(all.scores, sample)
adonis <- adonis2(com ~ treatment * week, data=all.scores, permutations = perm, method = "bray")
adonis

# pairwise adonis
ecoall.ad.pw <- pairwise.adonis2(com ~ treatment * week, data = com_meta, strata = 'sample')
ecoall.ad.pw

#### GLMMs ####

#relevel the predictor variables to ACM
data_filt$pH <- as.numeric(data_filt$pH)
data_filt <- data_filt %>% arrange(desc(treatment))
data_filt$treatment <- relevel(data_filt$treatment, ref = "water")
decom_delta$treatment <- relevel(decom_delta$treatment, ref = "water")

#how does prey type influence microbial functions
#chitinase
summary(aov(cumulative_chit ~ plant_block, data=data_filt))#sig
m1 <- brm(cumulative_chit ~ treatment + (1|plant_block),data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
m1b <- brm(cumulative_chit ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(m1)
summary(m1b)

m1e <- brm(sqrt_chit ~ treatment*day,data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))

m1d <- brm(sqrt_chit ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(m1d)
summary(m1e)

mcmc_areas(m1d, regex_pars="b_",
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray") + ggtitle("sqrt chitinase ~ treatment")

mcmc_areas(m1e, regex_pars="b_",
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray") + ggtitle("sqrt chitinase ~ treatment")


#protease
summary(aov(cumulative_prot ~ plant_block, data=data_filt))#sig
m2 <- brm(cumulative_prot ~ treatment + (1|plant_block),data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 20))

m2b <- brm(cumulative_prot ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(m2)
summary(m2b)

#bacterial abundance (no difference between block so removed)
summary(aov(flow_cyt_livingevents_uL ~ plant_block, data=data_filt))#not sig
#m3 <- brm(flow_cyt_livingevents_uL ~ treatment+ (1| plant),data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
#          control = list(adapt_delta = 0.999, max_treedepth = 20))
m3b <- brm(flow_cyt_livingevents_uL ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(m3)#similar effects with and without the random effect
summary(m3b)

#ph (no difference between block so removed)
summary(aov(pH ~ plant_block, data=data_filt))#not sig
#m4 <- brm(pH ~ treatment+ (1| plant),data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
#         control = list(adapt_delta = 0.999, max_treedepth = 20))
m4b <- brm(pH ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(m4)#exact same effects with and without the random effect
summary(m4b)

#change in biomass
decom_delta$treatment <- as.factor(decom_delta$treatment)
decom_delta <- as.data.frame(decom_delta)
ggplot(data_filt, aes(x=plant_block, y=delta_prey_percent, color=plant_block)) + geom_boxplot()
summary(aov(delta_prey_percent ~ plant_block, data=decom_delta))
#no effect of plant block so removed random effect
m5 <- brm(delta_prey_percent ~ treatment,data=decom_delta, family=Gamma(link="log"), iter = 100000, chains = 4, cores = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 20))

##############
saveRDS(m1b, file = "brms_m1_chit.RDS")
saveRDS(m2b, file = "brms_m2_prot.RDS")
saveRDS(m3b, file = "brms_m3_propliving.RDS")
saveRDS(m4b, file = "brms_m4_ph.RDS")
saveRDS(m5, file = "brms_m5_biomass.RDS")

##############
m1b <- readRDS("brms_m1_chit.RDS")
m2b <- readRDS("brms_m2_prot.RDS")
m3b <- readRDS("brms_m3_propliving.RDS")
m4b <- readRDS("brms_m4_ph.RDS")
m5 <- readRDS("brms_m5_biomass.RDS")

##############
summary(m1b)
pp_check(m1b)

summary(m2b)
pp_check(m2b)

summary(m3b)
pp_check(m3b)

summary(m4b)
pp_check(m4b)

summary(m5)
pp_check(m5)

##############
mcmc_areas(m1b, pars=c("b_Intercept", "b_treatmentant","b_treatmentbeetle", "b_treatmentfly", "b_treatmentpos", "b_treatmentlm", "b_treatmentcb"),
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray") + ggtitle("cumulative chitinase ~ treatment")

mcmc_areas(m2b, pars=c("b_Intercept", "b_treatmentant","b_treatmentbeetle", "b_treatmentfly", "b_treatmentpos", "b_treatmentlm", "b_treatmentcb"),
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ ggtitle("cumulative protease ~ treatment ")

mcmc_areas(m3b, pars=c("b_Intercept", "b_treatmentant","b_treatmentbeetle", "b_treatmentfly", "b_treatmentpos", "b_treatmentlm", "b_treatmentcb"),
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ ggtitle("living events ~ treatment ")

mcmc_areas(m4b, pars=c("b_Intercept", "b_treatmentant","b_treatmentbeetle", "b_treatmentfly", "b_treatmentpos", "b_treatmentlm", "b_treatmentcb"),
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ ggtitle("ph ~ treatment 
                                                                         ")

mcmc_areas(m5, regex_pars=("b_"),
           prob = 0.90, 
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ ggtitle("delta biomass % ~ treatment")
m5
##############
pm2 <- ggeffects::ggpredict(m2, "treatment")
colnames(pm2)[1] <- "treatment"
colnames(pm2)[2] <- "cumulative_prot"

ggplot(data_filt, aes(x = treatment, y = cumulative_prot, fill=treatment, color=treatment)) +geom_boxplot(outlier.shape=NA) + geom_jitter(size=2)+
  theme(legend.position="none") + geom_point(data=pm2, aes(x=treatment, y=cumulative_prot), color="black", size=3) +
  geom_linerange(data=pm2,aes(ymin=conf.low, ymax=conf.high), size=1, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", y = "protease") + theme_classic()+ scale_fill_manual(values=c(color_palette))+
  scale_color_manual(values=c(color_palette))


pm1 <- ggeffects::ggpredict(m1, "treatment")
colnames(pm1)[1] <- "treatment"
colnames(pm1)[2] <- "cumulative_chit"

ggplot(data_filt, aes(x = treatment, y = cumulative_chit, color=treatment)) + geom_jitter(size=2)+
  theme(legend.position="none") + geom_point(data=pm1, aes(x=treatment, y=chitinase_uM_min), color="black", size=1) +
  geom_linerange(data=pm1,aes(ymin=conf.low, ymax=conf.high), size=1, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", y = "chitinase") + theme_classic()

posterior_m1 <- as.data.frame(m1)
posterior_m1_melt <- posterior_m1[,1:5]
posterior_m1_melt <- reshape2::melt(posterior_m1_melt)
ggplot(posterior_m1_melt, aes(x = value, y=variable,
                                    fill = stat(x < 0))) +
  stat_halfeye() +
  scale_fill_manual(values=c("#F7D95C", "gray"))+
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none") + xlab("Chit")

posterior_m2 <- as.data.frame(m2)
posterior_m2_melt <- posterior_m2[,1:5]
posterior_m2_melt <- reshape2::melt(posterior_m2_melt)
ggplot(posterior_m2_melt, aes(x = value, y=variable,
                              fill = stat(x < 0))) +
  stat_halfeye() +
  scale_fill_manual(values=c("#F7D95C", "gray"))+
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none") + xlab("Prot")

posterior_m3 <- as.data.frame(m3)
posterior_m3_melt <- posterior_m3[,1:5]
posterior_m3_melt <- reshape2::melt(posterior_m3_melt)
ggplot(posterior_m3_melt, aes(x = value, y=variable,
                              fill = stat(x < 0))) +
  stat_halfeye() +
  scale_fill_manual(values=c("#F7D95C", "gray"))+
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none") + xlab("Bac Abun")

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
saveRDS(Exp4.physeq1_raw, "Exp4.physeq1_raw.RDS")

#### DECONTAM PACKAGE for identifying contaminants####
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
## Read in data and metadata as a phyloseq object (non-rarified data, non-filtered)
#what percentage of the asv are only present in one sample before prev filt
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

#clean up meta??

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
saveRDS(Exp4.physeq2, "Exp4.physeq2.RDS")
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

#write.csv(md16s, "exp4_md16s.csv", row.names=TRUE)
md16s <- read.csv("exp4_md16s.csv", header=TRUE)
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

#### alpha diveristy GLMM
md16s$treatment <- as.factor(md16s$treatment)
md16s$treatment <- relevel(md16s$treatment, ref = "water")

glm1 <- brm(ef.16s ~ treatment + (1|plant), family = poisson, data=md16s, iter=10000)
mcmc_areas(glm1, regex_pars=("b"),
           prob = 0.9, 
           prob_outer = 1, 
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="black")+ggtitle("Effective number of ASVs ~ treatment +(1| plantID)")

glm2 <- brm(richness ~ treatment+ day +(1|plant_block), family = poisson, data=md16s, iter=10000)
mcmc_areas(glm2, regex_pars=("b_"),
           prob = 0.9, 
           prob_outer = 1, 
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="black")+ggtitle("ASV Richness ~ treatment+day +(1| plantblock)")

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
saveRDS(Exp4.physeq3, "Exp4.physeq3.RDS")
Exp4.physeq3 <- readRDS("Exp4.physeq3.RDS")

#no water
Exp4.nowater.physeq3 = subset_samples(Exp4.physeq3, treatment != "water") 

#no controls (pos or water)
Exp4.nocontrols.physeq3 = subset_samples(Exp4.nowater.physeq3, treatment != "pos") 

#split by time
Exp4.day7.physeq3 = subset_samples(Exp4.physeq3, day == "7") 
Exp4.day14.physeq3 = subset_samples(Exp4.physeq3, day == "14")
Exp4.day35.physeq3 = subset_samples(Exp4.physeq3, day == "35") 

#### BETA DIVERSITY ####
## Calculate weighted Unifrac distance and run NMDS
wu.dist.16s<- distance(Exp4.physeq3,"wUniFrac")

set.seed(123)
wu.nmds.16s <- metaMDS(wu.dist.16s,
                             k = 2, 
                             trymax = 1000,
                             wascores = TRUE)

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta_phys <- sample_data(Exp4.physeq3)
meta_phys$ID <- rownames(meta_phys)
meta_phys <- data.frame(meta_phys)
data_merge2 <- merge(data.scores2, meta_phys, by = c("ID"))
data_merge2$treatment <- as.factor(data_merge2$treatment)
data_merge2$day <- as.factor(data_merge2$day)

plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[data_merge2$treatment],
     pch=19)
legend("topright", 
       legend=c("ant", "beetle","fly","LM", "CB", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=1,
       bty = "n")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "ant", col = "#c24000")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "beetle", col = "#008388")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "fly", col = "#5d004f")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "lm", col = "#00446a")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "cb", col = "#d88a00")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "pos", col = "#006a50")
ordispider(wu.nmds.16s,groups = data_merge2$treatment, show.groups = "water", col = "slategray")

plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[data_merge2$treatment],
     pch=19, main="wUniFrac all time points")
ordiellipse(wu.nmds.16s, groups = as.factor(data_merge2$treatment),col = c(color_palette))
legend("topright", 
       legend=c("ant", "beetle","fly","LM", "CB", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=1,
       bty = "n")

#not plotting the day correctly
plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c("#9b5fe0","#16a4d8","#8bd346")[data_merge2$day],
     pch=19)
legend("topleft", 
       legend=c("day07", "day14","day35"),
       col= c("#9b5fe0","#16a4d8","#8bd346"),
       pch=19,
       cex=1,
       bty = "n")

ggplot(data=data_merge2, aes(x=MDS1, y=MDS2, color=treatment)) + geom_point()+
  theme_classic()

ggplot(data=data_merge2, aes(x=MDS1, y=MDS2, color=day)) + geom_point()+
  theme_classic()


#### day 7
## Calculate weighted Unifrac distance and run NMDS
wu.dist.16s<- distance(Exp4.day7.physeq3,"wUniFrac")

set.seed(123)
wu.nmds.16s <- metaMDS(wu.dist.16s,
                       k = 2, 
                       trymax = 1000,
                       wascores = TRUE)

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta_phys <- sample_data(Exp4.day7.physeq3)
meta_phys$ID <- rownames(meta_phys)
meta_phys <- data.frame(meta_phys)
data_merge2 <- merge(data.scores2, meta_phys, by = c("ID"))
data_merge2$treatment <- as.factor(data_merge2$treatment)
data_merge2$day <- as.factor(data_merge2$day)

plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[data_merge2$treatment],
     pch=19, main="wUniFrac Day 7")
ordiellipse(wu.nmds.16s, groups = as.factor(data_merge2$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)
legend("bottomright", 
       legend=c("ant", "beetle","fly","LM", "CB", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=1,
       bty = "n")



#### day 14
wu.dist.16s<- distance(Exp4.day14.physeq3,"wUniFrac")

set.seed(123)
wu.nmds.16s <- metaMDS(wu.dist.16s,
                       k = 2, 
                       trymax = 1000,
                       wascores = TRUE)

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta_phys <- sample_data(Exp4.day14.physeq3)
meta_phys$ID <- rownames(meta_phys)
meta_phys <- data.frame(meta_phys)
data_merge2 <- merge(data.scores2, meta_phys, by = c("ID"))
data_merge2$treatment <- as.factor(data_merge2$treatment)
data_merge2$day <- as.factor(data_merge2$day)

plot(wu.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[data_merge2$treatment],
     pch=19, main="wUniFrac Day 14")
ordiellipse(wu.nmds.16s, groups = as.factor(data_merge2$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)
legend("bottomleft", 
       legend=c("ant", "beetle","fly","LM", "CB", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=1,
       bty = "n")



#### day 35
wu35.dist.16s<- distance(Exp4.day35.physeq3,"wUniFrac")

set.seed(123)
wu35.nmds.16s <- metaMDS(wu35.dist.16s,
                       k = 2, 
                       trymax = 1000,
                       wascores = TRUE)

## Plot NMDS
data.scores35 <- as.data.frame(scores(wu35.nmds.16s$points[,1:2]))
data.scores35$ID <- rownames(data.scores35)
meta_phys35 <- sample_data(Exp4.day35.physeq3)
meta_phys35$ID <- rownames(meta_phys35)
meta_phys35 <- data.frame(meta_phys35)
data_merge35 <- merge(data.scores35, meta_phys35, by = c("ID"))
data_merge35$treatment <- as.factor(data_merge35$treatment)
data_merge35$day <- as.factor(data_merge35$day)

prey_order <- c("ant", "beetle", "fly", "lm", "cb", "pos", "water")

# Use the factor function to set the order of levels
data_merge35$treatment <- factor(data_merge35$treatment, levels = prey_order)

plot(wu35.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[data_merge35$treatment],
     pch=19, main="wUniFrac Day 35")
ordiellipse(wu35.nmds.16s, groups = as.factor(data_merge35$treatment),kind = "se",col = c(color_palette), label = TRUE, conf=.95)
legend("topleft", 
       legend=c("ant", "beetle","fly","LM", "CB", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=1,
       bty = "n")

plot(wu35.nmds.16s$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c(color_palette)[data_merge35$treatment],
     pch=19, main="wUniFrac Day 35")
legend("topleft", 
       legend=c("ant", "beetle","fly","LM", "CB", "pos", "water"),
       col= c(color_palette),
       pch=19,
       cex=1,
       bty = "n")

ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "ant", col = "#c24000")
ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "beetle", col = "#008388")
ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "fly", col = "#5d004f")
ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "lm", col = "#00446a")
ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "cb", col = "#d88a00")
ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "pos", col = "#006a50")
ordispider(wu35.nmds.16s,groups = data_merge35$treatment, show.groups = "water", col = "slategray")

#### Beta Diversity Analysis
# beta disper to check dispersion between treatments day 35
wu.day35.bd <- betadisper(wu35.dist.16s, data_merge35$treatment)
wu.day35.bd
boxplot(wu.day35.bd)
anova(wu.day35.bd)# p=0.5574, no significant differences in the dispersion within each group between groups...


# PERMANOVA for categorical variables (factors)
#day 35
# set number of permutations
perm <- how(nperm = 999)
#specify a random variable (plant_block).
setBlocks(perm) <- with(data_merge35, plant_block)
ad35.16s.treat <- adonis2(wu35.dist.16s ~ data_merge35$treatment, permutations = perm, by="margin")
ad35.16s.treat #0.18

# pairwise adonis
wu.ad.pw <- pairwise.adonis2(wu.dist.16s ~ treatment * day, data = data_merge2, strata = 'plant_block')
wu.ad.pw


#MANTELL TESTS: community composition correlated with functional differences?


#### ANCOMBC ####



#### RELATIVE ABUNDANCE PLOTS ####
tax <- data.frame(tax_table(Exp4.physeq3))
samp <- data.frame(sample_data(Exp4.physeq3))
otu <- data.frame(otu_table(Exp4.physeq3))
rownames(otu)==rownames(tax)

#some different color palettes to try
library(viridis)
color_palette21 <- viridis_pal()(21)

# Plotting the relative abundance at Genus level
otu.g <- data.frame(Genus=tax$Genus,otu)
otu.g$Genus[is.na(otu.g$Genus)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Genus, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])

otu.gaop <- otu.gao[c(1:5,8:22),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
other <- colSums(otu.gao[c(6,7,23:nrow(otu.gao)),])#put the all the other genera in a "other" category
otu.gaopo <- rbind(otu.gaop, other)
data_percentage <- apply(otu.gaopo, 2, function(x){x*100/sum(x,na.rm=T)})#calculate relative abundance in each sample
barplot(data_percentage,col=color_palette21,legend.text=T,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.05, -0.05), cex=0.4))

#plot using ggplot
otu.gaopo_gg <- as.data.frame(otu.gaopo)
otu.gaopo_ggp <- adorn_percentages(otu.gaopo_gg,, 1:ncol(otu.gaopo_gg), denominator="col")#number of columns (samples), calculate RA
other <- otu.gaopo_ggp[21,]
otu.gaopo_ggp <- otu.gaopo_ggp[-21,]
otu.gaopo_ggp <- otu.gaopo_ggp[order(rowSums(otu.gaopo_ggp),decreasing = F),]
otu.gaopo_ggp <- rbind(other, otu.gaopo_ggp)
otu.gaopo_ggp[ "genus" ] <- rownames(otu.gaopo_ggp)
otu.gaopo_ggp$genus <- factor(otu.gaopo_ggp$genus, levels = otu.gaopo_ggp$genus)
otu.gaopo_ggp.rmelt<- reshape2::melt(otu.gaopo_ggp, id.vars="genus", value.name="Relative_Abundance", variable.name="Sample")
df3 <- left_join(otu.gaopo_ggp.rmelt, samp, by=c("Sample"="X"))
df3$day <- as.factor(df3$day)
ggplot(df3, aes(x=day, y=Relative_Abundance, fill = genus)) + 
  geom_bar( position = "fill", stat = "identity", width=1) + 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values = color_palette21) + facet_wrap(~treatment, nrow=2)+
  guides(fill = guide_legend(ncol = 1)) 

### make individual relative abundance plots for each treatment
Exp4.ant.physeq3 = subset_samples(Exp4.physeq3, treatment == "ant") 
Exp4.beetle.physeq3 = subset_samples(Exp4.physeq3, treatment == "beetle") 
Exp4.fly.physeq3 = subset_samples(Exp4.physeq3, treatment == "fly") 
Exp4.cb.physeq3 = subset_samples(Exp4.physeq3, treatment == "cb") 
Exp4.lm.physeq3 = subset_samples(Exp4.physeq3, treatment == "lm") 
Exp4.pos.physeq3 = subset_samples(Exp4.physeq3, treatment == "pos") 
Exp4.water.physeq3 = subset_samples(Exp4.physeq3, treatment == "water") 
jewel_tone_palette <- c("#000000","#a9a19c","#493829","#117a65","#bdd09f", "#0c2461", "#d35400",
                        "#1f8a13","#8f3b1b","#855723","#e67e22","#740058",
                        "palegoldenrod","#3498db","#f39c12", "#28b463","#af7ac5",
                        "#1f618d","#d35400","#154360","darkmagenta")
###ant
tax <- data.frame(tax_table(Exp4.ant.physeq3))
samp <- data.frame(sample_data(Exp4.ant.physeq3))
otu <- data.frame(otu_table(Exp4.ant.physeq3))
rownames(otu)==rownames(tax)

# Plotting the relative abundance at family level
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]

otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])
otu.gaop <- otu.gao[c(1:16,18:21),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
other <- colSums(otu.gao[c(17,23:nrow(otu.gao)),])#put the all the other genera in a "other" category
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Ant")


###beetle
tax <- data.frame(tax_table(Exp4.beetle.physeq3))
samp <- data.frame(sample_data(Exp4.beetle.physeq3))
otu <- data.frame(otu_table(Exp4.beetle.physeq3))
rownames(otu)==rownames(tax)
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Beetle")

###fly
tax <- data.frame(tax_table(Exp4.fly.physeq3))
samp <- data.frame(sample_data(Exp4.fly.physeq3))
otu <- data.frame(otu_table(Exp4.fly.physeq3))
rownames(otu)==rownames(tax)
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Fly")

###Les Mosses
tax <- data.frame(tax_table(Exp4.lm.physeq3))
samp <- data.frame(sample_data(Exp4.lm.physeq3))
otu <- data.frame(otu_table(Exp4.lm.physeq3))
rownames(otu)==rownames(tax)
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Les Mosses")

###Champ Buet
tax <- data.frame(tax_table(Exp4.cb.physeq3))
samp <- data.frame(sample_data(Exp4.cb.physeq3))
otu <- data.frame(otu_table(Exp4.cb.physeq3))
rownames(otu)==rownames(tax)
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])
otu.gaop <- otu.gao[c(1:13, 15:21),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
other <- colSums(otu.gao[c(14,21:nrow(otu.gao)),])#put the all the other genera in a "other" category
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Champ Buet")


###Positive control
tax <- data.frame(tax_table(Exp4.pos.physeq3))
samp <- data.frame(sample_data(Exp4.pos.physeq3))
otu <- data.frame(otu_table(Exp4.pos.physeq3))
rownames(otu)==rownames(tax)
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Positive Control")


###water control
tax <- data.frame(tax_table(Exp4.water.physeq3))
samp <- data.frame(sample_data(Exp4.water.physeq3))
otu <- data.frame(otu_table(Exp4.water.physeq3))
rownames(otu)==rownames(tax)
otu.g <- data.frame(Family=tax$Family,otu)
otu.g$Family[is.na(otu.g$Family)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Family, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.ga <- otu.ga[rowSums(otu.ga) != 0, ]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])
otu.gaop <- otu.gao[c(1:13, 15:21),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
other <- colSums(otu.gao[c(14,22:nrow(otu.gao)),])#put the all the other genera in a "other" category
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
df3$day <- factor(df3$day, levels = c("7", "14", "35"))
ggplot(df3, aes(x = reorder(sample_id, as.numeric(factor(day))), y = Relative_Abundance, fill = family)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_fill_manual(values = jewel_tone_palette) +
  guides(fill = guide_legend(ncol = 1))+ ggtitle("Water")





#### Protist ####
protist <- read.csv("protist_matrix.csv", header=TRUE)
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
ggplot(data=prot_meta, aes(x=day, y=ef.prot, group=treatment, color=treatment)) +
  geom_jitter() + ylab("Effective Number of Species") + theme_classic()+
  geom_smooth(se=FALSE)

ggplot(data=prot_meta, aes(x=day, y=richness, group=treatment, color=treatment)) +
  geom_jitter() + ylab("Richness") + theme_classic()+
  geom_smooth(se=FALSE)

ggplot(data=prot_meta, aes(x=treatment, y=richness, group=treatment, color=treatment)) +
  geom_boxplot() + ylab("Richness") + theme_classic()+geom_jitter()+
  facet_wrap(~day, nrow=2)

#### BETA DIVERSITY ####
## Calculate Jaccard distance and run NMDS
prot_com.filt <- prot_com[rowSums(prot_com) != 0, ]
prot_meta.filt <- prot_meta[match(rownames(prot_com.filt), rownames(prot_meta)), ]
prot.dist<- vegdist(prot_com.filt, method="jaccard", binary=TRUE)
set.seed(123)
prot.nmds <- metaMDS(prot.dist,
                       k = 2, 
                       trymax = 1000,
                       wascores = TRUE)

## Plot NMDS
ds <- as.data.frame(scores(prot.nmds$points[,1:2]))
ds$ID <- rownames(ds)
prot_meta.filt$ID <- rownames(prot_meta.filt)
data_merge2 <- merge(ds, prot_meta.filt, by = c("ID"))
data_merge2$treatment <- as.factor(data_merge2$treatment)
data_merge2$day <- as.factor(data_merge2$day)

ggplot(data_merge2, aes(x=MDS1, y=MDS2, color=treatment))+
  geom_point(position = position_jitterdodge(jitter.width = 8000, dodge.width = 5))+
  theme_classic()


## Calculate BC distance and run NMDS
prot.dist2<- vegdist(prot_com.filt, method="bray")
set.seed(123)
prot.nmds2 <- metaMDS(prot.dist2,
                     k = 2, 
                     trymax = 1000,
                     wascores = TRUE)

## Plot NMDS
ds <- as.data.frame(scores(prot.nmds2$points[,1:2]))
ds$ID <- rownames(ds)
prot_meta.filt$ID <- rownames(prot_meta.filt)
data_merge2 <- merge(ds, prot_meta.filt, by = c("ID"))
data_merge2$treatment <- as.factor(data_merge2$treatment)
data_merge2$day <- as.factor(data_merge2$day)

ggplot(data_merge2, aes(x=MDS1, y=MDS2, color=treatment))+
  geom_point(position = position_jitterdodge(jitter.width = 8000, dodge.width = 5))+
  theme_classic()

#try PCoA
pcoa <- wcmdscale(d = prot.dist, eig = TRUE)
a <- as.data.frame(pcoa$points)
ab <- cbind(a, prot_meta.filt)
ab$day <- as.factor(ab$day)
ggplot(data = ab,
       aes(x = Dim1, y = Dim2, color=treatment)) + geom_point() + theme_bw()
ggplot(data = ab,
       aes(x = Dim1, y = Dim2, color=day)) + geom_point() + theme_bw()
ordiplot(pcoa, display = 'sites', type = 'text') 


#### 16S METACODER DIFFERENTIAL ABUNDANCE HEAT TREES ####




#### 16S METACODER DIFFERENTIAL ABUNDANCE HEAT TREES ####
#using cleaned up and rarified data
#subset to week 1 and week 8
asv.35 <- asv16s.r[, grep(pattern="35$", colnames(asv16s.r))]
tax.35 <- subset(tax.16s2.r, row.names(tax.16s2.r) %in% row.names(asv.35)) #filter tax table to match asvs
asv.tax.35 <- cbind(tax.35,asv.35)
meta.35 <- subset(meta.r, row.names(meta.r) %in% colnames(asv.35)) #filter tax table to match asvs
meta.35 <- rownames_to_column(meta.35, var = "sampleID")

asv.tax.35 <- subset(asv.tax.35, select = -Confidence)
asv.tax.35 <- subset(asv.tax.35, Taxon != "Unassigned")

## the taxonomy section is hard to handle because of formatting, we can use the `taxa` package to parse this and process the abundance data too
obj <- parse_tax_data(asv.tax.35,
                      class_cols = "Taxon",
                      class_sep = ";",
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

#Calculate proportions
obj$data$tax_data <- metacoder::calc_obs_props(obj, "tax_data",
                                               cols = meta.35$sampleID)


# now we need to calculate abundances based on taxon not ASV
obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "tax_data",
                                                  cols = meta.35$sampleID)

#16S heatmap matrix by treatment
obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
                                                 cols = meta.35$sampleID, # What columns of sample data to use
                                                 groups = meta.35$treatment) # What category each sample is assigned to

obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

#lets look at the p-values and see if there is any significance
range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# [1] 0.6874437 1.00000000

## Focus only on significant taxa
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

set.seed(1115)
#this takes a while
diff_heattree <- metacoder::heat_tree_matrix(obj, data = "diff_table",
                                                   node_size = n_obs, 
                                                  node_label = taxon_names,
                                                  node_color = log2_median_ratio,
                                                  node_color_range = diverging_palette(),
                                                  node_color_trans = "linear", 
                                                  node_color_interval = c(-3, 3), 
                                                  edge_color_interval = c(-3, 3), 
                                                  node_size_axis_label = "Number of ASVs",
                                                  node_color_axis_label = "Log2 ratio median proportions",
                                                  layout = "davidson-harel", 
                                                  initial_layout = "reingold-tilford")

#saveRDS(diff_heattree, file = "diff_heattree.RDS")
diff_heattree <- readRDS("diff_heattree.RDS")
print(diff_heattree) ## Show taxonomic heat tree 



