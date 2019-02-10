#####Miombo Soils - Pilot Depth Profiles#####

library(dplyr)
library(tidyr)
library(ggplot2)
library(agricolae)

setwd("C:/Users/mthom/Documents/Miombo/thesisrevisions_manuscripts/Chp1_NatGeosubmission/Data")


#reading in soils dataset from pilot sites to 2 m
soils1 <- read.csv("TZSfullpilot_depthsamples.csv")
names(soils1)

#extract agriculture and old-growth forest sites for contextual enhancement to chronosequence sites.  No aggregation by depth or class.
soils1_fa<- soils1 %>%
  filter(LANDCOVER == "FOR.OG" | LANDCOVER == "AGR") %>%
  select(-NH4mgkg_avg, -NH4mgkg_stderr, -NH4kgha_avg, -NH4kgha_stderr, -NO3mgkg_avg, -NO3mgkg_stderr, -NO3kgha_avg, -NO3kgha_stderr, -NO3toNH4, -TotNkgha_avg, -TotNkgha_stderr) %>%
  droplevels()
write.csv(soils1_fa, "soils1_fa.csv")

#including regrowth forest sites into the mix...
soils1_far<- soils1 %>%
  filter(LANDCOVER == "FOR.OG" | LANDCOVER == "AGR" | LANDCOVER == "FOR.RG") %>%
  select(-NH4mgkg_avg, -NH4mgkg_stderr, -NH4kgha_avg, -NH4kgha_stderr, -NO3mgkg_avg, -NO3mgkg_stderr, -NO3kgha_avg, -NO3kgha_stderr, -NO3toNH4, -TotNkgha_avg, -TotNkgha_stderr) %>%
  droplevels()
write.csv(soils1_far, "soils1_far.csv")


#aggregation of agriculture and old-growth forest soil site properties by land-use type and depth:
soils1_fa_sumstats <- soils1_fa %>%
  group_by(LANDCOVER, DEPTHM) %>%
  summarize(sand = mean(SAND),
            sand_se = sd(SAND)/sqrt(length(SAND)),
            siltclay = mean(SILT+CLAY),
            siltclay_se = sd(SILT+CLAY)/sqrt(length(SILT+CLAY)),
            PH = mean(pH),
            PH_se = sd(pH)/sqrt(length(pH)),
            C_gkg = mean(Cgkg),
            C_gkg_se = sd(Cgkg)/sqrt(length(Cgkg)),
            N_gkg = mean(Ngkg),
            N_gkg_se = sd(Ngkg)/sqrt(length(Ngkg)),
            SOC = mean(C_Mgha), 
            SOC_se = sd(C_Mgha)/sqrt(length(C_Mgha)),
            TN = mean(N_kgha),
            TN_sd = sd(N_kgha)/sqrt(length(N_kgha)),
            CNratio = mean(CtoN),
            CNratio_sd =sd(CtoN)/sqrt(length(CtoN)),
            D15N = mean(d15N),
            D15N_se = sd(d15N)/sqrt(length(d15N)))
write.csv(soils1_fa_sumstats, "soils1_fa_sumstats.csv")

soils1_far_sumstats <- soils1_far %>%
  group_by(LANDCOVER, DEPTHM) %>%
  summarize(sand = mean(SAND),
            sand_se = sd(SAND)/sqrt(length(SAND)),
            siltclay = mean(SILT+CLAY),
            siltclay_se = sd(SILT+CLAY)/sqrt(length(SILT+CLAY)),
            PH = mean(pH),
            PH_se = sd(pH)/sqrt(length(pH)),
            C_gkg = mean(Cgkg),
            C_gkg_se = sd(Cgkg)/sqrt(length(Cgkg)),
            N_gkg = mean(Ngkg),
            N_gkg_se = sd(Ngkg)/sqrt(length(Ngkg)),
            SOC = mean(C_Mgha), 
            SOC_se = sd(C_Mgha)/sqrt(length(C_Mgha)),
            TN = mean(N_kgha),
            TN_sd = sd(N_kgha)/sqrt(length(N_kgha)),
            CNratio = mean(CtoN),
            CNratio_sd =sd(CtoN)/sqrt(length(CtoN)),
            D15N = mean(d15N),
            D15N_se = sd(d15N)/sqrt(length(d15N)))
write.csv(soils1_far_sumstats, "soils1_far_sumstats.csv")
names(soils1_far)

#summary statistics for soils of all 
mean_pctsand = mean(soils1_far$SAND)
sd_pctsand = mean(soils1_far$SAND)/length(soils1_far$SAND)
mean_pctsilt = mean(soils1_far$SILT)
sd_pctsilt = mean(soils1_far$SILT)/length(soils1_far$SILT)
mean_pctclay = mean(soils1_far$CLAY)
sd_pctclay = mean(soils1_far$CLAY)/length(soils1_far$CLAY)
pilotcores_sumstats <- as.list(mean_pctsand, sd_pctsand, mean_pctsilt, sd_pctsilt, mean_pctclay, sd_pctclay)
names(pilotcores_sumstats) <- c("PCTSAND", "SE_SAND", "PCTSILT", "SE_SILT", "PCTCLAY", "SE_CLAY")
mean_pH = mean(soils1_far$pH)
se_pH = sd(soils1_far$pH)/length(soils1_far$pH)


#check contents of summary soils
hist(soils1_fa$C_Mgha)
hist(log(soils1_fa$C_Mgha))
shapiro.test(log(soils1_fa$C_Mgha))
#log-transformed organic C stocks meet normality assumption..

hist(soils1_fa$Cgkg)
hist(log(soils1_fa$Cgkg+1))
hist(soils1_fa$Cgkg^(1/6))
shapiro.test(soils1_fa$Cgkg^(1/6))
#soil C concentrations raised to the 1/6 power are close to a normal distribution (best as I can find).
hist(soils1_far$Cgkg)
hist(log(soils1_far$Cgkg+1))
hist(soils1_far$Cgkg^(1/6))
shapiro.test(soils1_far$Cgkg^(1/6))


hist(soils1_fa$N_kgha)
hist(log(soils1_fa$N_kgha))
shapiro.test(log(soils1_fa$N_kgha+1))
#log-transformed total N meets normality assumption.

hist(soils1_fa$Ngkg)
hist(soils1_fa$Ngkg^(1/4))
shapiro.test(soils1_fa$Ngkg^(1/4))
#soil N concentrations raised to the 1/4 power meet normality assumptions.
hist(soils1_far$Ngkg)
hist(soils1_far$Ngkg^(1/4))
shapiro.test(soils1_far$Ngkg^(1/4))
#soil N concentrations raised to 1/4 power meet normality assumptions.

hist(soils1_far$CtoN)
hist(log(soils1_far$CtoN))
hist(log(soils1_far$CtoN)^(1/2))
shapiro.test(log(soils1_far$CtoN))
shapiro.test(log(soils1_far$CtoN)^(1/2))

#####ANOVA: site classs effects on C and N

#dataset1 = mature forest and ag pilot soil cores to 2 m
dataset1 <- soils1_fa
dataset1$DEPTHR <- factor(dataset1$DEPTHR)
dataset2 <- soils1_far
dataset2$DEPTHR <- factor(dataset2$DEPTHR)
dataset3 <- soils1_far2

names(dataset1)
Cstocks_class.lm <- lm(log(C_Mgha)~LANDCOVER+DEPTHR, dataset1)
summary(Cstocks_class.lm)
anova(Cstocks_class.lm)
#Cstocks_class.aov <- aov(lm(log(C_Mgha)~LANDCOVER+DEPTHR, dataset1))
#HSD.test(LANDCOVER ~ DEPTHR, data=lm(l) group=TRUE)
#plot(TukeyHSD(aov(Cstocks_class.lm)))
#wilcox.test(C_Mgha~LANDCOVER, dataset1)
#SUMMARY:
#SIGNIFICANT VARIABILIY OF C stocks by depth, but not land cover class.

#FOR.OG and AG: analysis of land cover and depth effects on C concetrations
Cconcs_classdepth.lm <- lm(Cgkg^(1/6)~LANDCOVER+DEPTHR, dataset1)
summary(Cconcs_classdepth.lm)
anova(Cconcs_classdepth.lm)
Cconcs_classdepth <- aov(Cgkg^(1/6)~LANDCOVER+DEPTHR, dataset1)
print(HSD.test(Cconcs_classdepth, "DEPTHR", group=TRUE))
anova(Cconcs_class.lm)
#Significant effects of depth but not land cover class, on C concentration.

#FOR.RG, FOR.OG and AG:
Cconcs2_classdepth.lm <- lm(Cgkg^(1/6)~LANDCOVER+DEPTHR, dataset2)
summary(Cconcs2_classdepth.lm)
anova(Cconcs2_classdepth.lm)
Cconcs2_classdepth <- aov(Cgkg^(1/6)~LANDCOVER+DEPTHR, dataset2)
print(HSD.test(Cconcs2_classdepth, "DEPTHR", group=TRUE))
print(HSD.test(Cconcs2_classdepth, "LANDCOVER", group=TRUE))
anova(Cconc2s_class.lm)

Nconcs_classdepth.lm <- lm(Ngkg^(1/4)~LANDCOVER+DEPTHR, dataset1)
summary(Nconcs_classdepth.lm)
anova(Nconcs_classdepth.lm)
Nconcs_classdepth <- aov(Ngkg^(1/4)~LANDCOVER+DEPTHR, dataset1)
print(HSD.test(Nconcs_classdepth, "DEPTHR", group=TRUE))
print(HSD.test(Nconcs_classdepth, "LANDCOVER", group=TRUE))
#Significant effects of depth and land cover class on N - with higher N concs in agriculture than old forests.

Nconcs2_classdepth.lm <- lm(Ngkg^(1/4)~LANDCOVER+DEPTHR, dataset2)
summary(Nconcs2_classdepth.lm)
anova(Nconcs2_classdepth.lm)
Nconcs2_classdepth <- aov(Ngkg^(1/4)~LANDCOVER+DEPTHR, dataset2)
print(HSD.test(Nconcs2_classdepth, "DEPTHR", group=TRUE))
print(HSD.test(Nconcs2_classdepth, "LANDCOVER", group=TRUE))
kruskal.test(Ngkg~LANDCOVER, dataset2)

#significant effect of land cover class on texture, but not depth.
kruskal.test(SAND~LANDCOVER, dataset1)
kruskal.test(SAND~DEPTHR, dataset1)
##Ag sites not abandoned - finer soil texture, more fertile.  

#summary stats of %sand across pilot soil core land cover groups
kruskal.test(SAND~LANDCOVER, dataset2)
soils1_far_lcgroup <- soils1_far %>%
  group_by(LANDCOVER) %>%
  summarize(pctsand = mean(SAND),
            se = sd(SAND)/length(SAND), 
            Cpct = mean(Cgkg),
            Cpct_se = sd(Cgkg)/length(Cgkg),
            Npct = mean(Ngkg),
            Npct_se = sd(Ngkg)/length(Ngkg))
write.csv(soils1_far_lcgroup, "soils1_far_lcgroup.csv")

kruskal.test(SAND~DEPTHR, dataset2)


#try kruskal-wallis test: effects of land cover class and depth on C:N
kruskal.test(CtoN~LANDCOVER, dataset2)
kruskal.test(CtoN~DEPTHR, dataset2)
###RESULT:  CtoN ratio does not differ significantly between ag, regrowth and mature forests; but does decrease significantly with depth.
