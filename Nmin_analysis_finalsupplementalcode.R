#####Miombo N Cycling with Regrowth: Data analysis #####

#####set working directory#####
getwd()
setwd("E:/TP_thesisprocessing_transfer/TP/Nanalysis")
#if laptop

setwd("D:/TP_thesisprocessing_transfer/TP/Nanalysis")
#####R packages used#####
library(plyr)
library(dplyr)
library(ggplot2)
list.files()
library(car)
library(nortest)
library(AID)
library(agricolae)

#####Input data files#####
#NMINF1 - Nmin data; LEAFCHEN - canopy leaf data; CANOPY - scaled canopy data

NMIN <- read.csv("NMINF1_correctSOCTN_nodepthsums_newSOCTN_newcompile.csv")
LEAF <- read.csv("LEAFCHEM2.csv")
CANOPY <- read.csv("CANOPY.csv")


####Foliar and canopy data##### 
names(LEAF)
names(CANOPY)
#LAIWETEDRYMIN = difference, wet-dry season LAI

par(mfrow=c(2,4))
hist(LEAF$d15N, main=paste("d15N"), xlab = c("Foliar d15N"), ylab = c("Num. samples"))
hist((LEAF$PCT_N*10), main=paste("Foliar N conc."), xlab = c("Foliar N g kg-1"), ylab = c("Num. samples"))
hist(CANOPY$LAIWETEDRYMIN, main=paste("Wet-DryLAI"), xlab = c("Wet-Dry Season LAI"), ylab = c("Num. sites"))
hist(CANOPY$CanopyN_kgha, main=paste("N Demand of Leaf Prod"), xlab = c("N kg ha-1"), ylab = c ("Num sites"))

hist(log(LEAF$d15N+3), main=paste("d15N"), xlab = c("(log Foliar d15N+3)"), ylab = c("Num. samples"))
hist(log(LEAF$PCT_N*10), main=paste("log Foliar N"), xlab = c("log Foliar N g kg-1"), ylab = c("Num. samples"))
hist(CANOPY$LAIWETEDRYMIN, main=paste("Wet-DryLAI"), xlab = c("Wet-Dry Season LAI"), ylab = c("Num. sites"))
hist(CANOPY$CanopyN_kgha, main=paste("N Demand of Leaf Prod"), xlab = c("N kg ha-1"), ylab = c ("Num sites"))


shapiro.test(log(LEAF$d15N+3))
shapiro.test((LEAF$PCT_N*10))
shapiro.test(CANOPY$LAIWETEDRYMIN)
shapiro.test(CANOPY$CanopyN_kgha)

#non-parametric tests for leaf d15N and PctN, not transformable to normal.
#raw data distributions not different than normal (shapiro-wilk p > 0.05) for canopy LAI and N_kgha, ok for parametric tests.


######Data tables: foliar and canopy data#####

leafchem_treeobs <- LEAF %>%
  group_by(SITE.TREE.N, SITE, CLASS, LABEL, TREE, SITE.TREE, CLASS.TREE) %>%
  summarise(Nleaf = length(PCT_N),
            PctC = mean(PCT_C),
            PctC_se = sd(PCT_C)/sqrt(length(PCT_N)),
            PctN = mean(PCT_N),
            PctN_se = sd(PCT_N)/sqrt(length(PCT_N)),
            CNx = mean(CN),
            CNx_se = sd(CN)/sqrt(length(PCT_N)),
            d15Nx = mean(d15N),
            d15N_se = sd(d15N)/sqrt(length(PCT_N)),
            SLA_gm2 = mean(SLAgm2),
            SLA_se = sd(SLAgm2, na.rm=T)/sqrt(length(PCT_N)))
leafchem_treeobs
write.csv(leafchem_treeobs, "leafchem_treeobs.csv")

leafchem_treeobs_Nbudget2 <- LEAF %>%
  group_by(SITE.TREE, CLASS.TREE) %>%
  summarise(Nleaf = length(PCT_N),
            PctC = mean(PCT_C),
            PctC_se = sd(PCT_C)/sqrt(length(PCT_N)),
            PctN = mean(PCT_N),
            PctN_se = sd(PCT_N)/sqrt(length(PCT_N)),
            CNx = mean(CN),
            CNx_se = sd(CN)/sqrt(length(PCT_N)),
            d15Nx = mean(d15N),
            d15N_se = sd(d15N)/sqrt(length(PCT_N)),
            SLA_gm2 = mean(SLAgm2),
            SLA_se = sd(SLAgm2, na.rm=T)/sqrt(length(PCT_N)))
write.csv(leafchem_treeobs_Nbudget2, "leafchem_treeobs_Nbudget2.csv")  

#TABLE SX Part I.  Leaf chemistry across classes (trees as observational units) - Foliar N, Foliar C:N.  Canopy data in following table.
names(leafchem_treeobs)
leafchem_summary_CLASS <- leafchem_treeobs %>%
  group_by(CLASS) %>%
  summarise(Nleaf = length(PctN),
            PCTC = mean(PctC),
            PCTC_se = sd(PctC)/sqrt(length(PctC)),
            PCTN = mean(PctN),
            PCTN_se = sd(PctN)/sqrt(length(PctN)),
            CN = mean(CNx),
            CN_se = sd(CNx)/sqrt(length(CNx)),
            d15N = mean(d15Nx),
            d15N_se = sd(d15Nx)/sqrt(length(d15Nx)),
            SLAgm2 = mean(SLA_gm2),
            SLAse = sd(SLA_gm2, na.rm=T)/sqrt(length(SLA_gm2)))
leafchem_summary_CLASS
write.csv(leafchem_summary_CLASS, "leafchem_summary_CLASS.csv")


#TABLE SX Part II.  Canopy LAI and N demand classes (sites as observational units) - Foliar data in preceding table.
names(CANOPY)
canopy_summary <- CANOPY %>%
  group_by(CLASS) %>%
  summarize(Nclass = length(LAIWETEDRYMIN),
            dLAI = mean(LAIWETEDRYMIN),
            dLAIse = sd(LAIWETEDRYMIN)/sqrt(length(LAIWETEDRYMIN)),
            CanopyCkghax = mean(CanopyC_kgha),
            CanopyCkghase = sd(CanopyC_kgha)/sqrt(length(LAIWETEDRYMIN)),
            CANOPYNkghax = mean(CanopyN_kgha),
            CANOPYNkghase = sd(CanopyN_kgha)/sqrt(length(LAIWETEDRYMIN))
            )
canopy_summary
write.csv(canopy_summary, "canopy_summary.csv")

#TABLE SX.  Leaf chemistry across trees by class.  Complete table.
#Leaf chemistry by tree at each class (trees as observational units)
leafchem_summary_CLASSTREE <- LEAF %>%
  group_by(CLASS.TREE) %>%
  summarise(Nleaf = length(PCT_N),
            PctC = mean(PCT_C),
            PctC_se = sd(PCT_C)/sqrt(length(PCT_N)),
            PctN = mean(PCT_N),
            PctN_se = sd(PCT_N)/sqrt(length(PCT_N)),
            CNx = mean(CN),
            CNx_se = sd(CN)/sqrt(length(PCT_N)),
            d15Nx = mean(d15N),
            d15N_se = sd(d15N)/sqrt(length(PCT_N)),
            SLA_gm2 = mean(SLAgm2),
            SLA_se = sd(SLAgm2)/sqrt(length(PCT_N)))
leafchem_summary_CLASSTREE
write.csv(leafchem_summary_CLASSTREE, "leafchem_summary_CLASSTREE.csv")

#for writing: leaf chemistry differences by tree species:
names(leafchem_treeobs)
leafchem_summary_TREE <- leafchem_treeobs %>%
  group_by(TREE) %>%
  summarise(Nleaf = length(PctN),
            PCTC = mean(PctC),
            PCTC_se = sd(PctC)/sqrt(length(PctC)),
            PCTN = mean(PctN),
            PCTN_se = sd(PctN)/sqrt(length(PctN)),
            CN = mean(CNx),
            CN_se = sd(CNx)/sqrt(length(CNx)),
            d15N = mean(d15Nx),
            d15N_se = sd(d15Nx)/sqrt(length(d15Nx)),
            SLAgm2 = mean(SLA_gm2),
            SLAse = sd(SLA_gm2, na.rm=T)/sqrt(length(SLA_gm2)))
write.csv(leafchem_summary_TREE, "leafchem_summary_tree.csv")

  


#####Analysis: Foliar and canopy data#####

#for "without mature forests" tests: dataset without mature forest
names(CANOPY)
levels(CANOPY$CLASS)
CANOPY_nomat <- CANOPY %>%
  filter(CLASS != "C4") %>%
  droplevels()

#dLAI
#dLAI One-factor ANOVA: Effect of site age class on dLAI.
LAIwetdrydiff_model1 <- lm(CANOPY$LAIWETEDRYMIN~CANOPY$CLASS)
LAIwetdrydiff_model1b <- aov(CANOPY$LAIWETEDRYMIN~CANOPY$CLASS)
anova(LAIwetdrydiff_model1)
summary(LAIwetdrydiff_model1)
TukeyHSD(LAIwetdrydiff_model1b)

#dLAI log age model
names(CANOPY)
LAIwetdrydiff_modelcont <- lm(CANOPY$LAIWETEDRYMIN~CANOPY$natlogage)
summary(LAIwetdrydiff_modelcont)

#dLAI log age model without mature forests
LAIwetdrydiff_modelcont_nomat <- lm(CANOPY_nomat$LAIWETEDRYMIN~CANOPY_nomat$natlogage)
summary(LAIwetdrydiff_modelcont_nomat)

#Canopy N demand
#Canopy N demand ANOVA: Effect of site age class on Canopy N demand.
CanopyNdemand_model1 <- lm(CANOPY$CanopyN_kgha~CANOPY$CLASS)
CanopyNdemand_modelaov <- aov(CANOPY$CanopyN_kgha~CANOPY$CLASS)
anova(CanopyNdemand_model1)
summary(CanopyNdemand_model1)
TukeyHSD(CanopyNdemand_modelaov)

#Canopy N demand log age model
CanopyNdemand_logagemodel <- lm(CANOPY$CanopyN_kgha~CANOPY$natlogage)
summary(CanopyNdemand_logagemodel)

#Canopy N demand log age model without mature forests
CanopyNdemand_logagemodel_nomat <- lm(CANOPY_nomat$CanopyN_kgha~CANOPY_nomat$natlogage)
summary(CanopyNdemand_logagemodel_nomat)


#Foliar data: trees as observational units

#  foliar dataset without mature sites.
leafchem_treeobs_nomat <- leafchem_treeobs %>%
  filter(CLASS != "VIL") %>%
  droplevels()

#add site ages to leafobs dataset#
write.csv(leafchem_treeobs, "leafchem_treeobs.csv")
leafchem_treeobs_agesin <- read.csv("leafchem_treeobs_agesin.csv")

#foliar dataset, ages in without mature sites
leafchem_treeobs_agesin_nomat <- leafchem_treeobs_agesin %>%
  filter(CLASS != "VIL") %>%
  droplevels()

#procedure: code block for foliar %N and d15N repeated below each section's analysis with the main data table called changed from the one WITH mature forests (leafchem_treeobs) to the one WITHOUT mature forests (leafchem_treeobs_nomat).  Results of various tests will be overwritten.

# Foliar %N
names(leafchem_treeobs)
levels(leafchem_treeobs$CLASS)
par(mfrow=c(2,3))
hist(leafchem_treeobs$PctN)
hist((leafchem_treeobs$PctN)^1.38)

shapiro.test(leafchem_treeobs$PctN)
shapiro.test((leafchem_treeobs$PctN)^1.38)
#data MARGINALLY different from normal; will do non-parametric stats.

# Foliar %N - species effects
boxplot(PctN~TREE, leafchem_treeobs)
kruskal.test(leafchem_treeobs$PctN~leafchem_treeobs$TREE)
# %N - species differences by group across all classes
MTvsMZcomp <- leafchem_treeobs %>% filter(TREE == "MT" | TREE == "MZ") %>% droplevels()
wilcox.test(MTvsMZcomp$PctN~MTvsMZcomp$TREE); #different

MTvsMNcomp <- leafchem_treeobs %>% filter(TREE == "MT" | TREE == "MN") %>% droplevels()
wilcox.test(MTvsMNcomp$PctN~MTvsMNcomp$TREE); #different

MNvsMZcomp <- leafchem_treeobs %>% filter(TREE == "MN" | TREE == "MZ") %>% droplevels()
wilcox.test(MNvsMZcomp$PctN~MNvsMZcomp$TREE)


# Foliar %N - age effects
boxplot(PctN~CLASS, leafchem_treeobs)
kruskal.test(leafchem_treeobs$PctN~leafchem_treeobs$CLASS)
C01vsVIL <- leafchem_treeobs %>% filter(CLASS == "C01" | CLASS == "VIL") %>%droplevels()
wilcox.test(C01vsVIL$PctN~C01vsVIL$CLASS)

C01vsC02 <- leafchem_treeobs %>% filter(CLASS == "C01" | CLASS == "C02") %>%droplevels()
wilcox.test(C01vsC02$PctN~C01vsC02$CLASS)

C01vsC03 <- leafchem_treeobs %>% filter(CLASS == "C01" | CLASS == "C03") %>%droplevels()
wilcox.test(C01vsC03$PctN~C01vsC03$CLASS)

C02vsC03 <- leafchem_treeobs %>% filter(CLASS == "C02" | CLASS == "C03") %>% droplevels() 
wilcox.test(C02vsC03$PctN~C02vsC03$CLASS)

C02vsVIL <- leafchem_treeobs %>% filter(CLASS == "C02" | CLASS == "VIL") %>% droplevels()
wilcox.test(C02vsVIL$PctN~C02vsVIL$CLASS)

C03vsVIL <- leafchem_treeobs %>% filter(CLASS == "C03" | CLASS == "VIL") %>% droplevels()
wilcox.test(C03vsVIL$PctN~C03vsVIL$CLASS)

#Foliar %N - age effects, excluding mature forests
kruskal.test(leafchem_treeobs_nomat$PctN~leafchem_treeobs_nomat$CLASS)
#Foliar $N - tree effects, excluding mature forests
kruskal.test(leafchem_treeobs_nomat$PctN~leafchem_treeobs_nomat$TREE)

##foliar %N comparison - mature vs. intermediate age sites grouped.
leafchem_treeobs_intvsmatcomp <- read.csv("leafchem_treeobs_intvsmatforestcomp.csv")
names(leafchem_treeobs_intvsmatcomp)
levels(leafchem_treeobs_intvsmatcomp$COMP2)

leafchem_treeobs_intvsmatcomp_C01out <- leafchem_treeobs_intvsmatcomp %>%
  filter(COMP2 != "C01") %>%
  droplevels()

kruskal.test(leafchem_treeobs_intvsmatcomp_C01out$PctN~leafchem_treeobs_intvsmatcomp_C01out$COMP2)
wilcox.test(leafchem_treeobs_intvsmatcomp_C01out$PctN~leafchem_treeobs_intvsmatcomp_C01out$COMP2)

# Foliar %N age effects by tree species
boxplot(PctN~CLASS.TREE, leafchem_treeobs)
kruskal.test(leafchem_treeobs$PctN~leafchem_treeobs$CLASS.TREE)

#mninga PctN variation by class
C01MNvsC02MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MN" | CLASS.TREE == "C02.MN") %>% droplevels()
wilcox.test(C01MNvsC02MN$PctN~C01MNvsC02MN$CLASS)

C01MNvsC03MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MN" | CLASS.TREE == "C03.MN") %>% droplevels()
wilcox.test(C01MNvsC03MN$PctN~C01MNvsC03MN$CLASS)

C01MNvsC04MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MN" | CLASS.TREE == "VIL.MN") %>% droplevels()
wilcox.test(C01MNvsC04MN$PctN~C01MNvsC04MN$CLASS)

C02MNvsC03MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MN" | CLASS.TREE == "C03.MN") %>% droplevels()
wilcox.test(C02MNvsC03MN$PctN~C02MNvsC03MN$CLASS)

C02MNvsC04MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MN" | CLASS.TREE == "VIL.MN") %>% droplevels()
wilcox.test(C02MNvsC04MN$PctN~C02MNvsC04MN$CLASS)

C03MNvsC04MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C03.MN" | CLASS.TREE == "VIL.MN") %>% droplevels()
wilcox.test(C03MNvsC04MN$PctN~C03MNvsC04MN$CLASS)

#mzima Pct15N variation by class
C01MZvsC02MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MZ" | CLASS.TREE == "C02.MZ") %>% droplevels()
wilcox.test(C01MZvsC02MZ$PctN~C01MZvsC02MZ$CLASS)

C01MZvsC03MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MZ" | CLASS.TREE == "C03.MZ") %>% droplevels()
wilcox.test(C01MZvsC03MZ$PctN~C01MZvsC03MZ$CLASS)

C01MZvsC04MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MZ" | CLASS.TREE == "VIL.MZ") %>% droplevels()
wilcox.test(C01MZvsC04MZ$PctN~C01MZvsC04MZ$CLASS)

C02MZvsC03MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MZ" | CLASS.TREE == "C03.MZ") %>% droplevels()
wilcox.test(C02MZvsC03MZ$PctN~C02MZvsC03MZ$CLASS)

C02MZvsC04MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MZ" | CLASS.TREE == "VIL.MZ") %>% droplevels()
wilcox.test(C02MZvsC04MZ$PctN~C02MZvsC04MZ$CLASS)

C03MZvsC04MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C03.MZ" | CLASS.TREE == "VIL.MZ") %>% droplevels()
wilcox.test(C03MZvsC04MZ$PctN~C03MZvsC04MZ$CLASS)

#mtundu Pct15N variation by class
C01MTvsC02MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MT" | CLASS.TREE == "C02.MT") %>% droplevels()
wilcox.test(C01MTvsC02MT$PctN~C01MTvsC02MT$CLASS)

C01MTvsC03MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MT" | CLASS.TREE == "C03.MT") %>% droplevels()
wilcox.test(C01MTvsC03MT$PctN~C01MTvsC03MT$CLASS)

C01MTvsC04MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MT" | CLASS.TREE == "VIL.MT") %>% droplevels()
wilcox.test(C01MTvsC04MT$PctN~C01MTvsC04MT$CLASS)

C02MTvsC03MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MT" | CLASS.TREE == "C03.MT") %>% droplevels()
wilcox.test(C02MTvsC03MT$PctN~C02MTvsC03MT$CLASS)

C02MTvsC04MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MT" | CLASS.TREE == "VIL.MT") %>% droplevels()
wilcox.test(C02MTvsC04MT$PctN~C02MTvsC04MT$CLASS)

C03MTvsC04MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C03.MT" | CLASS.TREE == "VIL.MT") %>% droplevels()
wilcox.test(C03MTvsC04MT$PctN~C03MTvsC04MT$CLASS)



# Foliar d15N - trees as observations
hist(leafchem_treeobs$d15Nx)
hist(log(leafchem_treeobs$d15Nx+3))
shapiro.test(log(leafchem_treeobs$d15Nx+3))
#not easily transformable because of species differences.

boxplot(d15Nx~TREE, leafchem_treeobs)
kruskal.test(leafchem_treeobs$d15Nx~leafchem_treeobs$TREE)

# Foliar d15N - species differences by group across all classes
MTvsMZcomp <- leafchem_treeobs %>% filter(TREE == "MT" | TREE == "MZ") %>% droplevels()
wilcox.test(MTvsMZcomp$d15Nx~MTvsMZcomp$TREE); #different

MTvsMNcomp <- leafchem_treeobs %>% filter(TREE == "MT" | TREE == "MN") %>% droplevels()
wilcox.test(MTvsMNcomp$d15Nx~MTvsMNcomp$TREE); #different

MNvsMZcomp <- leafchem_treeobs %>% filter(TREE == "MN" | TREE == "MZ") %>% droplevels()
wilcox.test(MNvsMZcomp$d15Nx~MNvsMZcomp$TREE)


# Foliar d15N - age effects
boxplot(d15Nx~CLASS, leafchem_treeobs)
kruskal.test(leafchem_treeobs$d15Nx~leafchem_treeobs$CLASS)
C01vsVIL <- leafchem_treeobs %>% filter(CLASS == "C01" | CLASS == "VIL") %>%droplevels()
wilcox.test(C01vsVIL$d15Nx~C01vsVIL$CLASS)

C01vsC02 <- leafchem_treeobs %>% filter(CLASS == "C01" | CLASS == "C02") %>%droplevels()
wilcox.test(C01vsC02$d15Nx~C01vsC02$CLASS)

C01vsC03 <- leafchem_treeobs %>% filter(CLASS == "C01" | CLASS == "C03") %>%droplevels()
wilcox.test(C01vsC03$d15Nx~C01vsC03$CLASS)

C02vsC03 <- leafchem_treeobs %>% filter(CLASS == "C02" | CLASS == "C03") %>% droplevels() 
wilcox.test(C02vsC03$d15Nx~C02vsC03$CLASS)

C02vsVIL <- leafchem_treeobs %>% filter(CLASS == "C02" | CLASS == "VIL") %>% droplevels()
wilcox.test(C02vsVIL$d15Nx~C02vsVIL$CLASS)

C03vsVIL <- leafchem_treeobs %>% filter(CLASS == "C03" | CLASS == "VIL") %>% droplevels()
wilcox.test(C03vsVIL$d15Nx~C03vsVIL$CLASS)

#Foliar d15N - age effects, excluding mature forests
kruskal.test(leafchem_treeobs_nomat$d15Nx~leafchem_treeobs_nomat$CLASS)
#Foliar d15N - species effects, excluding mature forests
kruskal.test(leafchem_treeobs_nomat$d15Nx~leafchem_treeobs_nomat$TREE)


# Foliar d15N age effects by tree species
boxplot(d15Nx~CLASS.TREE, leafchem_treeobs)
kruskal.test(leafchem_treeobs$d15Nx~leafchem_treeobs$CLASS.TREE)

#mninga d15Nx variation by class
C01MNvsC02MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MN" | CLASS.TREE == "C02.MN") %>% droplevels()
wilcox.test(C01MNvsC02MN$d15Nx~C01MNvsC02MN$CLASS)

C01MNvsC03MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MN" | CLASS.TREE == "C03.MN") %>% droplevels()
wilcox.test(C01MNvsC03MN$d15Nx~C01MNvsC03MN$CLASS)

C01MNvsC04MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MN" | CLASS.TREE == "VIL.MN") %>% droplevels()
wilcox.test(C01MNvsC04MN$d15Nx~C01MNvsC04MN$CLASS)

C02MNvsC03MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MN" | CLASS.TREE == "C03.MN") %>% droplevels()
wilcox.test(C02MNvsC03MN$d15Nx~C02MNvsC03MN$CLASS)

C02MNvsC04MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MN" | CLASS.TREE == "VIL.MN") %>% droplevels()
wilcox.test(C02MNvsC04MN$d15Nx~C02MNvsC04MN$CLASS)

C03MNvsC04MN <- leafchem_treeobs %>% filter(CLASS.TREE == "C03.MN" | CLASS.TREE == "VIL.MN") %>% droplevels()
wilcox.test(C03MNvsC04MN$d15Nx~C03MNvsC04MN$CLASS)

#mzima d15N variation by class
C01MZvsC02MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MZ" | CLASS.TREE == "C02.MZ") %>% droplevels()
wilcox.test(C01MZvsC02MZ$d15Nx~C01MZvsC02MZ$CLASS)

C01MZvsC03MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MZ" | CLASS.TREE == "C03.MZ") %>% droplevels()
wilcox.test(C01MZvsC03MZ$d15Nx~C01MZvsC03MZ$CLASS)

C01MZvsC04MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MZ" | CLASS.TREE == "VIL.MZ") %>% droplevels()
wilcox.test(C01MZvsC04MZ$d15Nx~C01MZvsC04MZ$CLASS)

C02MZvsC03MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MZ" | CLASS.TREE == "C03.MZ") %>% droplevels()
wilcox.test(C02MZvsC03MZ$d15Nx~C02MZvsC03MZ$CLASS)

C02MZvsC04MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MZ" | CLASS.TREE == "VIL.MZ") %>% droplevels()
wilcox.test(C02MZvsC04MZ$d15Nx~C02MZvsC04MZ$CLASS)

C03MZvsC04MZ <- leafchem_treeobs %>% filter(CLASS.TREE == "C03.MZ" | CLASS.TREE == "VIL.MZ") %>% droplevels()
wilcox.test(C03MZvsC04MZ$d15Nx~C03MZvsC04MZ$CLASS)

#mtundu d15N variation by class
C01MTvsC02MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MT" | CLASS.TREE == "C02.MT") %>% droplevels()
wilcox.test(C01MTvsC02MT$d15Nx~C01MTvsC02MT$CLASS)

C01MTvsC03MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MT" | CLASS.TREE == "C03.MT") %>% droplevels()
wilcox.test(C01MTvsC03MT$d15Nx~C01MTvsC03MT$CLASS)

C01MTvsC04MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C01.MT" | CLASS.TREE == "VIL.MT") %>% droplevels()
wilcox.test(C01MTvsC04MT$d15Nx~C01MTvsC04MT$CLASS)

C02MTvsC03MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MT" | CLASS.TREE == "C03.MT") %>% droplevels()
wilcox.test(C02MTvsC03MT$d15Nx~C02MTvsC03MT$CLASS)

C02MTvsC04MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C02.MT" | CLASS.TREE == "VIL.MT") %>% droplevels()
wilcox.test(C02MTvsC04MT$d15Nx~C02MTvsC04MT$CLASS)

C03MTvsC04MT <- leafchem_treeobs %>% filter(CLASS.TREE == "C03.MT" | CLASS.TREE == "VIL.MT") %>% droplevels()
wilcox.test(C03MTvsC04MT$d15Nx~C03MTvsC04MT$CLASS)


# Specific leaf area 
# SLA - trees as observations
hist(leafchem_treeobs$SLA_gm2)
hist(log(leafchem_treeobs$SLA_gm2))
shapiro.test(log(leafchem_treeobs$SLA_gm2))
#YES!!  Do parametric stats for SLA.  Normal at the tree level.

#Age class effects on SLA
boxplot(SLA_gm2~CLASS, leafchem_treeobs)
SLA_model1 <- lm(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$CLASS)
SLA_model1b <- aov(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$CLASS)
anova(SLA_model1)
summary(SLA_model1)

#Species effects on SLA
boxplot(SLA_gm2~TREE, leafchem_treeobs)
SLA_model2 <- lm(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$TREE)
SLA_model2b <- aov(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$TREE)
anova(SLA_model2)
summary(SLA_model2)
TukeyHSD(SLA_model2b)

#Two-factor ANOVA: age + species effects
SLA_model3 <- lm(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$TREE+leafchem_treeobs$CLASS)
SLA_model3b <- aov(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$TREE+leafchem_treeobs$CLASS)
anova(SLA_model3)
summary(SLA_model3)
TukeyHSD(SLA_model3b, bonferroni = T)


#Two-factor ANOVA with interactions: age*species
SLA_model4 <- lm(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$TREE*leafchem_treeobs$CLASS)
SLA_model4b <- aov(log(leafchem_treeobs$SLA_gm2)~leafchem_treeobs$TREE*leafchem_treeobs$CLASS)
anova(SLA_model4)
summary(SLA_model4)
TukeyHSD(SLA_model4b, bonferroni = T)

AIC(SLA_model3,SLA_model4)
#difference in AIC is pretty modest <10.  Go with additive effects.


#log age analysis
names(leafchem_treeobs)
#add site ages to leafobs dataset#
write.csv(leafchem_treeobs, "leafchem_treeobs.csv")
leafchem_treeobs_agesin <- read.csv("leafchem_treeobs_agesin.csv")

#age alone
SLA_logagemodela <- lm(log(leafchem_treeobs_agesin$SLA_gm2)~log(leafchem_treeobs_agesin$AGE))
summary(SLA_logagemodela)

#age + species
SLA_logagemodel1 <- lm(log(leafchem_treeobs_agesin$SLA_gm2)~log(leafchem_treeobs_agesin$AGE)+leafchem_treeobs_agesin$TREE)
summary(SLA_logagemodel1)
anova(SLA_logagemodel1)

SLA_logagemodel2 <- lm(log(leafchem_treeobs_agesin$SLA_gm2)~log(leafchem_treeobs_agesin$AGE)*leafchem_treeobs_agesin$TREE)
summary(SLA_logagemodel2)
anova(SLA_logagemodel2)

AIC(SLA_logagemodel1, SLA_logagemodel2)
#equivalent, more parsimonious one is without interactions.



#log age analysis, excluding mature forests
#age alone
SLA_logagemodel1_nomat <- lm(log(leafchem_treeobs_agesin_nomat$SLA_gm2)~log(leafchem_treeobs_agesin_nomat$AGE))
anova(SLA_logagemodel1_nomat)
summary(SLA_logagemodel1_nomat)

#age + species
SLA_logagemodel2_nomat <- lm(log(leafchem_treeobs_agesin_nomat$SLA_gm2)~log(leafchem_treeobs_agesin_nomat$AGE)+leafchem_treeobs_agesin_nomat$TREE)
anova(SLA_logagemodel2_nomat)
summary(SLA_logagemodel2_nomat)


#####Soil data and analysis:  end-stage code (update file with wrangle later)#####

#ENH4NO3_mgkgsoil_t0

#*use AN(C)OVA and 1/(1+) transform for mgNNH4kgsoil for the statistics.
#basic linear model - lm
?lm
?lme
?AIC
bartlett.test((1/(NMIN_adry_t0_D1andD2$EmgNNH4NO3kgsoil+1)~NMIN_adry_t0_D1andD2$CLASS))
#Sum passes Bartlett tests for class, depth

#dry season soil Nmin
hist(NMIN_adry_t0_D1andD2$EmgNNH4NO3kgsoil)
hist(1/NMIN_adry_t0_D1andD2$EmgNNH4NO3kgsoil)
shapiro.test(1/NMIN_adry_t0_D1andD2$EmgNNH4NO3kgsoil)
#not different from normal.

#CLASS by itself is not significant model
ENH4NO3_adry_t0_D1andD2_test <- lm((1/(EmgNNH4NO3kgsoil+1))~CLASS, NMIN_adry_t0_D1andD2)
summary(ENH4NO3_adry_t0_D1andD2_test)
anova(ENH4NO3_adry_t0_D1andD2_test)

#DEPTH by itself is a significant model
ENH4NO3_adry_t0_D1andD2_test1b <- lm((1/(EmgNNH4NO3kgsoil+1))~DEPTH, NMIN_adry_t0_D1andD2)
summary(ENH4NO3_adry_t0_D1andD2_test1b)
anova(ENH4NO3_adry_t0_D1andD2_test1b)

#SITE is not a significant model
ENH4NO3_adry_t0_D1andD2_test1c <- lm((1/(EmgNNH4NO3kgsoil+1))~SITE, NMIN_adry_t0_D1andD2)
summary(ENH4NO3_adry_t0_D1andD2_test1c)
anova(ENH4NO3_adry_t0_D1andD2_test1c)

#CLASS+DEPTH model significant; 
ENH4NO3_adry_t0_D1andD2_test2 <- lm((1/(EmgNNH4kgsoil+1))~CLASS+DEPTH, NMIN_adry_t0_D1andD2)
anova(ENH4NO3_adry_t0_D1andD2_test2)
?TukeyHSD
TukeyHSD(aov(lm((1/(mgNNH4kgsoil+1))~CLASS+DEPTH, NMIN_adry_t0_D1andD2)), bonferroni = T)
summary(ENH4NO3_adry_t0_D1andD2_test2)

ENH4NO3_adry_t0_D1andD2_test2 <- lm((1/(mgNNH4kgsoil+1))~CLASS*DEPTH, NMIN_adry_t0_D1andD2)
anova(ENH4NO3_adry_t0_D1andD2_test2)
TukeyHSD(aov(lm((1/(mgNNH4kgsoil+1))~CLASS*DEPTH, NMIN_adry_t0_D1andD2)))
summary(ENH4NO3_adry_t0_D1andD2_test2)



#CLASS+DEPTH+SOILTEMP+MOIST
ENH4NO3_adry_t0_D1andD2_test5 <- lm((1/(EmgNNH4NO3kgsoil+1))~CLASS+DEPTH+MOIST+SOILTEMP, NMIN_adry_t0_D1andD2)
ENH4NO3_adry_t0_D1andD2_test5b <- aov(lm((1/(EmgNNH4NO3kgsoil+1))~CLASS+DEPTH+MOIST+SOILTEMP, NMIN_adry_t0_D1andD2))
anova(ENH4NO3_adry_t0_D1andD2_test5)
summary(ENH4NO3_adry_t0_D1andD2_test5)
TukeyHSD(ENH4NO3_adry_t0_D1andD2_test5b, bonferonni=T)

#log age model - age alone
ENH4NO3_adry_t0_D1andD2_logagetest1 <- lm((1/EmgNNH4NO3kgsoil+1)~logage,NMIN_adry_t0_D1andD2)
anova(ENH4NO3_adry_t0_D1andD2_logagetest1)
summary(ENH4NO3_adry_t0_D1andD2_logagetest1)

#log age model - age, depth, temp, soil moiustre
ENH4NO3_adry_t0_D1andD2_logagetest2 <- lm((1/EmgNNH4NO3kgsoil+1)~logage+DEPTH,NMIN_adry_t0_D1andD2)
anova(ENH4NO3_adry_t0_D1andD2_logagetest2)
summary(ENH4NO3_adry_t0_D1andD2_logagetest2)


####repeat analysis, but exclude mature forests from the data called, NMIN_adry_t0_D1andD2.
names(NMIN_adry_t0_D1andD2)
levels(NMIN_adry_t0_D1andD2$CLASS)

NMIN_adry_t0_D1andD2_nomat <- NMIN_adry_t0_D1andD2 %>%
  filter(CLASS != "C4") %>%
  droplevels()

#log age model - age alone
ENH4NO3_adry_t0_D1andD2_nomat_logagetest1 <- lm((1/EmgNNH4NO3kgsoil+1)~logage,NMIN_adry_t0_D1andD2_nomat)
anova(ENH4NO3_adry_t0_D1andD2_nomat_logagetest1)
summary(ENH4NO3_adry_t0_D1andD2_nomat_logagetest1)

#log age model - age + depth
ENH4NO3_adry_t0_D1andD2_nomat_logagetest2 <- lm((1/EmgNNH4NO3kgsoil+1)~logage+DEPTH,NMIN_adry_t0_D1andD2_nomat)
anova(ENH4NO3_adry_t0_D1andD2_nomat_logagetest2)
summary(ENH4NO3_adry_t0_D1andD2_logagetest2)


#Soil total nitrogen analysis
hist(NMIN_adry_t0_D1andD2$N_Mghax)
shapiro.test(NMIN_adry_t0_D1andD2$N_Mghax)
#Data normal as they are.  No transform.

#model with age effects alone, by class.
TN_class_model1 <- lm(N_Mghax~CLASS, NMIN_adry_t0_D1andD2)
anova(TN_class_model1)
summary(TN_class_model1)

#other models
TN_class_model2 <- lm(N_Mghax~CLASS+DEPTH, NMIN_adry_t0_D1andD2)
summary(TN_class_model2)
anova(TN_class_model2)

TN_class_model3 <- lm(N_Mghax~CLASS+DEPTH+MOIST+SOILTEMP, NMIN_adry_t0_D1andD2)
summary(TN_class_model3)
anova(TN_class_model3)

TN_class_model4 <- lm(N_Mghax~DEPTH, NMIN_adry_t0_D1andD2)
summary(TN_class_model4)

##Log age model - none others needed, no clear pattern.
TN_class_model5 <- lm(N_Mghax~logage,NMIN_adry_t0_D1andD2)
summary(TN_class_model5)

##log age model + depth as co-variate
TN_class_model6 <- lm(N_Mghax~logage+DEPTH,NMIN_adry_t0_D1andD2)
summary(TN_class_model6)

#Log age model excluding mature forests 
TN_age_model7 <- lm(N_Mghax~logage, NMIN_adry_t0_D1andD2_nomat)
summary(TN_age_model7)

#Log age model exlcuding mature forests, age+ depth
TN_age_model8 <- lm(N_Mghax~logage+DEPTH, NMIN_adry_t0_D1andD2_nomat)
summary(TN_age_model8)


#Soil NO3:NH4, excluding mature forests
kruskal.test(NMIN_adry_t0_D1andD2_nomat$NO3toNH4mol~NMIN_adry_t0_D1andD2_nomat$CLASS)
kruskal.test(NMIN_adry_t0_D1andD2_nomat$NO3toNH4mol~NMIN_adry_t0_D1andD2_nomat$DEPTH)


#soil d15N.
names(NMIN_adry_t0_D1andD2)
#soil d15N analysis by class, including mature forests.
hist(NMIN_adry_t0_D1andD2$d15N)
hist(log(NMIN_adry_t0_D1andD2$d15N))
shapiro.test(log(NMIN_adry_t0_D1andD2$d15N))
#great.  d15N log-transformed is good for AN(C)OVA
hist(NMIN_adry_t0_D1andD2$MOIST)
hist(1/(NMIN_adry_t0_D1andD2$MOIST+1))
shapiro.test(1/(NMIN_adry_t0_D1andD2$MOIST+1))

#Significant effects of site age and moisture on d15N.

d15N_age <- lm(log(d15Nx)~CLASS, NMIN_adry_t0_D1andD2)
anova(d15N_age)
summary(d15N_age)

##d15N - ANOVA
d15N_age_model1 <- lm((log(d15Nx))~CLASS+DEPTH+SOILTEMP+MOIST, NMIN_adry_t0_D1andD2)
d15N_age_model1b <- aov(lm((log(d15Nx))~CLASS+DEPTH+SOILTEMP+MOIST, NMIN_adry_t0_D1andD2))
anova(d15N_age_model1)
summary(d15N_age_model1)
TukeyHSD(d15N_age_model1b)

plot(d15Nx~MOIST, NMIN_adry_t0_D1andD2)
plot(log(d15Nx)~SOILTEMP, NMIN_adry_t0_D1andD2)
plot(log(d15Nx)~SOILTEMP, NMIN_awet_t0_D1andD2)

plot(d15Nx~logage, NMIN_adry_t0_D1andD2)

##d15N - log age (mature forests included)
d15N_logage_model1 <- lm(log(d15Nx)~logage, NMIN_adry_t0_D1andD2)
summary(d15N_logage_model1)

d15N_logage_model2dry <- lm(log(d15Nx)~logage+DEPTH+SOILTEMP, NMIN_adry_t0_D1andD2)
summary(d15N_logage_model2dry)

##d15N - log age (mature forests excluded

d15N_logage_nomat_model1 <- lm(log(d15Nx)~logage, NMIN_adry_t0_D1andD2_nomat)
summary(d15N_logage_nomat_model1)

d15N_logage_nomat_model2 <- lm(log(d15Nx)~logage+DEPTH+SOILTEMP, NMIN_adry_t0_D1andD2_nomat)
summary(d15N_logage_nomat_model2)


###soil organic C analysis: just do along with the soil N, forms part of your argument that SOC is not building up with forest regrowth.
hist(NMIN_adry_t0_D1andD2$C_Mghax)
shapiro.test(NMIN_adry_t0_D1andD2$C_Mghax)
#Data normal as they are.  No transform.

#model with age effects alone, by class.
SOC_class_model1 <- lm(C_Mghax~CLASS, NMIN_adry_t0_D1andD2)
anova(SOC_class_model1)
summary(SOC_class_model1)

#other models
SOC_class_model2 <- lm(C_Mghax~CLASS+DEPTH, NMIN_adry_t0_D1andD2)
summary(SOC_class_model2)
anova(SOC_class_model2)

SOC_class_model3 <- lm(C_Mghax~CLASS+DEPTH+MOIST+SOILTEMP, NMIN_adry_t0_D1andD2)
summary(SOC_class_model3)
anova(SOC_class_model3)

##Log age model - none others needed, no clear pattern.
SOC_class_model5 <- lm(C_Mghax~logage,NMIN_adry_t0_D1andD2)
summary(SOC_class_model5)

##log age model + depth as co-variate
SOC_class_model6 <- lm(C_Mghax~logage+DEPTH,NMIN_adry_t0_D1andD2)
summary(SOC_class_model6)

#Log age model excluding mature forests 
SOC_age_model7 <- lm(C_Mghax~logage, NMIN_adry_t0_D1andD2_nomat)
summary(SOC_age_model7)

#Log age model exlcuding mature forests, age+ depth
SOC_age_model8 <- lm(C_Mghax~logage+DEPTH, NMIN_adry_t0_D1andD2_nomat)
summary(SOC_age_model8)



#####Soil Nmin analyses#####
names(NMIN_afull_t2mt0_DSUM_FM)
levels(NMIN_afull_t2mt0_DSUM_FM$SEASON)
levels(NMIN_afull_t2mt0_DSUM_FM$CLASS)
hist(NMIN_adry_t2mt0_DSUM_FM$EkgNNH4NO3ha)
shapiro.test(NMIN_adry_t2mt0_DSUM_FM$EkgNNH4NO3ha)
#data not different from normal untransformed.
hist(log(NMIN_adry_t2mt0_DSUM_FM$EkgNNH4NO3ha+1))
shapiro.test(log(NMIN_adry_t2mt0_DSUM_FM$EkgNNH4NO3ha+1))

#filter NMIN_afull to exclude mature forests
levels(NMIN_afull_t2mt0_DSUM_FM$CLASS)

NMIN_afull_t2mt0_DSUM_FM_nomat <- NMIN_afull_t2mt0_DSUM_FM %>%
  filter(CLASS != "C4") %>%
  droplevels()
NMIN_afull_t2mt0_DSUM_FM_nomat
write.csv(NMIN_afull_t2mt0_DSUM_FM_nomat, "NMIN_afull_t2mt0_DSUM_FM_nomat.csv")

#Nmin, both seasons: SEASON+AGE in Two-factor ANOVA
NMIN_seasonalvsclassvar <- lm(EkgNNH4NO3ha~SEASON+CLASS, NMIN_afull_t2mt0_DSUM_FM)
NMIN_seasonalvsclassvarb <- aov(EkgNNH4NO3ha~SEASON+CLASS, NMIN_afull_t2mt0_DSUM_FM)
anova(NMIN_seasonalvsclassvar)
summary(NMIN_seasonalvsclassvar)
TukeyHSD(NMIN_seasonalvsclassvarb, bonferroni = T)

#NMIN, both seasons: SEASON_logAGE for log age model
NMIN_seasonvsclass_logage <- lm(EkgNNH4NO3ha~logage+SEASON, NMIN_afull_t2mt0_DSUM_FM)
summary(NMIN_seasonvsclass_logage)

#Nmin, both seasons, log age model excluding mature forests
NMIN_seasonvsclass_logage_nomat <- lm(EkgNNH4NO3ha~logage+SEASON, NMIN_afull_t2mt0_DSUM_FM_nomat)
summary(NMIN_seasonvsclass_logage_nomat)


##Nmin, wet season
names(NMIN_awet_t2mt0_DSUM_FM)
hist(NMIN_awet_t2mt0_DSUM_FM$EkgNNH4NO3ha)
shapiro.test(NMIN_awet_t2mt0_DSUM_FM$EkgNNH4NO3ha)
#not different than normal, analyses done with raw data.

#filter NMIN_awet to exclude mature forests
levels(NMIN_afull_t2mt0_DSUM_FM$CLASS)

NMIN_awet_t2mt0_DSUM_FM_nomat <- NMIN_awet_t2mt0_DSUM_FM %>%
  filter(CLASS != "C4") %>%
  droplevels()


#Nmin, wet season: AGE in One-factor ANOVA
NMIN_wetvsclassvar <- lm(EkgNNH4NO3ha~CLASS, NMIN_awet_t2mt0_DSUM_FM)
NMIN_wetvsclassvarb <- aov(EkgNNH4NO3ha~CLASS, NMIN_awet_t2mt0_DSUM_FM)
anova(NMIN_wetvsclassvar)
summary(NMIN_wetvsclassvar)
TukeyHSD(NMIN_wetvsclassvarb, bonferroni = T)

#NMIN, both seasons: SEASON_logAGE for log age model
NMIN_wetvsclass_logage <- lm(EkgNNH4NO3ha~logage, NMIN_awet_t2mt0_DSUM_FM)
summary(NMIN_wetvsclass_logage)

#Nmin, both seasons, log age model excluding mature forests
NMIN_seasonvsclass_logage_nomat <- lm(EkgNNH4NO3ha~logage+SEASON, NMIN_afull_t2mt0_DSUM_FM_nomat)
summary(NMIN_seasonvsclass_logage_nomat)

#####New N budget analyses for 4/25/2016####

NBUD <- read.csv("Nbudgetin_20160425.csv")
names(NBUD)

NBUD_sumbyclass <- NBUD %>%
  group_by(CLASS) %>%
  summarize(N = length(N_kgha),
            SoilToTNkgha = mean(N_kgha),
            SoilToTNkgha_se = sd(N_kgha)/sqrt(length(N_kgha)),
            SoilNmin_kgha = mean(AnnNmin50pctwdscale_kgha),
            SoilNmin_kgha_se = sd(AnnNmin50pctwdscale_kgha)/sqrt(length(AnnNmin50pctwdscale_kgha)),
            WdBiomassNkgha = mean(WdBiomassN),
            WdBiomassNkgha_se = sd(WdBiomassN)/sqrt(length(WdBiomassN)),
            CanopyNkgha = mean(CanopyN_kgha),
            CanopyNkgha_se = sd(CanopyN_kgha)/sqrt(length(CanopyN_kgha)),
            TotalTreeNkgha = mean(TotalTreeN_kgha),
            TotalTreeNkgha_se = sd(TotalTreeN_kgha)/sqrt(length(TotalTreeN_kgha)))
            
write.csv(NBUD_sumbyclass, "NBUD_sumbyclass.csv")

NBUD_nomat <- NBUD %>%
  filter(CLASS != "C4") %>%
  droplevels()


hist(NBUD$N_kgha)
shapiro.test(NBUD$N_kgha)
#soil total N = normal.
hist(NBUD$AnnNmin50pctwdscale_kgha)
shapiro.test(NBUD$AnnNmin50pctwdscale_kgha)
#soil Nmin-annual scale = normal.
hist(log(NBUD$WdBiomassN))
shapiro.test(log(NBUD$WdBiomassN))
#WdBiomassN = log-transform for normality.
hist(NBUD$CanopyN_kgha)
shapiro.test(NBUD$CanopyN_kgha)
#CanopyN = normal.
hist(NBUD$TotalTreeN_kgha)
hist(log(NBUD$TotalTreeN_kgha))
shapiro.test(log(NBUD$TotalTreeN_kgha))
#TotalTreeN = log-normal.

##ANOVA models for site age effects (by class) on new aggregated N

##super-efficient table-making model
#Nstembyageclass.lm = lm(NSTEMHA~CLASS, TSV1)
#Nstemageclassmeanscomp <- aov(NSTEMHA~CLASS, TSV1)
#anova(Nstembyageclass.lm)
#summary(Nstembyageclass.lm)
#significant effects of class on stem density (p = 0.001)
#TukeyHSD(aov(lm(NSTEMHA~CLASS, TSV1)), ordered=T))
#HSD.test(Nstemageclassmeanscomp, "CLASS", group=TRUE)
#print(HSD.test(Nstemageclassmeanscomp, "CLASS", group=TRUE))
names(NBUD)
#soil total N
Nbud1_soiltotalN.lm <- lm(N_kgha~CLASS, NBUD)
Nbud1_soiltotalNageclassmeanscomp <- aov(N_kgha~CLASS, NBUD)
anova(Nbud1_soiltotalNageclassmeanscomp)
summary(Nbud1_soiltotalN.lm)
TukeyHSD(aov(N_kgha~CLASS, NBUD), ordered=T)
HSD.test(Nbud1_soiltotalNageclassmeanscomp, "CLASS", group=TRUE)
print(HSD.test(Nbud1_soiltotalNageclassmeanscomp, "CLASS", group=TRUE))
#No significant effects of age on soil total N,

#annualized soil Nmin
Nbud1_soilannualNmin.lm <- lm(AnnNmin50pctwdscale_kgha~CLASS, NBUD)
Nbud1_soilannualNminageclassmeanscomp <- aov(AnnNmin50pctwdscale_kgha~CLASS, NBUD)
anova(Nbud1_soilannualNminageclassmeanscomp)
summary(Nbud1_soilannualNmin.lm)
HSD.test(Nbud1_soilannualNminageclassmeanscomp, "CLASS", group=TRUE)
print(HSD.test(Nbud1_soilannualNminageclassmeanscomp, "CLASS", group=TRUE))
#No significant effects of age on annualized soil Nmin.
names(NBUD)

Nbud2_soilannualNmin_logage.lm <- lm(AnnNmin50pctwdscale_kgha~logage, NBUD)
summary(Nbud2_soilannualNmin_logage.lm)

Nbud3_soilannualNmin_logage.lm <- lm(AnnNmin50pctwdscale_kgha~logage, NBUD_nomat)
summary(Nbud3_soilannualNmin_logage.lm)

#Tree total N:
Nbud1_treetotalN.lm <- lm(log(TotalTreeN_kgha)~CLASS, NBUD)
Nbud1_treetotalNageclassmeanscomp <- aov(log(TotalTreeN_kgha)~CLASS, NBUD)
anova(Nbud1_treetotalNageclassmeanscomp)
#significant age effects
summary(Nbud1_treetotalN.lm)
HSD.test(Nbud1_treetotalNageclassmeanscomp, "CLASS", group=TRUE)
print(HSD.test(Nbud1_treetotalNageclassmeanscomp, "CLASS", group=TRUE))
#significant age effects and trend of increasing total tree N with site age.

#log age model-mature sites included
Nbud2_treetotalN_logage.lm <- lm(log(TotalTreeN_kgha)~logage, NBUD)
summary(Nbud2_treetotalN_logage.lm)

#log age model-mature sites excluded
Nbud3_treetotalN_logage.lm <- lm(log(TotalTreeN_kgha)~logage, NBUD_nomat)
summary(Nbud3_treetotalN_logage.lm)

#Tree N in woody biomass:
Nbud1_woodybiomassN.lm <- lm(log(WdBiomassN)~CLASS, NBUD)
Nbud1_woodybiomassageclassmeanscomp <- aov(log(WdBiomassN)~CLASS, NBUD)
anova(Nbud1_woodybiomassageclassmeanscomp)
summary(Nbud1_woodybiomassN.lm)
HSD.test(Nbud1_woodybiomassageclassmeanscomp, "CLASS", group=TRUE)
print(HSD.test(Nbud1_woodybiomassageclassmeanscomp, "CLASS", group=TRUE))

#log age model-mature sites included
Nbud2_woodybiomassN_logage.lm <- lm(log(WdBiomassN)~logage, NBUD)
summary(Nbud2_woodybiomassN_logage.lm)

#log age model-mature sites excluded
Nbud3_woodybiomassN_logage.lm <- lm(log(WdBiomassN)~logage, NBUD_nomat)
summary(Nbud3_woodybiomassN_logage.lm)

#N demand for canopy leaf biomass
Nbud1_canopyN.lm <- lm(CanopyN_kgha~CLASS, NBUD)
Nbud1_canopyNageclassmeanscomp <- aov(CanopyN_kgha~CLASS, NBUD)
anova(Nbud1_canopyNageclassmeanscomp)
summary(Nbud1_canopyN.lm)
HSD.test(Nbud1_canopyNageclassmeanscomp, "CLASS", group=TRUE)
print(HSD.test(Nbud1_canopyNageclassmeanscomp, "CLASS", group=TRUE))

Nbud1_canopyNlog.lm <- lm(log(CanopyN_kgha)~CLASS, NBUD)
Nbud1_canopyNageclassmeanscomplog <- aov(log(CanopyN_kgha)~CLASS, NBUD)
anova(Nbud1_canopyNageclassmeanscomplog)
summary(Nbud1_canopyNlog.lm)
print(HSD.test(Nbud1_canopyNageclassmeanscomplog, "CLASS", group=TRUE))
