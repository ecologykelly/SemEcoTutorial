###############################################################
#  Sudnick, Brodie & Williams
#manuscript ECE-2021-05-00859 entitled "Nature versus nurture: Structural equation modeling indicates that parental care does not mitigate consequences of poor environmental conditions"
# Code for Structural equation modeling (SEM)
#  Author: Kelly A. Williams (KAW)
# test
###############################################################

# Packages needed

install.packages("semPlot", dependencies = TRUE)
install.packages("piecewiseSEM", dependencies = TRUE)

library(semPlot)
library(readr)
library(tidyr)
library(dplyr)
library(piecewiseSEM)
library(piecewiseSEM)
library(lme4)
library(nlme)
library(car)
##############################################################
#this needs fixed - use relative paths

#K laptop


#setwd("C:/Dropbox/kelly's shared files/Student Folders/Maddie/Thesis analysis")
full <- read_csv("Data Files/full.csv")
names(full)
View(full)
head(full$Survival)

#M laptop
setwd("C:/Users/mcsud/Dropbox/Maddie/Thesis analysis")
full <- read_csv("Data Files/full.csv")

#### data prep####
full$Blowflies<-factor(full$Blowflies)
levels(factor(full$Blowflies))
levels(full$Blowflies)[levels(full$Blowflies)=="NO"]<-"No"
levels(full$Blowflies)[levels(full$Blowflies)=="YES"]<-"Yes"
#make bf a binary 0 1 for blowflies. 0 = No, 1 = yes
full$bf<-full$Blowflies
levels(full$bf)[levels(full$bf)=="No"] <- 0
levels(full$bf)[levels(full$bf)=="Yes"] <- 1
levels(factor(full$bf))
full$Nestl.ID<-factor(full$Nestl.ID)
length(unique(full$Nestl.ID)) #348 nestlings
names(full)
levels(factor(full$Species.J))

# tidy data for use in analyses
#Because repeated measure by nestling need to get one row of data per nestling
by.nestling<-full%>%group_by(Nestl.ID) # group data set by nestling ID. multiple observations per nestling
#need one row of data per nestling. First get mean values for the following variables
nestling<-by.nestling%>%
  summarize_at(c("Vol.2","mass.A","Total.Trip.min", "T.perc","nestling.temp"), mean, na.rm=TRUE)

# check and compare the two data frames
head(nestling)
head(by.nestling)

#additional checks
nestling[15,] # check data for the 15th row
a.15<-subset(full, Nestl.ID=="160188558")
a.15$Vol.1
a.15$mass.A
a.15$Total.Trip.min # done checking

# use the following function to pull out 1 value per nestling for bf,k, hem. input 1 nestling ID.
nest.calc = function(nestl){
  val<-by.nestling%>%filter(Nestl.ID == nestl)%>% # filter by nestling by.nestling data
    select(Nestl.ID,bf,K,Hem,Nest.ID, Nestling.Max, Species.J, Vol.2)
  bf.k.hem<-data.frame("Nestl.ID"=val[1,1], "bf"=val[1,2],"K"=val[1,3],"Hem"=val[1,4],"Nest.ID"=val[1,5], "Nestling.Max"=val[1,6], "Species.J"=val[1,7], "Vol.2"=val[1,8]) # only keep 1 value
}

list.bf.k.hem<-lapply(nestling$Nestl.ID, nest.calc) #apply function to get list of data frames This will take a minute.
bf.k.hem<-bind_rows(list.bf.k.hem) # merge list of data frames
nestling<-cbind(nestling, bf.k.hem) # new data with 348 obs 22 variables
head(nestling)

by.nestling$K
str(nestling) #348 obs 22 vars
nestling$Species<-factor(nestling$Species.J)

#### EABL data (remove TRES data) ####
levels(factor(nestling$Species))

EABL<-subset(nestling, Species=="EABL") #235 obs of 23 vars variables
levels(factor(EABL$Species))
#some checks:
names(EABL)
head(EABL)
EABL$bf
str(EABL)

#rename variables
Food<-as.numeric(EABL$Vol.2) # Vol.2 is volume of arthropods 1-10mm
bf<-EABL$bf
Temp<-as.numeric(EABL$nestling.temp)
Prov<-as.numeric(EABL$Total.Trip.min)
Att<-as.numeric(EABL$T.perc)
Growth<-as.numeric(EABL$K)
Nest.ID<-EABL$Nest.ID
Biomass<-EABL$mass.A  # biomass of arthropods
Brood<-EABL$Nestling.Max # maximum brood size recorded
Hem<-EABL$Hem

#create data frame with the following variables
dat<-data.frame(Food,Biomass, Temp, Prov, Att, Growth, Nest.ID,bf,Brood, Hem)
str(dat) # check

dat$Nestl.ID<-EABL$Nestl.ID #  add nestling ID to data frame
dat$bf<-as.numeric(EABL$bf)##BF as numeric
str(dat)
head(dat) #235 obs

summary(lm(Growth~Hem, data=dat)) # check relationship between growth and hematocrit. note not repeated measure by nest
#K.Hem.lm<-lme(Growth~Hem, random = 1|Nest.ID, data = dat)

#check for missing data
table(is.na(dat$Growth)) #62 missing 235-62 = 173
table(is.na(dat$Hem)) # 13 missing, 235-13=222
table(is.na(dat$Biomass)) #26 missing biomass so n = 209
dat1<-na.omit(dat) #146 obs. if pull na from growth and hem, n = 136 if biomass, 148 with biomass exclude Hem
head(dat1)
str(dat1) #
plot(dat1) #
names(dat1)
str(dat1)
which(is.na(dat1)) # check
which(is.na(dat1$Att)) # check if still missing data - should not be
#dat1<-na.omit(dat1) # now 209
cov(dat1[,c(1,2,3,4,5)]) # covariance


#### K (growth) data with biomass####
k<-dat[,1:9] # exclude Hem Gives n = 235
head(dat)
head(k)
k.dat<-na.omit(k) #n = 149
names(k.dat)
table(k.dat$Nest.ID) # check number of nestlings per nest. Need to use nest.ID as random

#### K data without biomass####   #### DO I NEED TO KEEP THIS FOR PAPER?
k.food.dat<-dat[,c(1,3:9)] #235
head(k.food.dat)
k.food.dat<-na.omit(k.food.dat) # 159

#### data with hem without growth####
hem.dat<-dat[,c(1:5,7:10)]
head(hem.dat)
hem.dat<-na.omit(hem.dat) # n = 179


#### use random of nest ID try piecewiseSEM####

str(k.dat) # #error with bf as factor so made it numeric above (should be numeric here). Source for suggestionL jslefche.github.io/sem_book/categorical suggests running as numeric
names(k.dat)

#run full glm and check for multicollinearity (remember all variables are standardized)
mod<-lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat)
vif(mod) # all vif<2 so good. later check to see if food and biomass should be used together

# Make a psem model. no correlations
fit.psem<-psem(
  lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat)
  )

summary(fit.psem) #AIC = 98.95, Cp=0, df=6. need to add Att~Prov, Temp~ Prov
fisherC(fit.psem)

fit.psem1<-psem(
  lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
  )

summary(fit.psem1) # AIC= 56.24, p=0.33 DF = 2 Fisher C= 2.244 p = 0.326, 2 df. Food not sig with growth but is with temp
fisherC(fit.psem1)
plot(fit.psem1, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

# look at correlated errors: att~ prov and prov~brood not correlated. Remove

fit.psem1.1<-psem(
  lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Food%~~%Biomass
)

summary(fit.psem1.1) #  aic = 61.15, Fisher's C = 7.152 p = .13, df = 4
fisherC(fit.psem1.1)


#put correlations back in but remove direct effect of food and biomass
fit.psem2<-psem(
  lme(Growth~bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
  )
summary(fit.psem2) #fit C = 10.33, p = 0.11, Df = 6 aic= 60.33, Need Growth ~ biomass?
fisherC(fit.psem2)

fit.psem3<-psem(
  lme(Growth~bf +Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.psem3) #AIc = 55.93, C = 3.93, p = .42, df = 4
fisherC(fit.psem3)
plot(fit.psem3, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

# remove att~bf
fit.psem4<-psem(
  lme(Growth~bf +Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.psem4) #Aic 54.84, C = 4.84 p = .56, df = 6  (no change when Growth%~~%Brood so left out)
fisherC(fit.psem4)
# look at paths
plot(fit.psem4, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

# remove att~prov
fit.psem5<-psem(
  lme(Growth~bf +Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.psem5) # aic 59.75, C= 9.75, p = 0.28 df = 8 need to include att & prov
fisherC(fit.psem5)
#remove Growth~att
fit.psem6<-psem(
  lme(Growth~bf +Biomass+ Prov+Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.psem6) # aic 55.25 fit C = 7.25 p 0.51 8df
fisherC(fit.psem6)
#### best above is fit.psem4 so modify####
# remove food
#### best model of K with biomass and food ####
fit.psem4.1<-psem(
  lme(Growth~bf +Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Prov%~~%Att,
  Prov%~~%Brood
 )
summary(fit.psem4.1) # aic 46.32, C= 2.32 p 0.677 df = 4, note: also changed Prov~Att. did not change outcome. when temp ~ food included significant but fit not good & aic 60s
fisherC(fit.psem4.1)
plot(fit.psem4.1, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)
rsquared(fit.psem4.1)

# temp not sig. drop
fit.psem4.2<-psem(
  lme(Growth~bf +Biomass+ Prov+Att, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Biomass, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Att,
  Prov%~~%Brood
)
summary(fit.psem4.2) # aic 32.4,C = 0.4 p 0.82 df = 2, when growth~brood, aic=72 p=0, df= 4 brood not relate to growth
fisherC(fit.psem4.2)
plot(fit.psem4.2, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

# not sure how biomass can directly affect growth. Remove. not better though
fit.psem4.2.2<-psem(
  lme(Growth~bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood
)
summary(fit.psem4.2.2) # aic 50.72 fit C = 8.72 p = .19 df = 6; directed sep test sGrowthay I need growth ~ biomass,
fisherC(fit.psem4.2.2)

# remove growth~Att
fit.psem4.1.1<-psem(
  lme(Growth~bf +Biomass+ Prov+Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Prov%~~%Att,
  Prov%~~%Brood
)
summary(fit.psem4.1.1) # aic 51.31, fit C = 9.31 p = 0.16, df = 6
fisherC(fit.psem4.1.1)


#### Run sem on k without biomass data = k.food.dat ####  DID NOT USE IN MODEL SELECTION. REMOVE???###
fit.kpsem<-psem(
  lme(Growth~bf +Food+ Prov+Att +Temp, random = ~1|Nest.ID, data=k.food.dat),
  lme(Prov~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  lme(Att~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  lme(Temp~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  Temp%~~%Prov,
  Prov%~~%Att,
  Att%~~%Temp,
  Prov%~~%Brood,
  Temp%~~%Brood,
  Growth%~~%Brood
)
summary(fit.kpsem) #AIC = 46, fit C = 0 p = 1 df = 0

#remove prov~bf
fit.kpsem1<-psem(
  lme(Growth~bf +Food+ Prov+Att +Temp, random = ~1|Nest.ID, data=k.food.dat),
  lme(Prov~Food, random = ~1|Nest.ID, data=k.food.dat),
  lme(Att~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  lme(Temp~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  Temp%~~%Prov,
  Prov%~~%Att,
  Att%~~%Temp,
  Prov%~~%Brood,
  Temp%~~%Brood,
  Growth%~~%Brood
)
summary(fit.kpsem1) #AIC = 44.23, fit C = 0.24 p = 0.9, df = 2



#remove att & growth,
fit.kpsem2<-psem(
  lme(Growth~bf +Food+ Prov+Temp, random = ~1|Nest.ID, data=k.food.dat),
  lme(Prov~Food, random = ~1|Nest.ID, data=k.food.dat),
  lme(Att~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  lme(Temp~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  Temp%~~%Prov,
  Prov%~~%Att,
  Att%~~%Temp,
  Prov%~~%Brood,
  Temp%~~%Brood,
  Growth%~~%Brood
)
summary(fit.kpsem2) #AIC = 42.87, fit C= 0.87, p = 0.93 df = 4

#remove growth~brood,
#### "best" model with food instead of biomass####
k.food.dat$bf.f<-factor(k.food.dat$bf) #. Prov still neg. dont use bf as numeric
fit.kpsem3<-psem(
  lme(Growth~bf +Food+ Prov+Temp, random = ~1|Nest.ID, data=k.food.dat),
  lme(Prov~Food, random = ~1|Nest.ID, data=k.food.dat),
  lme(Att~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  lme(Temp~Food+bf, random = ~1|Nest.ID, data=k.food.dat),
  Temp%~~%Prov,
  Prov%~~%Att,
  Att%~~%Temp,
  Prov%~~%Brood,
  Temp%~~%Brood
)
summary(fit.kpsem3) #AIC = 42.59, p = 0.96, df = 4
plot(fit.kpsem3, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)
rsquared(fit.kpsem3)


#### models with hem use hem.dat 135 obs  #####
fit.hem.psem<-psem(
  lme(Hem~bf +Food + Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.hem.psem) # AIC 54.5, fit p = 0.79, 2 df #aic = 54.44, p = 0.80, df = 2
plot(fit.hem.psem, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

# remove food
fit.hem.psem1<-psem(
  lme(Hem~bf + Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood
)
summary(fit.hem.psem1) # aic 46.5 but I think I want to keep food # aic = 46.18, p = 92, df = 2
#use food instead of biomss
fit.hem.psem1.1<-psem(
  lme(Hem~bf + Food+ Prov+Att +Temp, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+bf, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood
)
summary(fit.hem.psem1.1) # aic 46.6 fit p = 0.73 df = 2


# remove att~bf, temp~bf, prov~ bf
fit.hem.psem2<-psem(
  lme(Hem~bf +Food + Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.hem.psem2) # aic 54.96 fit p = 0.54, 8df # Aic = 52.03, p = 0.85, df = 8

#remove Att~Prov worse aic 62.68.  remove Hem~ Att & Temp
fit.hem.psem2.1<-psem(
  lme(Hem~bf +Food + Biomass+ Prov, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.hem.psem2.1) #aic 51.98, fit p = 0.79, df = 12

# romove food
fit.hem.psem2.2<-psem(
  lme(Hem~bf + Biomass+ Prov, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.hem.psem2.2) #aic 51, fit 0.83 df = 14

# remove bf
#### good but improved in 2.3.1####
fit.hem.psem2.3<-psem(
  lme(Hem~ Biomass+ Prov, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.hem.psem2.3) #aic 41.85 p = 0.98, df = 8 # 44.25, p = .83, df = 8
plot(fit.hem.psem2.3, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05) # biomass sig, near 0.09 neg with prov.

# att Temp~~ATT
#### Best Hem so far####
fit.hem.psem2.3.1<-psem(
  lme(Hem~ Biomass+ Prov, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass,
  Temp%~~%Att
)
summary(fit.hem.psem2.3.1) #aic 41.77, fit p = 0.94 df=6 Checked with hem~Temp - no effect and model a bit worse #Aic = 44.23, p = 0.65, df = 6
rsquared(fit.hem.psem2.3.1)

plot(fit.hem.psem2.3.1, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05) #

# remove prov~food
fit.hem.psem2.4<-psem(
  lme(Hem~ Biomass+ Prov, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Prov%~~%Brood,
  Att%~~%Prov,
  Food%~~%Biomass
)
summary(fit.hem.psem2.4) #aic 43.99, p = 0.816, df = 10



# add growth cor with hem (can't run hem and growth as Y)
fit.hem.psem.growth<-psem(
  lme(Hem~bf +Food + Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass, random = ~1|Nest.ID, data=hem.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass,
  Hem%~~%Growth
)
summary(fit.hem.psem.growth) # 54.96

#######################################################################################
#### TRES ####
levels(factor(nestling$Species))

TRES<-subset(nestling, Species=="TRES") #113 obs of 19 vars variables
levels(factor(TRES$Species))
str(TRES)
#rename variables into new data
tres<-c()
tres$Food<-as.numeric(TRES$Vol.1)
tres$bf<-as.numeric(TRES$bf)
tres$Temp<-as.numeric(TRES$nestling.temp)
tres$Prov<-as.numeric(TRES$Total.Trip.min)
tres$Att<-as.numeric(TRES$T.perc)
tres$Growth<-as.numeric(TRES$K)
tres$Nest.ID<-TRES$Nest.ID
tres$Biomass<-TRES$mass.A
tres$Brood<-TRES$Nestling.Max
tres$Hem<-TRES$Hem
str(tres)
tres<-data.frame(tres)
head(tres)
str(tres)  #113 obs
table(tres$Nest.ID) # 30 nests
table(is.na(tres$Growth)) # 61 false, 52 true
table(is.na(tres$Hem)) # 62 false, 51 true
table(is.na(tres$Food)) # 22 missing
table(is.na(tres$Biomass)) # 30 missing
head(tres.dat)
tres.dat<-na.omit(tres) #50 obs
str(tres.dat)
table(tres.dat$Nest.ID) # only 12 nests
tres.k<-tres[,c(1:7,9)] # 113
tres.k<-na.omit(tres.k) #50
tres.hem<-tres[,c(1:5,6:9)] #
tres.hem<-na.omit(tres.hem) #50
names(tres.dat)
head(tres.dat)
str(tres.dat)
tres.dat$Nest.ID<-factor(tres.dat$Nest.ID)
#### model Growth####
# model singular with food and biomass. start with
summary(lme(Growth~Biomass, random = ~1|Nest.ID, data=tres.dat))

fit.k.tres.psem<-psem(
  lme(Growth~Prov,random = ~1|Nest.ID, data=tres.dat),
  lme(Prov~Biomass, random = ~1|Nest.ID, data=tres.dat),
  Growth%~~%Brood
)
summary(fit.k.tres.psem) #aic 17 no direct effect has 2 df, fit p = 0.579
plot(fit.k.tres.psem, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

#remove correlation
fit.k.tres.psem1<-psem(
  lme(Growth~Prov,random = ~1|Nest.ID, data=tres.dat),
  lme(Prov~Biomass, random = ~1|Nest.ID, data=tres.dat)
)
summary(fit.k.tres.psem1) # biomass not right

#try food
fit.k.tres.psem2<-psem(
  lme(Growth~Prov+Food,random = ~1|Nest.ID, data=tres.dat),
  lme(Prov~Food, random = ~1|Nest.ID, data=tres.dat)
)
summary(fit.k.tres.psem2) # aic 18 fit p = 1, 0 df
plot(fit.k.tres.psem2, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

#remove growth~food
fit.k.tres.psem3<-psem(
  lme(Growth~Prov,random = ~1|Nest.ID, data=tres.dat),
  lme(Prov~Food, random = ~1|Nest.ID, data=tres.dat)
)
summary(fit.k.tres.psem3) # aic 17.4 fit p = 0.486, 2 df
plot(fit.k.tres.psem3, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

#add Hem. doesn't work with Growth~Hem try correlation
fit.k.tres.psem4<-psem(
  lme(Growth~Prov,random = ~1|Nest.ID, data=tres.dat),
  lme(Prov~Food, random = ~1|Nest.ID, data=tres.dat),
  Growth%~~%Hem
)
summary(fit.k.tres.psem4) # aic 17.4, fit = 0.486, df=2
plot(fit.k.tres.psem4, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

tres.dat$bf # all 1.
