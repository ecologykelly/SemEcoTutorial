###############################################################
#  Sudnick, Brodie & Williams
#manuscript ECE-2021-05-00859 entitled "Nature versus nurture: Structural equation modeling indicates that parental care does not mitigate consequences of poor environmental conditions"
# Code for Structural equation modeling (SEM)
# Author: Kelly A. Williams (KAW)
###############################################################

# Install dependencies if you don't have them

install.packages("semPlot", dependencies = TRUE)
install.packages("piecewiseSEM", dependencies = TRUE)

# load packages (install if you don't have them)
library(semPlot)
library(readr)
library(tidyr)
library(dplyr)
library(piecewiseSEM)
library(lme4)
library(nlme)
library(car)
##############################################################

# load data

full <- read_csv("data/full.csv")


#### data prep ####

full$Blowflies<-factor(full$Blowflies)
levels(factor(full$Blowflies)) # check levels and make consistent
levels(full$Blowflies)[levels(full$Blowflies)=="NO"]<-"No"
levels(full$Blowflies)[levels(full$Blowflies)=="YES"]<-"Yes"
levels(factor(full$Blowflies)) # check that now two levels No and Yes

#make bf a binary 0 1 for blowflies. 0 = No, 1 = yes
full$bf<-full$Blowflies
levels(full$bf)[levels(full$bf)=="No"] <- 0
levels(full$bf)[levels(full$bf)=="Yes"] <- 1
levels(factor(full$bf)) # check levels
full$Nestl.ID<-factor(full$Nestl.ID)
length(unique(full$Nestl.ID)) #348 nestlings


#### tidy data for use in analyses ####
#Because repeated measures of nestlings, need to get one row of data per nestling. Several steps
#First group data set by nestling ID.
by.nestling<-full%>%group_by(Nestl.ID) # multiple observations per nestling
# Next need one row of data per nestling for variables that will be used in analyses.
#Get mean values for the following variables
nestling<-by.nestling%>%
  summarize_at(c("Vol.2","mass.A","Total.Trip.min", "T.perc","nestling.temp"), mean, na.rm=TRUE)

# check and compare the two data frames
head(nestling) # includes one line of data per nestling (n = 348) with only variables requested
head(by.nestling) #includes all variables with multiple lines per nestling

#Next use the following function to pull out 1 value per nestling for bf(blow fly presence absence),k (growth rate), hem (hematocrit) since these are the same values on multiple lines
nest.calc = function(nestl){
  val<-by.nestling%>%filter(Nestl.ID == nestl)%>% # filter by nestling by.nestling data
    select(Nestl.ID,bf,K,Hem,Nest.ID, Nestling.Max, Species.J, Vol.2)
  bf.k.hem<-data.frame("Nestl.ID"=val[1,1], "bf"=val[1,2],"K"=val[1,3],"Hem"=val[1,4],"Nest.ID"=val[1,5], "Nestling.Max"=val[1,6], "Species.J"=val[1,7], "Vol.2"=val[1,8]) # only keep 1 value
}

list.bf.k.hem<-lapply(nestling$Nestl.ID, nest.calc) #apply function to get list of data frames
#The above will take a minute.
bf.k.hem<-bind_rows(list.bf.k.hem) # merge list of data frames --> n = 348
head(bf.k.hem)
nestling1<-cbind(nestling, bf.k.hem) # data set with one row per nestling with 14 variables


#### subset data to get only EABL ####
nestling1$Species<-factor(bf.k.hem$Species.J) # Species needs to be a factor to subset data
levels(factor(nestling1$Species)) #note TRES data are included
EABL<-subset(nestling1, Species=="EABL") #235 obs
str(EABL)


#### name variables for path model ####
Food<-as.numeric(EABL$Vol.2) # Vol.2 is volume of arthropods 1-10mm
bf<-EABL$bf # will need to make numeric below
Temp<-as.numeric(EABL$nestling.temp)
Prov<-as.numeric(EABL$Total.Trip.min)
Att<-as.numeric(EABL$T.perc)
Growth<-as.numeric(EABL$K)
Nest.ID<-factor(EABL$Nest.ID)
Biomass<-EABL$mass.A  # biomass of arthropods
Brood<-EABL$Nestling.Max # maximum brood size recorded
Hem<-EABL$Hem

#create data frame with the following variables
dat<-data.frame(Food,Biomass, Temp, Prov, Att, Growth, Nest.ID,bf,Brood, Hem)
str(dat) # check
dat$Nestl.ID<-EABL$Nestl.ID #  add nestling ID to data frame
dat$bf<-as.numeric(EABL$bf) ##bf as numeric #problem in piecewiseSEM with bf as factor so made it numeric. Source for suggestion: jslefche.github.io/sem_book/categorical
str(dat) #235 obs

summary(lm(Growth~Hem, data=dat)) # check relationship between growth and hematocrit. note not repeated measure by nest so model not valid - just a check

#check for missing data and check sample sizes for different variables
table(is.na(dat$Growth)) #62 missing 235-62 = 173
table(is.na(dat$Hem)) # 13 missing, 235-13=222
table(is.na(dat$Biomass)) #26 missing biomass so n = 209


#### K (growth) data with biomass####
head(dat)
k<-dat[,1:9] # exclude Hem Gives n = 235
head(k)
k.dat<-na.omit(k) # n = 149
cov(k.dat[,c(1,2,3,4,5,6)]) # variance/covariance
table(k.dat$Nest.ID) # check number of nestlings per nest. Need to use nest.ID as random variable! can't use Lavaan. need to use piecewiseSEM
length(unique(k.dat$Nest.ID)) # 40 nests

#### data with hem without growth####
hem.dat<-dat[,c(1:5,7:10)] # n = 235
head(hem.dat)
hem.dat<-na.omit(hem.dat) # n = 179
length(unique(hem.dat$Nest.ID)) # 49 nests


#### use random = nest ID with piecewiseSEM!!!####

#run full glm and check for multicollinearity (all variables were already standardized)
mod<-lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat)
vif(mod) # all vif<2 so good.

##### psem on EABL growth rates ####
#Make a psem model with no correlations use k.dat
fit.psem<-psem(
  lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat)
  )

summary(fit.psem) # model output
fisherC(fit.psem) # fit. p = 0 so not a good fit
AIC(fit.psem)
#conditional independence claims to be used in determining the goodness-of-fit for piecewise structural equation models indicate we might need Att~Prov and need Temp~Prov

# add correlations
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

summary(fit.psem1) # Food not sig with growth but is with temp
fisherC(fit.psem1)
AIC(fit.psem1)
plot(fit.psem1, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

# above look at correlated errors: att~ prov and prov~brood not correlated. Remove and refit
fit.psem1.1<-psem(
  lme(Growth~Food+Biomass+ bf + Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Food%~~%Biomass
)

summary(fit.psem1.1) #
fisherC(fit.psem1.1)
AIC(fit.psem1.1)

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
summary(fit.psem2) # Need Growth ~ biomass? marginal sig test of directed selection
fisherC(fit.psem2)
AIC(fit.psem2)

# add Growth~biomass back in
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
summary(fit.psem3) #
fisherC(fit.psem3)
AIC(fit.psem3)
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
summary(fit.psem4) #also checked but found no change when Growth%~~%Brood so left out
fisherC(fit.psem4)
AIC(fit.psem4)
rsquared(fit.psem4)
# look at paths
plot(fit.psem4, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05)

#drop growth~biomass
fit.psem4a<-psem(
  lme(Growth~bf +Prov+Att +Temp, random = ~1|Nest.ID, data=k.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  lme(Att~Food+Biomass, random = ~1|Nest.ID, data=k.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=k.dat),
  Prov%~~%Temp,
  Temp%~~%Brood,
  Att%~~%Prov,
  Prov%~~%Brood,
  Food%~~%Biomass
)
summary(fit.psem4a) #also checked but found no change when Growth%~~%Brood so left out
fisherC(fit.psem4a)
AIC(fit.psem4a)


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
summary(fit.psem5) #  need to include att & prov
fisherC(fit.psem5)
AIC(fit.psem5)

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
summary(fit.psem6) #
fisherC(fit.psem6)
AIC(fit.psem6)



#### models with hematocrit use hem.dat  (Aug 2021 n = 179; original was 135 obs) #####
fit.hem<-psem(
  lme(Hem~bf +Food + Biomass+ Prov+Att +Temp, random = ~1|Nest.ID, data=hem.dat),
  lme(Prov~Food+Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Att~Food+Biomass+bf, random = ~1|Nest.ID, data=hem.dat),
  lme(Temp~Food+Biomass+bf, random = ~1|Nest.ID, data=hem.dat))
summary(fit.hem) #
AIC(fit.hem)
fisherC(fit.hem)


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
summary(fit.hem.psem) #
AIC(fit.hem.psem)
fisherC(fit.hem.psem)
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
summary(fit.hem.psem1) #
#use food instead of biomass?
fisherC(fit.hem.psem1)
AIC(fit.hem.psem1)


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
summary(fit.hem.psem1.1) #
fisherC(fit.hem.psem1.1)
AIC(fit.hem.psem1.1)

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
summary(fit.hem.psem2) #
AIC(fit.hem.psem2)
fisherC(fit.hem.psem2)

#if remove Att~Prov worse aic 62.68.
# keep correlations same. For fixed effects, remove Hem~ Att & Temp
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
summary(fit.hem.psem2.1) #
fisherC(fit.hem.psem2.1)
AIC(fit.hem.psem2.1)

# remove food
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
summary(fit.hem.psem2.2) #
fisherC(fit.hem.psem2.2)
AIC(fit.hem.psem2.2)

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
summary(fit.hem.psem2.3) #
AIC(fit.hem.psem2.3)
fisherC(fit.hem.psem2.3)
plot(fit.hem.psem2.3, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05) # biomass sig, near 0.09 neg with prov.

# att Temp~~ATT
#### Best Hem model####
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
summary(fit.hem.psem2.3.1) # Checked with hem~Temp - no effect and model a bit worse
fisherC(fit.hem.psem2.3.1)
AIC(fit.hem.psem2.3.1)
rsquared(fit.hem.psem2.3.1)

plot(fit.hem.psem2.3.1, node_attrs = list(shape="rectangle",color="black", show="std", x = 4), ns_dashed=TRUE, alpha=0.05) #


