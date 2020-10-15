
# Code used to examine climate, landscape, and butterfly trait effects on butterfly abundance and diversity trends

library(nlme)
library(MuMIn)
library(bbmle)
library(lme4)
library(car)
library(spdep)
library(viridis)
library(ggplot2)
library(DHARMa)
library(glmmTMB)

##############################################################
# Abundance trends with species*site (no pseudo-absences)

#####
# 6 traits, Jun-Aug set

data1 = read.table('./data/butterfly_traits_envars_trends_50km_NoMergeJunAug_m5_trim_6traits.txt',sep='\t',as.is=T,check.names=F,header=T); str(data1)
# Ensure categorical variables are factors
data1$Species = as.factor(data1$Species)
data1$Family = as.factor(data1$Family)
data1$NCGR = as.factor(data1$NCGR)
data1$LarvalColor = as.factor(data1$LarvalColor)
data1$LarvalHair = as.factor(data1$LarvalHair)
data1$AdultColor = as.factor(data1$AdultColor)
data1$grid_id = as.factor(data1$grid_id)

# Z-transform numeric covariates
data1$Cropland.2005_2015z = (data1$Cropland.2005_2015 - mean(data1$Cropland.2005_2015)) / sd(data1$Cropland.2005_2015)
data1$Built.2005_2015z = (data1$Built.2005_2015 - mean(data1$Built.2005_2015)) / sd(data1$Built.2005_2015)
data1$Cropland.trendz = (data1$Cropland.trend - mean(data1$Cropland.trend)) / sd(data1$Cropland.trend)
data1$Precip.1993_2018z = (data1$Precip.1993_2018 - mean(data1$Precip.1993_2018)) / sd(data1$Precip.1993_2018)
data1$Temp.1993_2018z = (data1$Temp.1993_2018 - mean(data1$Temp.1993_2018)) / sd(data1$Temp.1993_2018)
data1$Precip.trendz = (data1$Precip.trend - mean(data1$Precip.trend)) / sd(data1$Precip.trend)
data1$Temp.trendz = (data1$Temp.trend - mean(data1$Temp.trend)) / sd(data1$Temp.trend)
data1$Diet.breadth.familiesz = (data1$Diet.breadth.families - mean(data1$Diet.breadth.families)) / sd(data1$Diet.breadth.families)
data1$AdultSizez = (data1$AdultSize - mean(data1$AdultSize)) / sd(data1$AdultSize)

# Identify best spatial autocorrelation structure  https://stats.idre.ucla.edu/r/faq/how-do-i-model-a-spatially-autocorrelated-outcome/
mod1 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1)
mod.gau <- update(mod1, correlation = corGaus(1, form = ~ Lon + Lat), method = "ML")
mod.exp <- update(mod1, correlation = corExp(1, form = ~ Lon + Lat), method = "ML")
mod.spher <- update(mod1, correlation = corSpher(1, form = ~ Lon + Lat), method = "ML")
AIC(mod.gau,mod.exp,mod.spher) #exp lowest

# Hand-build models to find best model by parsimony
# Full +/- traits
mod1 = lme(Abundance.trend ~ 1, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod2 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod3 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod4 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod5 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod6 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod7 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod8 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod9 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod10 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod11 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod12 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# No traits, leave one env out
mod13 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod14 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod15 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod16 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod17 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod18 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# Landscape only
mod19 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z , random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod20 = lme(Abundance.trend ~ Cropland.trendz + Built.2005_2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# Climate only
mod21 = lme(Abundance.trend ~ Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod22 = lme(Abundance.trend ~ Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod23 = lme(Abundance.trend ~ Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod24 = lme(Abundance.trend ~ Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# All diffs +/- traits
mod25 = lme(Abundance.trend ~ Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod26 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod27 = lme(Abundance.trend ~ AdultSizez + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod28 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod29 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod30 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod31 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# One at a time
mod32 = lme(Abundance.trend ~ AdultSizez, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod33 = lme(Abundance.trend ~ Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod34 = lme(Abundance.trend ~ Cropland.2005_2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod35 = lme(Abundance.trend ~ NCGR, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod36 = lme(Abundance.trend ~ LarvalColor, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod37 = lme(Abundance.trend ~ LarvalHair, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod38 = lme(Abundance.trend ~ AdultColor, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod39 = lme(Abundance.trend ~ Cropland.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod40 = lme(Abundance.trend ~ Built.2005_2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod41 = lme(Abundance.trend ~ Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod42 = lme(Abundance.trend ~ Precip.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod43 = lme(Abundance.trend ~ Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod44 = lme(Abundance.trend ~ Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# Traits only
mod45 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod46 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod47 = lme(Abundance.trend ~ LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod48 = lme(Abundance.trend ~ NCGR + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod49 = lme(Abundance.trend ~ NCGR + LarvalColor + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod50 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod51 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod52 = lme(Abundance.trend ~ NCGR + AdultColor + AdultSizez, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# All contemporary +/- traits
mod53 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod54 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod55 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod56 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod57 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod58 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod59 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod60 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod61 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod62 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod63 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod64 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# No traits or diffs, leave one env out
mod65 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod66 = lme(Abundance.trend ~ Cropland.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod67 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod68 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod69 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod70 = lme(Abundance.trend ~ Cropland.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# No traits or contemporary, leave one env out
mod71 = lme(Abundance.trend ~ Cropland.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod72 = lme(Abundance.trend ~ Cropland.trendz + Precip.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')

AICctab(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,
	mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43,
	mod44,mod45,mod46,mod47,mod48,mod49,mod50,mod51,mod52,mod53,mod54,mod55,mod56,mod57,mod58,mod59,mod60,mod61,mod62,mod63,mod64,
	mod65,mod66,mod67,mod68,mod69,mod70,mod71,mod72,weights=TRUE)
modav = model.avg(list(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,
	mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43,
	mod44,mod45,mod46,mod47,mod48,mod49,mod50,mod51,mod52,mod53,mod54,mod55,mod56,mod57,mod58,mod59,mod60,mod61,mod62,mod63,mod64,
	mod65,mod66,mod67,mod68,mod69,mod70,mod71,mod72))
confint(modav, level = 0.90)

mod13 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod13)
mod3 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod3)
mod4 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod4)
mod17 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod17)
mod2 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod2)
mod5 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod5)

# For extrapolation map
mod13a = lme(Abundance.trend ~ Built.2005_2015 + Precip.1993_2018*Precip.trend + Temp.1993_2018*Temp.trend, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='REML')
mod13a


#####
# 6 traits, Jul set

data1 = read.table('./data/butterfly_traits_envars_trends_50km_NoMergeJul_m5_trim_6traits.txt',sep='\t',as.is=T,check.names=F,header=T); str(data1)
# Ensure categorical variables are factors
data1$Species = as.factor(data1$Species)
data1$Family = as.factor(data1$Family)
data1$NCGR = as.factor(data1$NCGR)
data1$LarvalColor = as.factor(data1$LarvalColor)
data1$LarvalHair = as.factor(data1$LarvalHair)
data1$AdultColor = as.factor(data1$AdultColor)
data1$grid_id = as.factor(data1$grid_id)

# Z-transform numeric covariates
data1$Cropland.2005_2015z = (data1$Cropland.2005_2015 - mean(data1$Cropland.2005_2015)) / sd(data1$Cropland.2005_2015)
data1$Built.2005_2015z = (data1$Built.2005_2015 - mean(data1$Built.2005_2015)) / sd(data1$Built.2005_2015)
data1$Cropland.trendz = (data1$Cropland.trend - mean(data1$Cropland.trend)) / sd(data1$Cropland.trend)
data1$Precip.1993_2018z = (data1$Precip.1993_2018 - mean(data1$Precip.1993_2018)) / sd(data1$Precip.1993_2018)
data1$Temp.1993_2018z = (data1$Temp.1993_2018 - mean(data1$Temp.1993_2018)) / sd(data1$Temp.1993_2018)
data1$Precip.trendz = (data1$Precip.trend - mean(data1$Precip.trend)) / sd(data1$Precip.trend)
data1$Temp.trendz = (data1$Temp.trend - mean(data1$Temp.trend)) / sd(data1$Temp.trend)
data1$Diet.breadth.familiesz = (data1$Diet.breadth.families - mean(data1$Diet.breadth.families)) / sd(data1$Diet.breadth.families)
data1$AdultSizez = (data1$AdultSize - mean(data1$AdultSize)) / sd(data1$AdultSize)

# Identify best spatial autocorrelation structure  https://stats.idre.ucla.edu/r/faq/how-do-i-model-a-spatially-autocorrelated-outcome/
mod1 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1)
mod.gau <- update(mod1, correlation = corGaus(1, form = ~ Lon + Lat), method = "ML")
mod.exp <- update(mod1, correlation = corExp(1, form = ~ Lon + Lat), method = "ML")
mod.spher <- update(mod1, correlation = corSpher(1, form = ~ Lon + Lat), method = "ML")
AIC(mod.gau,mod.exp,mod.spher) #exp lowest

# Hand-build models to find best model by parsimony
# Full +/- traits
mod1 = lme(Abundance.trend ~ 1, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod2 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod3 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod4 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod5 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod6 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod7 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod8 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod9 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod10 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod11 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod12 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# No traits, leave one env out
mod13 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod14 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod15 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod16 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod17 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod18 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# Landscape only
mod19 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z , random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod20 = lme(Abundance.trend ~ Cropland.trendz + Built.2005_2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# Climate only
mod21 = lme(Abundance.trend ~ Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod22 = lme(Abundance.trend ~ Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod23 = lme(Abundance.trend ~ Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod24 = lme(Abundance.trend ~ Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# All diffs +/- traits
mod25 = lme(Abundance.trend ~ Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod26 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod27 = lme(Abundance.trend ~ AdultSizez + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod28 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod29 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod30 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod31 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + Diet.breadth.familiesz + Cropland.trendz + Precip.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# One at a time
mod32 = lme(Abundance.trend ~ AdultSizez, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod33 = lme(Abundance.trend ~ Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod34 = lme(Abundance.trend ~ Cropland.2005_2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod35 = lme(Abundance.trend ~ NCGR, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod36 = lme(Abundance.trend ~ LarvalColor, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod37 = lme(Abundance.trend ~ LarvalHair, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod38 = lme(Abundance.trend ~ AdultColor, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod39 = lme(Abundance.trend ~ Cropland.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod40 = lme(Abundance.trend ~ Built.2005_2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod41 = lme(Abundance.trend ~ Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod42 = lme(Abundance.trend ~ Precip.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod43 = lme(Abundance.trend ~ Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod44 = lme(Abundance.trend ~ Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# Traits only
mod45 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod46 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod47 = lme(Abundance.trend ~ LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod48 = lme(Abundance.trend ~ NCGR + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod49 = lme(Abundance.trend ~ NCGR + LarvalColor + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod50 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod51 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod52 = lme(Abundance.trend ~ NCGR + AdultColor + AdultSizez, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# All contemporary +/- traits
mod53 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod54 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod55 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod56 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod57 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod58 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod59 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod60 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod61 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod62 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod63 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod64 = lme(Abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# No traits or diffs, leave one env out
mod65 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod66 = lme(Abundance.trend ~ Cropland.2005_2015z + Precip.1993_2018z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod67 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod68 = lme(Abundance.trend ~ Cropland.2005_2015z + Built.2005_2015z + Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod69 = lme(Abundance.trend ~ Built.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod70 = lme(Abundance.trend ~ Cropland.2005_2015z + Precip.1993_2018z + Temp.1993_2018z + Temp.1993_2018z*Precip.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
# No traits or contemporary, leave one env out
mod71 = lme(Abundance.trend ~ Cropland.trendz + Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
mod72 = lme(Abundance.trend ~ Cropland.trendz + Precip.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')

AICctab(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,
	mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43,
	mod44,mod45,mod46,mod47,mod48,mod49,mod50,mod51,mod52,mod53,mod54,mod55,mod56,mod57,mod58,mod59,mod60,mod61,mod62,mod63,mod64,
	mod65,mod66,mod67,mod68,mod69,mod70,mod71,mod72,weights=TRUE)
modav = model.avg(list(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,
	mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43,
	mod44,mod45,mod46,mod47,mod48,mod49,mod50,mod51,mod52,mod53,mod54,mod55,mod56,mod57,mod58,mod59,mod60,mod61,mod62,mod63,mod64,
	mod65,mod66,mod67,mod68,mod69,mod70,mod71,mod72))
confint(modav, level = 0.90)

mod5 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod5)
mod3 = lme(Abundance.trend ~ AdultSizez + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod3)
mod4 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod4)
mod9 = lme(Abundance.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Precip.1993_2018z*Temp.1993_2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat), method='ML')
summary(mod9)


###################################################
# GLS with richness, and evenness (per site) - minimum sample of 100 butterflies

library(nlme)
library(MuMIn)
library(bbmle)
library(lme4)
library(car)
library(spdep)
library(viridis)
library(ggplot2)
library(DHARMa)
library(glmmTMB)

###
# JunAug dataset

# Import data
con = file("./data/butterfly_diversity_envars_persite_JunAug_min100.txt",'r')
dat = readLines(con)
close(con)
data1 = c()
for (i in 1:length(dat)){
	d1 = strsplit(dat[i],'\t')[[1]]
	data1 = data.frame(rbind(data1,d1),stringsAsFactors=F)
}
colnames(data1) = data1[1,]
data1 = data1[-1,]
data1$Rarefied.richness.trend = as.numeric(data1$Rarefied.richness.trend)
data1 = data1[which(!is.na(data1$Rarefied.richness.trend)),]
data1$Rarefied.evenness.trend = as.numeric(data1$Rarefied.evenness.trend)
data1$Cropland.2005_2015 = as.numeric(data1$Cropland.2005_2015)
data1$Cropland.trend = as.numeric(data1$Cropland.trend)
data1$Built.2005_2015 = as.numeric(data1$Built.2005_2015)
data1$Precip.1993_2018 = as.numeric(data1$Precip.1993_2018)
data1$Precip.trend = as.numeric(data1$Precip.trend)
data1$Temp.1993_2018 = as.numeric(data1$Temp.1993_2018)
data1$Temp.trend = as.numeric(data1$Temp.trend)
data1$Lon = as.numeric(data1$Lon)
data1$Lat = as.numeric(data1$Lat)
summary(data1)

# Z-transform numeric covariates
data1$Cropland.2005_2015z = (data1$Cropland.2005_2015 - mean(data1$Cropland.2005_2015)) / sd(data1$Cropland.2005_2015)
data1$Built.2005_2015z = (data1$Built.2005_2015 - mean(data1$Built.2005_2015)) / sd(data1$Built.2005_2015)
data1$Cropland.trendz = (data1$Cropland.trend - mean(data1$Cropland.trend)) / sd(data1$Cropland.trend)
data1$Precip.1993_2018z = (data1$Precip.1993_2018 - mean(data1$Precip.1993_2018)) / sd(data1$Precip.1993_2018)
data1$Temp.1993_2018z = (data1$Temp.1993_2018 - mean(data1$Temp.1993_2018)) / sd(data1$Temp.1993_2018)
data1$Precip.trendz = (data1$Precip.trend - mean(data1$Precip.trend)) / sd(data1$Precip.trend)
data1$Temp.trendz = (data1$Temp.trend - mean(data1$Temp.trend)) / sd(data1$Temp.trend)

# Richness (per site)
data.spatialCor.glsSpher = gls(Rarefied.richness.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z:Precip.1993_2018z, data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="ML")
dd <- dredge(data.spatialCor.glsSpher) # takes a while
subset(dd, delta < 2)
plot(1,1,col='white',xaxt='n',yaxt='n',ann=F,bty='n'); legend('topleft',legend='Done!',cex=3)
# Intercept-only is the best model
#write.csv(dd,"AICcdredge_richness_rarefied_persite_JunAug.csv", row.names = FALSE)
write.csv(dd,"AICcdredge_richness_rarefied_persite_JunAug_min100.csv", row.names = FALSE)
# Best-supported model
mod1 = gls(Rarefied.richness.trend ~ Temp.trendz, data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")

# Evenness (per site)
data.spatialCor.glsSpher = gls(Rarefied.evenness.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z:Precip.1993_2018z, data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="ML")
dd <- dredge(data.spatialCor.glsSpher) # takes a while
subset(dd, delta < 2)
plot(1,1,col='white',xaxt='n',yaxt='n'); legend('topleft',legend='Done!',cex=3)
#write.csv(dd,"AICcdredge_evenness_rarefied_persite_JunAug.csv", row.names = FALSE)
write.csv(dd,"AICcdredge_evenness_rarefied_persite_JunAug_min100.csv", row.names = FALSE)
# Best-supported model
mod1 = gls(Rarefied.evenness.trend ~ Cropland.2005_2015z + Precip.1993_2018z + Precip.trendz + Temp.trendz,data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")


###
# Jul dataset

library(nlme)
library(MuMIn)
library(bbmle)
library(lme4)
library(car)
library(spdep)
library(viridis)
library(ggplot2)
library(DHARMa)
library(glmmTMB)

# Import data
con = file("./data/butterfly_diversity_envars_persite_Jul_min100.txt",'r')
dat = readLines(con)
close(con)
data1 = c()
for (i in 1:length(dat)){
	d1 = strsplit(dat[i],'\t')[[1]]
	data1 = data.frame(rbind(data1,d1),stringsAsFactors=F)
}
colnames(data1) = data1[1,]
data1 = data1[-1,]
data1$Rarefied.richness.trend = as.numeric(data1$Rarefied.richness.trend)
data1 = data1[which(!is.na(data1$Rarefied.richness.trend)),]
data1$Rarefied.evenness.trend = as.numeric(data1$Rarefied.evenness.trend)
data1$Cropland.2005_2015 = as.numeric(data1$Cropland.2005_2015)
data1$Cropland.trend = as.numeric(data1$Cropland.trend)
data1$Built.2005_2015 = as.numeric(data1$Built.2005_2015)
data1$Precip.1993_2018 = as.numeric(data1$Precip.1993_2018)
data1$Precip.trend = as.numeric(data1$Precip.trend)
data1$Temp.1993_2018 = as.numeric(data1$Temp.1993_2018)
data1$Temp.trend = as.numeric(data1$Temp.trend)
data1$Lon = as.numeric(data1$Lon)
data1$Lat = as.numeric(data1$Lat)
summary(data1)

# Z-transform numeric covariates
data1$Cropland.2005_2015z = (data1$Cropland.2005_2015 - mean(data1$Cropland.2005_2015)) / sd(data1$Cropland.2005_2015)
data1$Built.2005_2015z = (data1$Built.2005_2015 - mean(data1$Built.2005_2015)) / sd(data1$Built.2005_2015)
data1$Cropland.trendz = (data1$Cropland.trend - mean(data1$Cropland.trend)) / sd(data1$Cropland.trend)
data1$Precip.1993_2018z = (data1$Precip.1993_2018 - mean(data1$Precip.1993_2018)) / sd(data1$Precip.1993_2018)
data1$Temp.1993_2018z = (data1$Temp.1993_2018 - mean(data1$Temp.1993_2018)) / sd(data1$Temp.1993_2018)
data1$Precip.trendz = (data1$Precip.trend - mean(data1$Precip.trend)) / sd(data1$Precip.trend)
data1$Temp.trendz = (data1$Temp.trend - mean(data1$Temp.trend)) / sd(data1$Temp.trend)

# Richness (per site)
data.spatialCor.glsSpher = gls(Rarefied.richness.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z,
	data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="ML")
dd <- dredge(data.spatialCor.glsSpher) # takes a while
subset(dd, delta < 2)
plot(1,1,col='white',xaxt='n',yaxt='n'); legend('topleft',legend='Done!',cex=3)
# Intercept-only is the best model
write.csv(dd,"AICcdredge_richness_rarefied_persite_Jul_min100.csv", row.names = FALSE)

# Evenness (per site)
data.spatialCor.glsSpher = gls(Rarefied.evenness.trend ~ Cropland.2005_2015z*Cropland.trendz + Built.2005_2015z + Precip.1993_2018z*Precip.trendz + Temp.1993_2018z*Temp.trendz + Temp.1993_2018z*Precip.1993_2018z,
	data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="ML")
dd <- dredge(data.spatialCor.glsSpher) # takes a while
subset(dd, delta < 2)
plot(1,1,col='white',xaxt='n',yaxt='n'); legend('topleft',legend='Done!',cex=3)
write.csv(dd,"AICcdredge_evenness_rarefied_persite_Jul_min100.csv", row.names = FALSE)
