
setwd('<your/working/directory>')


# Merge trends with covariates

data0 = read.csv("./data/NABA_trends_envars_50km_wNAs.csv",as.is=T,check.names=F,header=TRUE)

# Add butterfly trait data
trends.summary = read.table('./butterfly_trends_summary_1993-2018_50km_wNAs.txt',sep='\t',as.is=T,check.names=F,header=T)
trends = trends.summary[,c(1,15,19:23)]
colnames(trends)[which(colnames(trends)=='mdn_Tau')] = 'Median.abundance.trend'
traits = read.table('./data/NABA_traits.txt',sep='\t',as.is=T,check.names=F,header=T); str(traits)
traits$Habitat[which(traits$Habitat=='')] = NA
traits$NCGR[which(traits$NCGR=='')] = NA
traits$OverwinteringStage[which(traits$OverwinteringStage=='')] = NA
traits$Clutches[which(traits$Clutches=='')] = NA
traits$LarvalColor[which(traits$LarvalColor=='')] = NA
traits$AdultColor[which(traits$AdultColor=='')] = NA
traits2 = data.frame('Species' = paste(traits$Genus,traits$Species,sep=' '),
	'Family' = as.factor(traits$Family),
	'Subfamily' = as.factor(traits$Subfamily),
	'Habitat' = as.factor(traits$Habitat),
	'NCGR' = as.factor(traits$NCGR),
	'OverwinteringStage' = as.factor(traits$OverwinteringStage),
	'MinSouthBroodsPerYear' = as.factor(traits$MinSouthBroodsPerYear),
	'MinNorthBroodsPerYear' = as.factor(traits$MinNorthBroodsPerYear),
	'Clutches' = as.factor(traits$Clutches),
	'LarvalColor' = as.factor(traits$LarvalColor),
	'AdultSize' = traits$AdultSize,
	'AdultColor' = as.factor(traits$AdultColor),
	stringsAsFactors=F)
	
test = match(trends$Species,traits2$Species)
trends$Species[which(is.na(test))]

trait.trends = merge(trends,traits2,by='Species',all.x=F,sort=F); str(trait.trends)
species = trait.trends$Species
diet = read.csv('./data/NABA_HostPlants.csv',as.is=T,check.names=F,header=T)
diet$Species = paste(diet$LEP_genus,diet$LEP_species,sep=' ')
diet.breadth.families = apply(array(species),1,function(x){length(unique(diet$HOST_family[which(diet$Species==x)]))})
diet.breadth.genera = apply(array(species),1,function(x){length(unique(diet$HOST_genus[which(diet$Species==x)]))})
trait.trends$Diet.breadth.families = diet.breadth.families
trait.trends$Diet.breadth.genera = diet.breadth.genera

data1 = merge(data0,trait.trends,by='Species',all.x=F,sort=F); str(data1)
data1 = data1[which(!is.na(data1$AdultSize)),] #9 species lack size data
write.table(data1,'./data/butterfly_traits_envars_trends_50km_wNAs.txt',sep='\t',row.names=F,quote=F)

# Prune dataframe
data2 = data1[,c(1:13,25:27,30,39,41)]
species = unique(data2$Species) #545 species remain
sp.grid.count = apply(array(species),1,function(x){length(unique(data2$grid_id[which(data2$Species==x)]))})
perc1 = length(unique(data2$grid_id))*0.01
species.keep = species[which(sp.grid.count>perc1)] #394 out of 554 species occupied at least 10% of grid cells
keep.pos = match(data2$Species,species.keep)
data3 = data2[which(!is.na(keep.pos)),]

abundance.sd = sd(data3$Abundance.trend)
length(which(abs(data3$Abundance.trend) > (abundance.sd*5))) #46 trends removed because > 5 standard deviations

data4 = data3[which(abs(data3$Abundance.trend) < (abundance.sd*5)),]
hist(data4$Abundance.trend)
summary(data4$Abundance.trend)

write.table(data4,'./data/butterfly_traits_envars_trends_50km_wNAs_trim.txt',sep='\t',row.names=F,quote=F)


###############################################################
# Explore some models

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

data1 = read.table('./data/butterfly_traits_envars_trends_50km_wNAs_trim.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Species = as.factor(data1$Species)
data1$grid_id = as.factor(data1$grid_id)

#Standardize the landscape data using z-score to get to same scale for AIC comparison
attach(data1,warn.conflicts=F)

Cropland.2015 <- as.vector(data1$Cropland.2015); mCropland.2015 <- mean(Cropland.2015, na.rm=TRUE); sdCropland.2015 <- sd(Cropland.2015, na.rm=TRUE)
Cropland.2015z <- as.vector( (Cropland.2015-mCropland.2015) / sdCropland.2015 )

Cropland.diff <- as.vector(data1$Cropland.diff); mCropland.diff <- mean(Cropland.diff, na.rm=TRUE);	sdCropland.diff <- sd(Cropland.diff, na.rm=TRUE)
Cropland.diffz <- as.vector( (Cropland.diff-mCropland.diff) / sdCropland.diff )

Built.2015 <- as.vector(data1$Built.2015); mBuilt.2015 <- mean(Built.2015, na.rm=TRUE);	sdBuilt.2015 <- sd(Built.2015, na.rm=TRUE)
Built.2015z <- as.vector( (Built.2015-mBuilt.2015) / sdBuilt.2015 )

Built.diff <- as.vector(data1$Built.diff); mBuilt.diff <- mean(Built.diff, na.rm=TRUE);	sdBuilt.diff <- sd(Built.diff, na.rm=TRUE)
Built.diffz <- as.vector( (Built.diff-mBuilt.diff) / sdBuilt.diff )

Precip.2018 <- as.vector(data1$Precip.2018); mPrecip.2018 <- mean(Precip.2018, na.rm=TRUE);	sdPrecip.2018 <- sd(Precip.2018, na.rm=TRUE)
Precip.2018z <- as.vector( (Precip.2018-mPrecip.2018) / sdPrecip.2018 )

Precip.diff <- as.vector(data1$Precip.diff); mPrecip.diff <- mean(Precip.diff, na.rm=TRUE);	sdPrecip.diff <- sd(Precip.diff, na.rm=TRUE)
Precip.diffz <- as.vector( (Precip.diff-mPrecip.diff) / sdPrecip.diff )

Temp.2018 <- as.vector(data1$Temp.2018); mTemp.2018 <- mean(Temp.2018, na.rm=TRUE);	sdTemp.2018 <- sd(Temp.2018, na.rm=TRUE)
Temp.2018z <- as.vector( (Temp.2018-mTemp.2018) / sdTemp.2018 )

Temp.diff <- as.vector(data1$Temp.diff); mTemp.diff <- mean(Temp.diff, na.rm=TRUE); sdTemp.diff <- sd(Temp.diff, na.rm=TRUE)
Temp.diffz <- as.vector( (Temp.diff-mTemp.diff) / sdTemp.diff )

Diet.breadth.families <- as.vector(data1$Diet.breadth.families); mDiet.breadth.families <- mean(Diet.breadth.families, na.rm=TRUE); sdDiet.breadth.families <- sd(Diet.breadth.families, na.rm=TRUE)
Diet.breadth.familiesz <- as.vector( (Diet.breadth.families-mDiet.breadth.families) / sdDiet.breadth.families )

AdultSize <- as.vector(data1$AdultSize); mAdultSize <- mean(AdultSize, na.rm=TRUE); sdAdultSize <- sd(AdultSize, na.rm=TRUE)
AdultSizez <- as.vector( (AdultSize-mAdultSize) / sdAdultSize )

data1$Cropland.2015z = Cropland.2015z
data1$Built.2015z = Built.2015z
data1$Cropland.diffz = Cropland.diffz
data1$Built.diffz = Built.diffz
data1$Precip.2018z = Precip.2018z
data1$Temp.2018z = Temp.2018z
data1$Precip.diffz = Precip.diffz
data1$Temp.diffz = Temp.diffz
data1$Diet.breadth.familiesz = Diet.breadth.familiesz
data1$AdultSizez = AdultSizez

write.table(data1,'./data/butterfly_traits_envars_trends_50km_wNAs_trim_forModelling.txt',sep='\t',row.names=F,quote=F)


##############################################
# Run models and find best correlation structure


mod1 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1)
summary(mod1)
AIC(mod1) #125801.5

soil.gau <- update(mod1, correlation = corGaus(1, form = ~ Lon + Lat), method = "ML")
summary(soil.gau)

soil.exp <- update(mod1, correlation = corExp(1, form = ~ Lon + Lat), method = "ML")
summary(soil.exp)

soil.spher <- update(mod1, correlation = corSpher(1, form = ~ Lon + Lat), method = "ML")
summary(soil.spher)

AIC(soil.gau,soil.exp,soil.spher)


#####################
# Hand-build models to find best model by parsimony

data1 = read.table("./data/butterfly_traits_envars_trends_50km_wNAs_trim_forModelling.txt",sep='\t',as.is=T,check.names=F,header=TRUE)

# Full +/- traits
mod0 = lme(Abundance.trend ~ 1, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod1 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod2 = lme(Abundance.trend ~ AdultSizez + Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod3 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod4 = lme(Abundance.trend ~ Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# No traits, leave one env out
mod5 = lme(Abundance.trend ~ Built.2015z*Built.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod6 = lme(Abundance.trend ~ Cropland.2015z*Cropland.diffz + Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod7 = lme(Abundance.trend ~ Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod8 = lme(Abundance.trend ~ Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz + Precip.2018z*Precip.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# Landscape only
mod9 = lme(Abundance.trend ~ Cropland.2015z*Cropland.diffz + Built.2015z*Built.diffz , random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod10 = lme(Abundance.trend ~ Cropland.diffz + Built.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# Climate only
mod11 = lme(Abundance.trend ~ Precip.2018z*Precip.diffz + Temp.2018z*Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod12 = lme(Abundance.trend ~ Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# All diffs +/- traits
mod13 = lme(Abundance.trend ~ Cropland.diffz + Built.diffz + Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod14 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.diffz + Built.diffz + Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod15 = lme(Abundance.trend ~ AdultSizez + Cropland.diffz + Built.diffz + Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod16 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.diffz + Built.diffz + Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# One at a time
mod17 = lme(Abundance.trend ~ AdultSizez, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod18 = lme(Abundance.trend ~ Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod19 = lme(Abundance.trend ~ Cropland.2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod20 = lme(Abundance.trend ~ Cropland.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod21 = lme(Abundance.trend ~ Built.2015z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod22 = lme(Abundance.trend ~ Built.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod23 = lme(Abundance.trend ~ Precip.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod24 = lme(Abundance.trend ~ Precip.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod25 = lme(Abundance.trend ~ Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod26 = lme(Abundance.trend ~ Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# Traits only
mod27 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# All contemporary +/- traits
mod28 = lme(Abundance.trend ~ Cropland.2015z + Built.2015z + Precip.2018z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod29 = lme(Abundance.trend ~ AdultSizez + Diet.breadth.familiesz + Cropland.2015z + Built.2015z + Precip.2018z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod30 = lme(Abundance.trend ~ AdultSizez + Cropland.2015z + Built.2015z + Precip.2018z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod31 = lme(Abundance.trend ~ Diet.breadth.familiesz + Cropland.2015z + Built.2015z + Precip.2018z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# No traits or diffs, leave one env out
mod32 = lme(Abundance.trend ~ Built.2015z + Precip.2018z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod33 = lme(Abundance.trend ~ Cropland.2015z + Precip.2018z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod34 = lme(Abundance.trend ~ Cropland.2015z + Built.2015z + Temp.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod35 = lme(Abundance.trend ~ Cropland.2015z + Built.2015z + Precip.2018z, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
# No traits or contemporary, leave one env out
mod36 = lme(Abundance.trend ~ Built.diffz + Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod37 = lme(Abundance.trend ~ Cropland.diffz + Precip.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod38 = lme(Abundance.trend ~ Cropland.diffz + Built.diffz + Temp.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
mod39 = lme(Abundance.trend ~ Cropland.diffz + Built.diffz + Precip.diffz, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))

AICctab(mod0.mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,
	mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39, weights=TRUE)

modav = model.avg(list(mod32,mod23,mod11))

confint(modav, level=0.90)
importance(modav)


##############################################################
# Total abundance (per grid cell)

library(nlme)
library(MuMIn)
library(bbmle)

data = read.table("./data/butterfly_total_abundance_envars_50km.txt",sep='\t',as.is=T,check.names=F,header=TRUE)
attach(data)
summary(data)

# Z-transform covariates
Cropland.2015 <- as.vector(data$Cropland.2015)
mCropland.2015 <- mean(Cropland.2015, na.rm=TRUE)
sdCropland.2015 <- sd(Cropland.2015, na.rm=TRUE)
Cropland.2015z <- as.vector( (Cropland.2015-mCropland.2015) / sdCropland.2015 )

Cropland.diff <- as.vector(data$Cropland.diff)
mCropland.diff <- mean(Cropland.diff, na.rm=TRUE)
sdCropland.diff <- sd(Cropland.diff, na.rm=TRUE)
Cropland.diffz <- as.vector( (Cropland.diff-mCropland.diff) / sdCropland.diff )

Built.2015 <- as.vector(data$Built.2015)
mBuilt.2015 <- mean(Built.2015, na.rm=TRUE)
sdBuilt.2015 <- sd(Built.2015, na.rm=TRUE)
Built.2015z <- as.vector( (Built.2015-mBuilt.2015) / sdBuilt.2015 )

Built.diff <- as.vector(data$Built.diff)
mBuilt.diff <- mean(Built.diff, na.rm=TRUE)
sdBuilt.diff <- sd(Built.diff, na.rm=TRUE)
Built.diffz <- as.vector( (Built.diff-mBuilt.diff) / sdBuilt.diff )

Precip.2018 <- as.vector(data$Precip.2018)
mPrecip.2018 <- mean(Precip.2018, na.rm=TRUE)
sdPrecip.2018 <- sd(Precip.2018, na.rm=TRUE)
Precip.2018z <- as.vector( (Precip.2018-mPrecip.2018) / sdPrecip.2018 )

Precip.diff <- as.vector(data$Precip.diff)
mPrecip.diff <- mean(Precip.diff, na.rm=TRUE)
sdPrecip.diff <- sd(Precip.diff, na.rm=TRUE)
Precip.diffz <- as.vector( (Precip.diff-mPrecip.diff) / sdPrecip.diff )

Temp.2018 <- as.vector(data$Temp.2018)
mTemp.2018 <- mean(Temp.2018, na.rm=TRUE)
sdTemp.2018 <- sd(Temp.2018, na.rm=TRUE)
Temp.2018z <- as.vector( (Temp.2018-mTemp.2018) / sdTemp.2018 )

Temp.diff <- as.vector(data$Temp.diff)
mTemp.diff <- mean(Temp.diff, na.rm=TRUE)
sdTemp.diff <- sd(Temp.diff, na.rm=TRUE)
Temp.diffz <- as.vector( (Temp.diff-mTemp.diff) / sdTemp.diff )


##make global model 
data.spatialCor.glsSpher = gls(Abundance.trend ~  Cropland.2015z + Cropland.diffz + Built.2015z + Built.diffz + Precip.2018z + Precip.diffz + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz + Precip.2018z:Precip.diffz + Temp.2018z:Precip.2018z,
	data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")

dd <- dredge(data.spatialCor.glsSpher) # takes a while
write.csv(dd,"AICcdredge_Total_abundance_grid.csv", row.names = FALSE)

subset(dd, delta < 4)

#models with delta.aicc < 4
modav = model.avg(get.models(dd, subset = delta < 4)) # get averaged coefficients

confint(modav, level = 0.90)

mod1 = gls(Abundance.trend ~  Temp.2018z + Temp.diffz , data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod2 = gls(Abundance.trend ~  Cropland.2015z + Temp.2018z + Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod3 = gls(Abundance.trend ~  Precip.2018z + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz + Temp.2018z:Precip.2018z,	data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod4 = gls(Abundance.trend ~  Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod5 = gls(Abundance.trend ~  Precip.2018z + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod6 = gls(Abundance.trend ~  Precip.2018z + Temp.2018z + Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod7 = gls(Abundance.trend ~  Cropland.2015z + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod8 = gls(Abundance.trend ~  Precip.2018z + Temp.2018z + Temp.diffz + Temp.2018z:Precip.2018z, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod9 = gls(Abundance.trend ~  Temp.2018z, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod10 = gls(Abundance.trend ~  Cropland.2015z + Precip.2018z + Temp.2018z + Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod11 = gls(Abundance.trend ~  Precip.diffz + Temp.2018z + Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod12 = gls(Abundance.trend ~  Precip.2018z + Precip.diffz + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz + Temp.2018z:Precip.2018z, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod13 = gls(Abundance.trend ~  Cropland.2015z + Precip.2018z + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz + Temp.2018z:Precip.2018z, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod14 = gls(Abundance.trend ~  Precip.2018z + Precip.diffz + Temp.2018z + Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")
mod15 = gls(Abundance.trend ~  Cropland.2015z + Precip.2018z + Temp.2018z + Temp.diffz + Temp.2018z:Temp.diffz, data=data,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)
summary(mod8)
summary(mod9)
summary(mod10)
summary(mod11)
summary(mod12)
summary(mod13)
summary(mod14)
summary(mod15)


###################################################
# GLS with rarefied richness and evenness (per site)

library(nlme)
library(MuMIn)
library(bbmle)

con = file("./data/butterfly_diversity_envars_persite.txt",'r')
dat = readLines(con)
close(con)
data1 = c()
for (i in 1:length(dat)){
	d1 = strsplit(dat[i],'\t')[[1]]
	data1 = data.frame(rbind(data1,d1),stringsAsFactors=F)
}
colnames(data1) = data1[1,]
data1 = data1[-1,]
data1$Abundance.trend = as.numeric(data1$Abundance.trend)
data1$Richness.rarefied.trend = as.numeric(data1$Richness.rarefied.trend)
data1$Evenness.evar.rarefied.trend = as.numeric(data1$Evenness.evar.rarefied.trend)
data1$Cropland.2015 = as.numeric(data1$Cropland.2015)
data1$Cropland.diff = as.numeric(data1$Cropland.diff)
data1$Built.2015 = as.numeric(data1$Built.2015)
data1$Built.diff = as.numeric(data1$Built.diff)
data1$Precip.2018 = as.numeric(data1$Precip.2018)
data1$Precip.diff = as.numeric(data1$Precip.diff)
data1$Temp.2018 = as.numeric(data1$Temp.2018)
data1$Temp.diff = as.numeric(data1$Temp.diff)
data1$Lon = as.numeric(data1$Lon)
data1$Lat = as.numeric(data1$Lat)
data1$Elevation = as.numeric(data1$Elevation)

data1 = data1[which(!is.na(data1$Richness.rarefied.trend)),]
attach(data1)
summary(data1)

# Z-transform covariates
Cropland.2015 <- as.vector(data1$Cropland.2015)
mCropland.2015 <- mean(Cropland.2015, na.rm=TRUE)
sdCropland.2015 <- sd(Cropland.2015, na.rm=TRUE)
Cropland.2015z <- as.vector( (Cropland.2015-mCropland.2015) / sdCropland.2015 )

Cropland.diff <- as.vector(data1$Cropland.diff)
mCropland.diff <- mean(Cropland.diff, na.rm=TRUE)
sdCropland.diff <- sd(Cropland.diff, na.rm=TRUE)
Cropland.diffz <- as.vector( (Cropland.diff-mCropland.diff) / sdCropland.diff )

Built.2015 <- as.vector(data1$Built.2015)
mBuilt.2015 <- mean(Built.2015, na.rm=TRUE)
sdBuilt.2015 <- sd(Built.2015, na.rm=TRUE)
Built.2015z <- as.vector( (Built.2015-mBuilt.2015) / sdBuilt.2015 )

Built.diff <- as.vector(data1$Built.diff)
mBuilt.diff <- mean(Built.diff, na.rm=TRUE)
sdBuilt.diff <- sd(Built.diff, na.rm=TRUE)
Built.diffz <- as.vector( (Built.diff-mBuilt.diff) / sdBuilt.diff )

Precip.2018 <- as.vector(data1$Precip.2018)
mPrecip.2018 <- mean(Precip.2018, na.rm=TRUE)
sdPrecip.2018 <- sd(Precip.2018, na.rm=TRUE)
Precip.2018z <- as.vector( (Precip.2018-mPrecip.2018) / sdPrecip.2018 )

Precip.diff <- as.vector(data1$Precip.diff)
mPrecip.diff <- mean(Precip.diff, na.rm=TRUE)
sdPrecip.diff <- sd(Precip.diff, na.rm=TRUE)
Precip.diffz <- as.vector( (Precip.diff-mPrecip.diff) / sdPrecip.diff )

Temp.2018 <- as.vector(data1$Temp.2018)
mTemp.2018 <- mean(Temp.2018, na.rm=TRUE)
sdTemp.2018 <- sd(Temp.2018, na.rm=TRUE)
Temp.2018z <- as.vector( (Temp.2018-mTemp.2018) / sdTemp.2018 )

Temp.diff <- as.vector(data1$Temp.diff)
mTemp.diff <- mean(Temp.diff, na.rm=TRUE)
sdTemp.diff <- sd(Temp.diff, na.rm=TRUE)
Temp.diffz <- as.vector( (Temp.diff-mTemp.diff) / sdTemp.diff )


#######
# Richness (per site)

data.spatialCor.glsSpher = gls(Richness.rarefied.trend ~ Cropland.2015z + Cropland.diffz + Built.2015z + Built.diffz + Precip.2018z + Precip.diffz + Temp.2018z + Temp.diffz + Temp.2018z:Precip.2018z,
	data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")

dd <- dredge(data.spatialCor.glsSpher) # takes a while
subset(dd, delta < 4)
# Intercept-only is the best model
write.csv(dd,"AICcdredge_richness_rarefied_persite.csv", row.names = FALSE)

######
# Evenness (per site)

data.spatialCor.glsSpher = gls(Evenness.evar.rarefied.trend ~  Cropland.2015z + Cropland.diffz + Built.2015z + Built.diffz + Precip.2018z + Precip.diffz + Temp.2018z + Temp.diffz + Temp.2018z:Precip.2018z,
	data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")

dd <- dredge(data.spatialCor.glsSpher) # takes a while
subset(dd, delta < 4)
write.csv(dd,"AICcdredge_evenness_rarefied_persite.csv", row.names = FALSE)

mod1 = gls(Evenness.evar.rarefied.trend ~ Temp.2018z,data=data1,correlation=corSpher(form = ~Lat + Lon,nugget=T),method="REML")

# Plot regression
source('./code/visreg.utilities.R')
Visreg = function(fit, xvar, by, breaks = 3, type="conditional", data=NULL, trans=I, 
	scale="linear", alpha=0.05, nn=101, cond=list(),
    jitter=FALSE, collapse=FALSE, plot = TRUE, ...){
    Data <- setupF(fit, xvar, parent.frame(), data)
    xvar <- attr(Data, "xvar")
    cond <- setupCond(cond, Data, by, breaks)
    yName <- makeYName(fit, scale, trans, type)
    v <- setupV(fit, Data, xvar, nn, cond, type, trans, alpha, by, yName)
	par(oma=c(0,0,0,0),mar=c(5,5,1,1))
	plot(v,bty='n',xaxt='n',yaxt='n',ylim=c(-0.5,0.4),cex.lab=2,pch=16,cex=2,lwd=3,main='',ylab='Evenness trend')
	axis(1,lwd=3,cex.axis=1.5)
	axis(2,lwd=3,cex.axis=1.5,at=seq(-0.5,0.4,0.15))
}
rwb <- colorRampPalette(colors = c("orange", "white", "blue"))(50)
Visreg(fit=mod1,xvar='Temp.2018z')



