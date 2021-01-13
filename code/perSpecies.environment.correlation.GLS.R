

# https://www.flutterbys.com.au/stats/tut/tut8.4a.html
# https://www.rdocumentation.org/packages/MuMIn/versions/1.3.6/topics/dredge

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
### Traits

trends.summary = read.table('./data/butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
trends = trends.summary[,c(1,15,19:23)]
colnames(trends)[which(colnames(trends)=='mdn_Tau')] = 'Median.abundance.trend'

traits = read.table('./data/NABA_traits.txt',sep='\t',as.is=T,check.names=F,header=T); str(traits)
traits$Habitat[which(traits$Habitat=='')] = NA
traits$NCGR[which(traits$NCGR=='')] = NA
traits$OverwinteringStage[which(traits$OverwinteringStage=='')] = NA
traits$Clutches[which(traits$Clutches=='')] = NA
traits$LarvalColor[which(traits$LarvalColor=='')] = NA
traits$Hair[which(traits$Hair=='')] = NA
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
	'LarvalHair' = as.factor(traits$Hair),
	'AdultSize' = traits$AdultSize,
	'AdultColor' = as.factor(traits$AdultColor),
	stringsAsFactors=F)
dim(traits2) #652 species
traits2 = traits2[which(!is.na(traits2$AdultSize) & !is.na(traits2$NCGR) & !is.na(traits2$LarvalHair) & !is.na(traits2$LarvalColor)),]
dim(traits2) #435 species with complete set of 6 traits

data1 = merge(trends,traits2,by='Species',all.x=F,sort=F); str(data1)
species = data1$Species
diet = read.csv('./data/NABA_HostPlants_UniqueGenera_wRefs_061820.updated.csv',as.is=T,check.names=F,header=T)
diet$Species = paste(diet$LEP_genus,diet$LEP_species,sep=' ')
diet.breadth.families = apply(array(species),1,function(x){length(unique(diet$HOST_family[which(diet$Species==x)]))})
data1$Diet.breadth.families = diet.breadth.families

# Ensure categorical variables are factors
data1$Species = as.factor(data1$Species)
data1$Family = as.factor(data1$Family)
data1$NCGR = as.factor(data1$NCGR)
data1$LarvalColor = as.factor(data1$LarvalColor)
data1$LarvalHair = as.factor(data1$LarvalHair)
data1$AdultColor = as.factor(data1$AdultColor)

# Z-transform numeric covariates
data1$Diet.breadth.familiesz = (data1$Diet.breadth.families - mean(data1$Diet.breadth.families)) / sd(data1$Diet.breadth.families)
data1$AdultSizez = (data1$AdultSize - mean(data1$AdultSize)) / sd(data1$AdultSize)

# Model median abundance trends as a function of butterfly traits
mod1 = lme(Median.abundance.trend ~ NCGR + LarvalColor + LarvalHair + AdultColor + AdultSizez + Diet.breadth.familiesz, random= ~ 1 | Family/Species, data=data1)
summary(mod1)
intervals(mod1,level=0.95,which='fixed')

families = sort(unique(as.character(data1$Family)))
fm.count = apply(array(families),1,function(x){length(which(as.character(data1$Family)==x))})
cbind(families,fm.count)


#####
### Environmental variables

data0 = read.table('./data/butterfly_traits_envars_trends_50km_NoMergeJunAug_m5_trim_6traits.txt',sep='\t',as.is=T,check.names=F,header=T); str(data1)
species = sort(unique(data0$Species))
sp.count = apply(array(species),1,function(x){length(which(data0$Species==x))})
species2 = species[which(sp.count>=10)]; length(species); length(species2) #178 species occupy at least 10 grid cells

out = data.frame('Species'=NA,'Variable'=NA,'Value'=-999,'Std.Error'=-999,'DF'=-999,'t-value'=-999,'p-value'=-999)
for (s in 1:length(species2)){
	print(noquote(species2[s]))
	data1 = data0[which(data0$Species==species2[s]),]
	data1$grid_id = as.factor(data1$grid_id)
	# Z-transform numeric covariates
	data1$Cropland.2005_2015z = (data1$Cropland.2005_2015 - mean(data1$Cropland.2005_2015)) / sd(data1$Cropland.2005_2015)
	data1$Built.2005_2015z = (data1$Built.2005_2015 - mean(data1$Built.2005_2015)) / sd(data1$Built.2005_2015)
	data1$Cropland.trendz = (data1$Cropland.trend - mean(data1$Cropland.trend)) / sd(data1$Cropland.trend)
	data1$Precip.1993_2018z = (data1$Precip.1993_2018 - mean(data1$Precip.1993_2018)) / sd(data1$Precip.1993_2018)
	data1$Temp.1993_2018z = (data1$Temp.1993_2018 - mean(data1$Temp.1993_2018)) / sd(data1$Temp.1993_2018)
	data1$Precip.trendz = (data1$Precip.trend - mean(data1$Precip.trend)) / sd(data1$Precip.trend)
	data1$Temp.trendz = (data1$Temp.trend - mean(data1$Temp.trend)) / sd(data1$Temp.trend)
	# Model abundance trends as a function of environmental variables
	mod1 = lme(Abundance.trend ~ Cropland.2005_2015z + Cropland.trendz + Built.2005_2015z + Precip.1993_2018z + Precip.trendz + Temp.1993_2018z + Temp.trendz, data=data1, random= ~ 1 | grid_id, method='ML')
	add1 = data.frame(summary(mod1)$tTable)
	add1$Variable = rownames(add1)
	add1$Species = species2[s]
	add1 = add1[,c(6,7,1:5)]
	out = data.frame(rbind(out,add1),stringsAsFactors=F)
}
out = out[-1,]

variables = unique(out$Variable)
out2 = data.frame('Variable'=variables,'Sig.up'=-999,'Sig.down'=-999)
for (v in 1:length(variables)){
	dat = out$Value[which(out$p.value<0.05 & out$Variable==variables[v])]
	out2$Sig.up[v] = length(which(dat>0))
	out2$Sig.down[v] = length(which(dat<0))
}

png('./plots/environment effects species vote count.png',res=300,width=480*4,height=480*4)
par(oma=c(0,0,0,0),mar=c(10,5,1,4))
boxplot(list(
	out$Value[which(out$p.value<0.05 & out$Variable==variables[2])],
	out$Value[which(out$p.value<0.05 & out$Variable==variables[3])],
	out$Value[which(out$p.value<0.05 & out$Variable==variables[4])],
	out$Value[which(out$p.value<0.05 & out$Variable==variables[5])],
	out$Value[which(out$p.value<0.05 & out$Variable==variables[6])],
	out$Value[which(out$p.value<0.05 & out$Variable==variables[7])],
	out$Value[which(out$p.value<0.05 & out$Variable==variables[8])]),
	xaxt='n',yaxt='n',frame=F,xlab='',ylab='Covariate effect',cex.lab=2,ylim=c(-5,5),
	boxcol='grey80',whisklty=1,whisklwd=2,staplecol='white',outpch=16,outcol='grey40',outcex=0.5
)
axis(1,lwd=3,labels=F);axis(2,lwd=3,cex.axis=1.5)
title(xlab='Variable',line=8,cex.lab=2)
text(x=1:7,y=-6,labels=c('Prop. cropland','Cropland trend','Prop. built','Precipitation','Precip. trend','Temperature','Temp. trend'),cex=1.5,srt=-45,xpd=T,adj=0)
abline(h=0,lty=2,col='grey50')
dev.off()

