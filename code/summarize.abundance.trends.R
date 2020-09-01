

##########################################################################
# Median trends across cells per species

library(sf)
library(rgdal)
library(ggplot2)
library(rgeos)
library(viridis)
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

############################
# JunAug set
butterfly_counts = read.table('./data/butterfly_data_gridded_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km_NoMergeJunAug_m5.shp',layer='butterfly_sites_50km_NoMergeJunAug_m5',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
butterfly_grid = readOGR(dsn='./shapefiles/butterfly_grid_50km_NoMergeJunAug_m5.shp',layer='butterfly_grid_50km_NoMergeJunAug_m5',verbose=F,stringsAsFactors=F) #shapefile with grid
butterfly_grid2 = butterfly_grid[match(unique(butterfly_circles@data$grid_id),butterfly_grid@data$grid_id),]

# summarize sites per grid
sites = unique(butterfly_counts$Site)
counts2 = butterfly_counts[match(sites,butterfly_counts$Site),]
gids = unique(counts2$grid_id)
gcount = apply(array(gids),1,function(x){length(which(counts2$grid_id==x))})
mean(gcount) #1.203
sd(gcount)/sqrt(length(gcount)) #0.023

files = list.files('C:/Users/mcros/Desktop/Postdoc UGA/NABA_butterflies/inla_model_output/50km_wNAs_noMerge/',full.names=T)
files = files[grep('JunAug',files)]

out = c()
for (f in 1:length(files)){
	print(noquote(f))
	file1 = read.table(files[f],sep='\t',as.is=T,check.names=F,header=T)
	file1$Species = strsplit(strsplit(files[f],'/')[[1]][9],'_')[[1]][1]
	out = data.frame(rbind(out,file1),stringsAsFactors=F)
}

nrow(out)
length(which(!is.na(out$tau_sig))) #1,908
length(which(!is.na(out$tau_sig))) / nrow(out) #1.0%
sigtrends = out$tau_sig[which(!is.na(out$tau_sig))]
length(which(sigtrends>0)) #517
length(which(sigtrends>0)) / length(sigtrends) #27.1%
length(which(sigtrends<0)) #1,391
length(which(sigtrends<0)) / length(sigtrends) #72.9%

# Do trends vary by ecoregion?
data1 = read.table('./data/butterfly_traits_envars_trends_50km_NoMergeJul_m5_trim_6traits.txt',sep='\t',as.is=T,check.names=F,header=T); str(data1)
data1$Species = as.factor(data1$Species)
data1$Family = as.factor(data1$Family)
data1$EcoI = NA
for (i in 1:nrow(data1)){
	data1$EcoI[i] = out$EcoI[which(out$grid_id==data1$grid_id[i] & out$Species==data1$Species[i])]
}
data2 = data1[which(data1$EcoI!=0),]
data2$EcoI = as.factor(data2$EcoI)
library(nlme)
mod1 = lme(Abundance.trend ~ EcoI, random= ~ 1 | Family/Species, data=data1, correlation=corExp(1, form=~Lon+Lat))
summary(mod1)
par(oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(data2$EcoI,data2$Abundance.trend,xlab='Ecoregion',ylab='Abundance trend (%/yr)',cex.lab=2,ylim=c(-20,30),yaxt='n',frame=F,boxcol='grey70',medcol='black',whiskcol='black',whisklty=1,col='grey70',outpch=16,outcol='grey70',staplecol='white',whisklwd=2,cex.axis=1.5)
abline(h=0,lwd=2,lty=2,col='pink')
axis(1,lwd=3,labels=F,at=1:9)
axis(2,lwd=3,at=seq(-20,30,10),cex.axis=1.5)

# Manuscript SI - mean abundance trends barplot
species = unique(out$Species)
means = apply(array(species),1,function(x){mean(out$tau[which(out$Species==x)],na.rm=T)})
t1 = t.test(as.numeric(means))

# Barplot
mean1 = t1[['estimate']]
uppers = t1[['conf.int']][2]
lowers = t1[['conf.int']][1]
png('./plots/t-test species abundance barplot JunAug.png',width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mean1,ylim=c(-4,2),names.arg='',ylab='Avg. abundance trend',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='Species',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/15)
dev.off()

# Summarize trends per species
species = unique(out$Species)
sout = data.frame('Species'=NA,'N.cells'=-999,'N.increasing'=-999,'N.decreasing'=-999)
for (s in 1:length(species)){
	sdat = out[which(out$Species==species[s]),]
	sout[s,1] = species[s]
	sout[s,2] = nrow(sdat)
	sout[s,3] = length(which(sdat$tau>1))
	sout[s,4] = length(which(sdat$tau<(-1)))
}
sout$Perc.increasing = 100 * sout$N.increasing / sout$N.cells
sout$Perc.decreasing = 100 * sout$N.decreasing / sout$N.cells
sout$Perc.nochange = 100 - (sout$Perc.decreasing + sout$Perc.increasing)
length(which(sout$Perc.increasing>50)) #108 species increasing in over 50% of range
length(which(sout$Perc.increasing>50)) / nrow(sout) #23.7%
length(which(sout$Perc.decreasing>50)) #210 species decreasing in over 50% of range
length(which(sout$Perc.decreasing>50)) / nrow(sout) #46.1%
length(which(sout$Perc.nochange>50)) #39 species exhibited no change in over 50% of range
length(which(sout$Perc.nochange>50)) / nrow(sout) #8.6%

length(which(sout$Perc.increasing>90)) #61 species increasing in over 90% of range
length(which(sout$Perc.increasing>90)) / nrow(sout) #13.4%
length(which(sout$Perc.decreasing>90)) #122 species decreasing in over 90% of range
length(which(sout$Perc.decreasing>90)) / nrow(sout) #26.8%
length(which(sout$Perc.nochange>90)) #4 species exhibited no change in over 90% of range
length(which(sout$Perc.nochange>90)) / nrow(sout) #0.9%


# Summarize species trends per grid cell
gids = unique(out$grid_id)
out2 = data.frame('grid_id'=gids,'decreasing'=-999,'increasing'=-999,'nochange'=-999,'n.species'=-999,'mdn_Tau'=-999,'avg_Tau_inv'=-999)
for (g in 1:length(gids)){
	outg = out[which(out$grid_id==gids[g]),]
	out2[g,2] = length(which(outg$tau<(-1)))
	out2[g,3] = length(which(outg$tau>1))
	out2[g,4] = length(which(outg$tau>(-1) & outg$tau<1))
	out2[g,5] = length(unique(outg$Species))
	out2[g,6] = median(outg$tau,na.rm=T)
	outg2 = outg[which(outg$Species=='Erebia theano' | outg$Species=='Pieris rapae' | outg$Species=='Thymelicus lineola'),]
	if (nrow(outg2)==0){
		out2[g,7] = NA
	} else {
		out2[g,7] = mean(outg2$tau,na.rm=T)
	}
}
out3 = merge(butterfly_grid2, out2, by="grid_id", all=F)
out3$perc.decreasing = out3$decreasing / length(gids)
out3$perc.increasing = out3$increasing / length(gids)
out3$perc.nochange = out3$nochange / length(gids)
out3$ratio = out3$decreasing / out3$increasing
writeOGR(out3,dsn='./shapefiles/butterfly_species_abundance_trend_summary_JunAug.shp',layer='butterfly_species_abundance_trend_summary_JunAug.shp',driver='ESRI Shapefile',overwrite=T)
out3_sf <- as(out3, "sf")

ecos = unique(out$EcoI)
ecosum = data.frame('EcoregionI'=ecos,'Decreasing10'=-999,'Increasing10'=-999,'Decreasing50'=-999,'Increasing50'=-999,'Decreasing90'=-999,'Increasing90'=-999,'N.cells'=-999,'N.species'=-999,'Avg.ncells.per.species'=-999)
secosum = c()
for (e in 1:length(ecos)){
	print(noquote(ecos[e]))
	edat = out[which(out$EcoI==ecos[e]),]
	species = unique(edat$Species)
	spout = data.frame('Species'=species,'ncells'=-999,'up'=-999,'down'=-999,'eco1'=NA)
	for (s in 1:length(species)){
		sedat = edat[which(edat$Species==species[s]),]
		spout[s,2] = nrow(sedat)
		spout[s,3] = length(which(sedat$tau>1))
		spout[s,4] = length(which(sedat$tau<(-1)))
		spout[s,5] = ecos[e]
	}
	secosum = data.frame(rbind(secosum,spout),stringsAsFactors=F)
	down.perc = spout$down / spout$ncells
	up.perc = spout$up / spout$ncells
	ecosum[e,2] = length(which(down.perc>0.1))
	ecosum[e,3] = length(which(up.perc>0.1))
	ecosum[e,4] = length(which(down.perc>0.5))
	ecosum[e,5] = length(which(up.perc>0.5))
	ecosum[e,6] = length(which(down.perc>0.9))
	ecosum[e,7] = length(which(up.perc>0.9))
	ecosum[e,8] = length(unique(edat$grid_id))
	ecosum[e,9] = length(unique(edat$Species))
	ecosum[e,10] = mean(spout$ncells)
}
# Add column for all ecoregions combined (full dataset)
e = e + 1
species = unique(out$Species)
	spout = data.frame('Species'=species,'ncells'=-999,'up'=-999,'down'=-999,'eco1'=NA)
for (s in 1:length(species)){
	sedat = out[which(out$Species==species[s]),]
	spout[s,2] = nrow(sedat)
	spout[s,3] = length(which(sedat$tau>1))
	spout[s,4] = length(which(sedat$tau<(-1)))
	spout[s,5] = 'All sites'
}
secosum = data.frame(rbind(secosum,spout),stringsAsFactors=F)

down.perc = spout$down / spout$ncells
up.perc = spout$up / spout$ncells
ecosum[e,1] = 'All sites'
ecosum[e,2] = length(which(down.perc>0.1))
ecosum[e,3] = length(which(up.perc>0.1))
ecosum[e,4] = length(which(down.perc>0.5))
ecosum[e,5] = length(which(up.perc>0.5))
ecosum[e,6] = length(which(down.perc>0.9))
ecosum[e,7] = length(which(up.perc>0.9))
ecosum[e,8] = length(unique(out$grid_id))
ecosum[e,9] = length(unique(out$Species))
ecosum[e,10] = mean(spout$ncells)
write.csv(ecosum,'ecoregion_abundance_summary_JunAug.csv',quote=F,row.names=F)

# Manuscript figure 2
png('./plots/fig2b_JunAug.png',res=300,width=480*4,height=480*6)
par(mfrow=c(3,2),oma=c(0,0,0,0),mar=c(5,5,2,1))
# Panel A
dat50 = read.table('butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
dat50$mdn_Tau[which(dat50$mdn_Tau>30)] = 35 #4 trends
dat50$mdn_Tau[which(dat50$mdn_Tau<(-30))] = -35 #20 trends
h = hist(dat50$mdn_Tau,xaxt='n',yaxt='n',xlab='Median abundance trend (%/yr)',ylim=c(0,240),ylab='No. species',cex.lab=2,main='',col='grey50',breaks=seq(-35,35,5))
axis(1,lwd=3,cex.axis=1.5,at=seq(-35,35,10),labels=c('-30<',seq(-25,25,10),'>30'))
axis(2,lwd=3,cex.axis=1.5,seq(0,240,40))
abline(v=0,lwd=3,lty=1,col='pink')
# Panel B
secosum$perc.up = secosum$up / secosum$ncells
secosum$perc.down = secosum$down / secosum$ncells
dat = secosum[which(secosum$eco1=='All sites'),]
h3 = hist(dat$perc.down*100,breaks=seq(0,100,5),xaxt='n',yaxt='n',xlab='% Grid cells decreasing',ylim=c(0,120),ylab='No. species',cex.lab=2,main='',col='grey50')
axis(1,lwd=3,cex.axis=1.5)
axis(2,lwd=3,cex.axis=1.5,at=seq(0,150,30))
# Panel C
trend.summary = read.table('./diversity_trends_JunAug.txt',sep='\t',as.is=T,check.names=F,header=T)
hist(trend.summary$Richness.rarefied.trend,xlab='Rarefied richness trend (sd/yr)',ylim=c(0,350),breaks=seq(-0.6,0.6,0.05),ylab='No. sites',cex.lab=2,col='grey50',yaxt='n',xaxt='n',main='')
axis(1,lwd=3,cex.axis=1.5,at=round(seq(-0.6,0.6,0.1),2))
axis(2,lwd=3,cex.axis=1.5)
abline(v=0,lwd=3,lty=1,col='pink')
# Panel D
hist(trend.summary$Evenness.evar.rarefied.trend,xlab='Rarefied evenness trend (sd/yr)',ylim=c(0,150),breaks=seq(-0.6,0.5,0.05),ylab='No. sites',cex.lab=2,col='grey50',yaxt='n',xaxt='n',main='')
axis(1,lwd=3,cex.axis=1.5,at=round(seq(-0.6,0.5,0.1),1))
axis(2,lwd=3,cex.axis=1.5)
abline(v=0,lwd=3,lty=1,col='pink')
# Panel E
res_sf = read.table('./inla_model_output/50km_wNAs_noMerge/Allspecies_trends_1993-2018_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
hist(res_sf$tau,main='',ylim=c(0,100),breaks=seq(-14,8,1),ylab='No. grid cells',xlab='Median total abundance trend (%/yr)',col='grey50',cex.lab=2,xaxt='n',yaxt='n')
axis(1,lwd=3,cex.axis=1.5,seq(-14,8,2))
axis(2,lwd=3,cex.axis=1.5)
abline(v=0,lwd=3,col='pink')
# Panel F
res_sf = read.table('./inla_model_output/50km_wNAs_noMerge/Allspecies_trends_1993-2018_50km_noMergeJunAug_m5_noTP.txt',sep='\t',as.is=T,check.names=F,header=T)
hist(res_sf$tau,main='',ylim=c(0,100),breaks=seq(-10,8,1),ylab='No. grid cells',xlab='Median total abundance trend (%/yr)',col='grey50',cex.lab=2,xaxt='n',yaxt='n')
axis(1,lwd=3,cex.axis=1.5,seq(-10,6,2))
axis(2,lwd=3,cex.axis=1.5)
abline(v=0,lwd=3,col='pink')
dev.off()



secosum$up2down = secosum$up / secosum$down
hist(secosum$up2down)
secosum$down2up = secosum$down / secosum$up
hist(secosum$down2up)
length(which(secosum$up2down<1)) / nrow(secosum)
length(which(secosum$up2down>1)) / nrow(secosum)

# Histograms per ecoregion
png('./plots/increasing or decreasing abundance ecoregion histograms.png',height=480*2,width=480*2)
par(mfrow=c(5,4),oma=c(0,0,0,0),mar=c(5,5,3,1))
ymaxs = c(120,120,80,120,50,50,60,120,40,50)
for (e in 1:length(ecos)){
	tdat = secosum[which(secosum$eco1==ecos[e]),]
	hist(tdat$perc.up*100,xlab='% Increasing',ylab='No. species',ylim=c(0,ymaxs[e]),cex.lab=2,col='cornflowerblue',yaxt='n',xaxt='n',main=ecos[e],cex.main=1.5)
	axis(1,lwd=3,cex.axis=1.5)
	axis(2,lwd=3,cex.axis=1.5)
	hist(tdat$perc.down*100,xlab='% Decreasing',ylab='No. species',ylim=c(0,ymaxs[e]),cex.lab=2,col='orange',yaxt='n',xaxt='n',main=ecos[e],cex.main=1.5)
	axis(1,lwd=3,cex.axis=1.5)
	axis(2,lwd=3,cex.axis=1.5)
}
dev.off()




