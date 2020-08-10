
setwd('<your/working/directory>')


##########################################################################
# Median trends across sites per species

library(sf)
library(rgdal)
library(ggplot2)
library(rgeos)
library(viridis)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km.shp',layer='butterfly_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
butterfly_grid = spTransform(readOGR(dsn='./shapefiles/butterfly_grid_50km.shp',layer='butterfly_grid_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
butterfly_grid2 = butterfly_grid[match(unique(butterfly_circles@data$grid_id),butterfly_grid@data$grid_id),]
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

eco = unique(eco1@data[,1])
epos = apply(array(eco),1,function(x){which(eco1@data[,1]==x)[1]})
cbind(eco,eco1@data[epos,2])

files = list.files('./inla_model_output/50km_wNAs',full.names=T)
files = files[-grep('Allspecies',files)]

out = c()
for (f in 1:length(files)){
	print(noquote(f))
	file1 = read.table(files[f],sep='\t',as.is=T,check.names=F,header=T)
	file1$Species = strsplit(strsplit(files[f],'/')[[1]][9],'_')[[1]][1]
	out = data.frame(rbind(out,file1),stringsAsFactors=F)
}

nrow(out)
length(which(!is.na(out$tau_sig))) #2,177
length(which(!is.na(out$tau_sig))) / nrow(out) #7.7%
sigtrends = out$tau_sig[which(!is.na(out$tau_sig))]
length(which(sigtrends>0)) #649
length(which(sigtrends>0)) / length(sigtrends) #30.0%
length(which(sigtrends<0)) #1,528
length(which(sigtrends<0)) / length(sigtrends) #70.2%

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
length(which(sout$Perc.increasing>50)) #153 species increasing in over 50% of range
length(which(sout$Perc.increasing>50)) / nrow(sout) #27.7%
length(which(sout$Perc.decreasing>50)) #221 species decreasing in over 50% of range
length(which(sout$Perc.decreasing>50)) / nrow(sout) #40.0%
length(which(sout$Perc.nochange>50)) #89 species exhibited no change in over 50% of range
length(which(sout$Perc.nochange>50)) / nrow(sout) #16.1%

length(which(sout$Perc.increasing>90)) #64 species increasing in over 90% of range
length(which(sout$Perc.increasing>90)) / nrow(sout) #11.6%
length(which(sout$Perc.decreasing>90)) #145 species decreasing in over 90% of range
length(which(sout$Perc.decreasing>90)) / nrow(sout) #26.2%
length(which(sout$Perc.nochange>90)) #22 species exhibited no change in over 90% of range
length(which(sout$Perc.nochange>90)) / nrow(sout) #4.0%


# Summarize species trends per grid cell
gids = unique(out$grid_id)
out2 = data.frame('grid_id'=gids,'decreasing'=-999,'increasing'=-999,'nochange'=-999,'n.species'=-999)
for (g in 1:length(gids)){
	outg = out[which(out$grid_id==gids[g]),]
	out2[g,2] = length(which(outg$tau<(-1)))
	out2[g,3] = length(which(outg$tau>1))
	out2[g,4] = length(which(outg$tau>(-1) & outg$tau<1))
	out2[g,5] = length(unique(outg$Species))
}
out3 = merge(butterfly_grid2, out2, by="grid_id", all=F)
out3$perc.decreasing = out3$decreasing / length(gids)
out3$perc.increasing = out3$increasing / length(gids)
out3$perc.nochange = out3$nochange / length(gids)
out3$ratio = out3$decreasing / out3$increasing
out3_sf <- as(out3, "sf")

ecos = unique(out$EcorgnI)
ecosum = data.frame('EcoregionI'=ecos,'Decreasing10'=-999,'Increasing10'=-999,'Decreasing50'=-999,'Increasing50'=-999,'Decreasing90'=-999,'Increasing90'=-999,'N.cells'=-999,'N.species'=-999,'Avg.ncells.per.species'=-999)
secosum = c()
for (e in 1:length(ecos)){
	print(noquote(ecos[e]))
	edat = out[which(out$EcorgnI==ecos[e]),]
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
write.csv(ecosum,'ecoregion_abundance_summary.csv',quote=F,row.names=F)

secosum$perc.up = secosum$up / secosum$ncells
secosum$perc.down = secosum$down / secosum$ncells
dat = secosum[which(secosum$eco1=='All sites'),]
dat50 = read.table('butterfly_trends_summary_1993-2018_50km_wNAs.txt',sep='\t',as.is=T,check.names=F,header=T)
dat50$mdn_Tau[which(dat50$mdn_Tau>30)] = 35
dat50$mdn_Tau[which(dat50$mdn_Tau<(-30))] = -35
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(5,5,2,1))
h4 = hist(out_$tau,xaxt='n',yaxt='n',xlab='Cell-specific abundance trend',ylab='No. species*sites',cex.lab=2,main='',col='grey50',breaks=seq(-24,24,1))
axis(1,lwd=3,cex.axis=1.5,at=seq(-24,24,4),labels=c('-20<',seq(-20,20,4),'>20'))
axis(2,lwd=3,cex.axis=1.5)
abline(v=0,lwd=3,lty=1,col='pink')
h = hist(dat50$mdn_Tau,xaxt='n',yaxt='n',xlab='Median abundance trend',ylim=c(0,240),ylab='No. species',cex.lab=2,main='',col='grey50',breaks=seq(-35,35,5))
axis(1,lwd=3,cex.axis=1.5,at=seq(-35,35,10),labels=c('-30<',seq(-25,25,10),'>30'))
axis(2,lwd=3,cex.axis=1.5,seq(0,240,40))
abline(v=0,lwd=3,lty=1,col='pink')
check = dat50[which(dat50$mdn_Tau<(-80) | dat50$mdn_Tau>80),]
out_ = out
out_$tau[which(out_$tau>20)] = 24
out_$tau[which(out_$tau<(-20))] = -24
h2 = hist(dat$perc.up*100,xaxt='n',yaxt='n',xlab='% Grid cells increasing',ylim=c(0,270),ylab='No. species',cex.lab=2,main='',col='grey50')
axis(1,lwd=3,cex.axis=1.5)
axis(2,lwd=3,cex.axis=1.5,at=seq(0,270,30))
h3 = hist(dat$perc.down*100,xaxt='n',yaxt='n',xlab='% Grid cells decreasing',ylim=c(0,150),ylab='No. species',cex.lab=2,main='',col='grey50')
axis(1,lwd=3,cex.axis=1.5)
axis(2,lwd=3,cex.axis=1.5,at=seq(0,150,30))

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

