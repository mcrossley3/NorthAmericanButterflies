
#####################################
# Rank abundance curve metrics

library(vegan)
library(sads)
library(viridis)
library(PerformanceAnalytics)

# Define function that generates rank abundance curves
calc.ranks = function(sp1,ab1){
	#sp1 = a vector of species names (species records do not have to be unique)
	#ab1 = a vector of corresponding abundances of species in sp1
	pos1 = which(ab1>0)
	if (length(pos1)<8){
		return(NA)
	} else {
		sp2 = sp1[pos1]
		ab2 = ab1[pos1]
		order1 = order(ab2,decreasing=T)
		sp3 = sp2[order1]
		ab3 = ab2[order1]
		mod = rad.preempt(ab3)
		decay.rate = coef(mod)[1]
		dominance = ab3[1] / sum(ab3) #counts of most abundant species divided by total counts across all species
		sf = fitsad(ab3, sad="ls")
		fisher.alpha = coef(sf)[2]
		out = c(decay.rate,dominance,fisher.alpha)
		names(out) = c('Decay.rate','Dominance','Fishers.alpha')
		return(out)
	}
}

# Define function that calculates evenness according to Smith & Wilson 1996 "evar") See http://esapubs.org/archive/ecol/E093/190/appendix-B.htm
get.evar = function(x1){
	x1 = x1[x1>0]
	S = length(x1)
	y1 = sum(log(x1)/S)
	y2 = sum(apply(array(x1),1,function(a){((log(a)-y1)^2)/S}))
	evar = 1 - ((2/pi) * atan(y2))
	return(evar)
}

# Define function that calculates Pielou's evenness
get.pei = function(px){
	px = px[px>0]
	pei = diversity(px) / log(length(px))
	return(pei)
}

# Define function that rarefies evenness
rarefy2 = function (xx, sample1, se=FALSE, MARGIN=1, divfun){
    xx <- as.matrix(xx)
    minsample <- min(apply(xx, MARGIN, sum))
    rarefun <- function(y, sample2) {
        y <- y[y > 0]
        J <- sum(y)
        ldiv <- lchoose(J, sample2)
        p1 <- ifelse(J - y < sample2, 0, exp(lchoose(J - y, sample2) - ldiv)) #p1 is larger for rarer individuals
#       out <- sum(1 - p1) #original output of rarefy function
		p2 = 1-p1
		evenness.rare = divfun(p2)
		return(evenness.rare)
	}
    E.rare <- apply(xx, MARGIN, rarefun, sample2 = sample1)
    attr(E.rare, "Subsample") <- sample1
    return(E.rare)
}

# Define function that converts species abundance to species occurrence probabilities
get.probs = function (xx, sample1, se=FALSE, MARGIN=1, divfun){
    xx <- as.matrix(xx)
    minsample <- min(apply(x, MARGIN, sum))
    rarefun <- function(y, sample2) {
        y <- y[y > 0]
        J <- sum(y)
        ldiv <- lchoose(J, sample2)
        p1 <- ifelse(J - y < sample2, 0, exp(lchoose(J - y, sample2) - ldiv)) #p1 is larger for rarer individuals
#       out <- sum(1 - p1) #original output of rarefy function
		p2 = 1-p1
		evenness.rare = divfun(p2)
		return(p2)
	}
    E.rare <- apply(xx, MARGIN, rarefun, sample2 = sample1)
    attr(E.rare, "Subsample") <- sample1
    return(E.rare)
}


##################################
# JunAug set

butterfly_counts = read.table('./data/butterfly_data_NoMergeJunAug_w0s_m5.txt',sep='\t',as.is=T,check.names=F,header=T)

# Rarefy richness and evenness
sites = unique(butterfly_counts$Site)
species = sort(unique(butterfly_counts$Species))
years = 1993:2018

# Create list of community data matrices for each year
comm.mat.list = list()
keep.list = list()
for (y in 1:length(years)){
	print(years[y])
	mat1 = matrix(NA,nrow=length(sites),ncol=length(species))
	for (i in 1:length(sites)){
		pos1 = which(butterfly_counts$Year==years[y] & butterfly_counts$Site==sites[i])
		butterfly_counts2 = butterfly_counts[pos1,]
		mat1[i,] = butterfly_counts2$N.butterflies[match(species,butterfly_counts2$Species)]
	}
	na.count = apply(mat1,1,function(x){length(which(is.na(x)))})
	keep = which(na.count<length(species)) #position of sites that had any species counts in year y
	comm.mat.list[[y]] = mat1[keep,]
	keep.list[[y]] = keep
}

richness.rarefied = matrix(NA,nrow=length(sites),ncol=length(years))
for (y in 1:length(years)){
	print(years[y])
	mat1 = comm.mat.list[[y]]
	mat1[which(is.na(mat1))] = 0
	raremax <- min(rowSums(mat1))
	Srare <- rarefy(mat1, raremax)
#	sp.count = apply(mat1,1,function(x){length(which(x>0))})
#	plot(sp.count,Srare)
	keep = keep.list[[y]]
	richness.rarefied[keep,y] = Srare
}

library(scales)
evenness.pei.rarefied = matrix(NA,nrow=length(sites),ncol=length(years))
evenness.evar.rarefied = matrix(NA,nrow=length(sites),ncol=length(years))
for (y in 1:length(years)){
	print(years[y])
	mat1 = comm.mat.list[[y]]
	mat1[which(is.na(mat1))] = 0
	raremax <- min(rowSums(mat1))
	Erare.pei <- rarefy2(xx=mat1, sample1=raremax, divfun=get.pei)
	Erare.evar <- rarefy2(xx=mat1, sample1=raremax, divfun=get.evar)
	keep = keep.list[[y]]
	evenness.pei.rarefied[keep,y] = Erare.pei
	evenness.evar.rarefied[keep,y] = Erare.evar
}
plot(c(evenness.pei.rarefied),c(evenness.evar.rarefied),ylim=c(0,1),bty='n',xlab="Pielou's evenness rarefied",ylab="Evar rarefied",cex.lab=1.5,cex.axis=1.2,pch=16,col=alpha('black',0.5))


# Explore difference between raw and rarefied abundance
#dout = data.frame()
#for (y in 1:length(years)){ #takes a LONG time
#	print(years[y])
#	keep = keep.list[[y]]
#	sites1 = sites[keep]
#	mat1 = comm.mat.list[[y]]
#	mat1[which(is.na(mat1))] = 0
#	raremax <- min(rowSums(mat1))
#	for (i in 1:nrow(mat1)){ #for each site
#		y1 = mat1[i,][mat1[i,]>1]
#		J <- sum(y1)
#        ldiv <- lchoose(J, raremax)
#        p1 <- ifelse(J - y1 < raremax, 0, exp(lchoose(J - y1, raremax) - ldiv)) #p1 is larger for rarer individuals
#		p2 = 1-p1
#		p2sum = sum(p2)
#		p3 = mat1[i,]
#		p3[mat1[i,]>1] = p2
#		add1 = data.frame('Site'=rep(sites1[i],ncol(mat1)),'Year'=rep(years[y],ncol(mat1)),'Species'=species,'Min.sample'=rep(raremax,ncol(mat1)),'Raw.abundance'=mat1[i,],'Prob.abundance'=p3,'Sum.raw.abundance'=rep(J,ncol(mat1)),'Sum.prob.abundance'=rep(p2sum,ncol(mat1)),stringsAsFactors=F)
#		dout = data.frame(rbind(dout,add1),stringsAsFactors=F)
#	}
#}
#write.table(dout,'butterfly_abundance_probabilities.txt',sep='\t',row.names=F)
#plot(dout$Min.sample,(dout$Prob.abundance/dout$Sum.prob.abundance)-(dout$Raw.abundance/dout$Sum.raw.abundance),xlab='Minimum sample size',ylab='Difference in relative abundance',cex.lab=1.5,cex.axis=1.2,bty='n',pch=16,col=alpha('black',0.5))


#############
# Get diversity metrics

maxes = c()
div.summary = data.frame('Site'=NA,'Year'=NA,'Decay.rate'=-999,'Dominance'=-999,'Fishers.alpha'=-999,'Richness'=-999,'Richness.rarefied'=-999,'Evenness.pei'=-999,'Evenness.evar'=-999,'Evenness.pei.rarefied'=-999,'Evenness.evar.rarefied'=-999,'Dominant.species1'=NA,'Dominant.species2'=NA,'Dominant.species3'=NA,'Party_Hours'=-999,'grid_id'=NA,'Latitude'=-999,'Longitude'=-999,'Julian.date'=NA,stringsAsFactors=F)
rx = 1
for (s in 1:length(sites)){
	print(s)
	dat = butterfly_counts[which(butterfly_counts$Site==sites[s]),]
	years = sort(as.numeric(unique(dat$Year)))
	max1 = c()
	for (y in 1:length(years)){
#		print(noquote(paste0('     ',years[y])))
		# Get curve metrics
		dat2 = dat[which(dat$Year==years[y]),]
		ranks = calc.ranks(dat2$Species,dat2$N.butterflies)
		div.summary[rx,1] = sites[s]
		div.summary[rx,2] = years[y]
		div.summary[rx,3] = ranks[1] #decay rate
		div.summary[rx,4] = ranks[2] #dominance
		div.summary[rx,5] = ranks[3] #fisher's alpha
		pos1 = which(dat2$N.butterflies>0)
		sp2 = dat2$Species[pos1]
		ab2 = dat2$N.butterflies[pos1]
		order1 = order(ab2,decreasing=T)
		sp3 = sp2[order1]
		ab3 = ab2[order1]
		ypos = which(1993:2018==years[y])
		div.summary[rx,6] = length(pos1) #richness
		div.summary[rx,7] = richness.rarefied[s,ypos] #rarefied richness
		div.summary[rx,8] = get.pei(ab2) #Pielou's evenness
		div.summary[rx,9] = get.evar(ab2) #Evar evenness
		div.summary[rx,10] = evenness.pei.rarefied[s,ypos] #rarefied Peilou's evenness
		div.summary[rx,11] = evenness.evar.rarefied[s,ypos] #rarefied Evar evenness
		div.summary[rx,12] = sp3[1] #1st dominant species
		div.summary[rx,13] = sp3[2] #2nd dominant species
		div.summary[rx,14] = sp3[3] #3rd dominant species
		div.summary[rx,15] = sum(unique(dat2$Party_Hours))
		div.summary[rx,16] = unique(dat2$grid_id)
		div.summary[rx,17] = unique(dat2$Latitude)
		div.summary[rx,18] = unique(dat2$Longitude)
		div.summary[rx,19] = paste0(unique(dat2$Julian.date),collapse=';')
		rx = rx + 1
		max1 = c(max1,ab3[1]/sum(ab3))
	}
	maxes = c(maxes,max(max1,na.rm=T))
}
write.table(div.summary,'site_diversity_summary_JunAug.txt',sep='\t',row.names=F,quote=F)

check = div.summary[grep('Anchorage',div.summary$Site),]

# Plot RACs
site.species.count = apply(array(sites),1,function(x){y=butterfly_counts[which(butterfly_counts$Site==x),];y2=y[which(y$N.butterflies>0),];length(unique(y2$Species))})
cols = viridis(26)
years1 = 1993:2018
for (s in 1:length(sites)){
	print(s)
	dat = butterfly_counts[which(butterfly_counts$Site==sites[s]),]
	years = sort(as.numeric(unique(dat$Year)))
	png(paste0('./plots/rank_abundance/RACs/JunAug/ranks_site',s,'.png'))
	par(oma=c(0,0,0,0),mar=c(5,6,2,1))
	plot(x=1:site.species.count[s],y=c(rep(0,site.species.count[s]-1),maxes[s]),col='white',xaxt='n',yaxt='n',xlab='Species rank',ylab='Relative abundance',cex.lab=2,bty='n',main=sites[s])
	axis(1,lwd=3,cex.axis=1.5,at=1:site.species.count[s])
	axis(2,lwd=3,cex.axis=1.5)
	legend('topright',legend=years1,col=cols,ncol=4,cex=1,bty='n',lwd=3,lty=1)
	rank2 = div.summary[which(div.summary$Site==sites[s]),]
	rich1.diff = rank2$Richness[nrow(rank2)] - rank2$Richness[1]
	rich2.diff = round(rank2$Richness.rarefied[nrow(rank2)] - rank2$Richness.rarefied[1],2)
	even.diff = round(rank2$Evenness[nrow(rank2)] - rank2$Evenness[1],2)
	legend(x=site.species.count[s]/2,y=maxes[s]/2,legend=c(paste0('Richness diff = ',rich1.diff),paste0('Richness rarefied diff = ',rich2.diff),paste0('Evenness diff = ',even.diff)),cex=1,bty='n')
	for (y in 1:length(years)){
		dat2 = dat[which(dat$Year==years[y]),]
		pos1 = which(dat2$N.butterflies>0)
		sp2 = dat2$Species[pos1]
		ab2 = dat2$N.butterflies[pos1]
		order1 = order(ab2,decreasing=T)
		sp3 = sp2[order1]
		ab3 = ab2[order1]
		lines(x=1:length(pos1),y=ab3/sum(ab3),lwd=4,col=cols[which(years==years[y])])
	}
	dev.off()
}

chart.Correlation(div.summary[,2:11], histogram=TRUE, pch=19)


# Plot RACs based on species probabilities
#dout = read.table('butterfly_abundance_probabilities.txt',sep='\t',as.is=T,check.names=F,header=T)
#str(dout)
#site.species.count = apply(array(sites),1,function(x){y=dout[which(dout$Site==x),];y2=y[which(y$Prob.abundance>0),];length(unique(y2$Species))})
#cols = viridis(26)
#years1 = 1993:2018
#for (s in 1:length(sites)){
#	print(s)
#	dat = dout[which(dout$Site==sites[s]),]
#	years = sort(as.numeric(unique(dat$Year)))
#	png(paste0('./plots/rank_abundance/RACs_prob/ranks_site',s,'.png'))
#	par(oma=c(0,0,0,0),mar=c(5,6,2,1))
#	plot(x=1:site.species.count[s],y=c(rep(0,site.species.count[s]-1),maxes[s]),col='white',xaxt='n',yaxt='n',xlab='Species rank',ylab='Relative abundance',cex.lab=2,bty='n',main=sites[s])
#	axis(1,lwd=3,cex.axis=1.5,at=1:site.species.count[s])
#	axis(2,lwd=3,cex.axis=1.5)
#	legend('topright',legend=years1,col=cols,ncol=4,cex=1,bty='n',lwd=3,lty=1)
#	rank2 = div.summary[which(div.summary$Site==sites[s]),]
#	rich1.diff = rank2$Richness[nrow(rank2)] - rank2$Richness[1]
#	rich2.diff = round(rank2$Richness.rarefied[nrow(rank2)] - rank2$Richness.rarefied[1],2)
#	even.diff = round(rank2$Evenness[nrow(rank2)] - rank2$Evenness[1],2)
#	legend(x=site.species.count[s]/2,y=maxes[s]/2,legend=c(paste0('Richness diff = ',rich1.diff),paste0('Richness rarefied diff = ',rich2.diff),paste0('Evenness diff = ',even.diff)),cex=1,bty='n')
#	for (y in 1:length(years)){
#		dat2 = dat[which(dat$Year==years[y]),]
#		pos1 = which(dat2$Prob.abundance>0)
#		sp2 = dat2$Species[pos1]
#		ab2 = dat2$Prob.abundance[pos1]
#		order1 = order(ab2,decreasing=T)
#		sp3 = sp2[order1]
#		ab3 = ab2[order1]
#		lines(x=1:length(pos1),y=ab3/sum(ab3),lwd=4,col=cols[which(years==years[y])])
#	}
#	dev.off()
#}


#####################################################
# Trends in rank abundance metrics (using AR_REML)

setwd('C:/Users/mcros/Desktop/Postdoc UGA/NABA_butterflies')

source('C:/Users/mcros/Desktop/Postdoc UGA/LTER_insect_diversity/code/AR_reml.R')
library(rgdal)
library(rgeos)
library(scales)
library(nlme)

# Import diversity metrics
con = file('site_diversity_summary_JunAug.txt','r')
add.lines = readLines(con)
close(con)
div.summary = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(div.summary) = strsplit(add.lines[1],'\t')[[1]]
div.summary$Decay.rate = as.numeric(div.summary$Decay.rate)
div.summary$Dominance = as.numeric(div.summary$Dominance)
div.summary$Fishers.alpha = as.numeric(div.summary$Fishers.alpha)
div.summary$Richness = as.numeric(div.summary$Richness)
div.summary$Evenness.pei = as.numeric(div.summary$Evenness.pei)
div.summary$Evenness.pei.rarefied = as.numeric(div.summary$Evenness.pei.rarefied)
div.summary$Evenness.evar = as.numeric(div.summary$Evenness.evar)
div.summary$Evenness.evar.rarefied = as.numeric(div.summary$Evenness.evar.rarefied)
div.summary$Richness.rarefied = as.numeric(div.summary$Richness.rarefied)
div.summary$Latitude = as.numeric(div.summary$Latitude)
div.summary$Longitude = as.numeric(div.summary$Longitude)
div.summary$Year = as.numeric(div.summary$Year)
div.summary$grid_id = as.integer(div.summary$grid_id)
div.summary$Party_Hours = as.numeric(div.summary$Party_Hours)

check = div.summary[grep('Anchorage',div.summary$Site),]

sites = unique(div.summary$Site)
years = 1993:2018
trend.summary = data.frame('Site'=NA,'Richness.trend'=-999,'Richness.rarefied.trend'=-999,'Evenness.pei.trend'=-999,'Evenness.pei.rarefied.trend'=-999,'Evenness.evar.trend'=-999,'Evenness.evar.rarefied.trend'=-999,'Dominance.trend'=-999,'Fishers.alpha.trend'=-999,'Decay.rate.trend'=-999,'Party_Hours.trend'=-999,'Series.length'=-999)
trend.summary2 = data.frame('Site'=NA,'Richness.trend'=-999,'Richness.rarefied.trend'=-999,'Evenness.pei.trend'=-999,'Evenness.pei.rarefied.trend'=-999,'Evenness.evar.trend'=-999,'Evenness.evar.rarefied.trend'=-999,'Dominance.trend'=-999,'Fishers.alpha.trend'=-999,'Decay.rate.trend'=-999,'Party_Hours.trend'=-999,'Series.length'=-999)
variables = c('Richness','Richness.rarefied','Evenness.pei','Evenness.pei.rarefied','Evenness.evar','Evenness.evar.rarefied','Dominance','Fishers.alpha','Decay.rate','Party_Hours')
thresh = 5
lost.sites = data.frame()
for (f in 1:length(sites)){
	print(f)
	dat1 = div.summary[which(div.summary$Site==sites[f]),]
	trend.summary[f,1] = trend.summary2[f,1] = sites[f]
	for (v in 1:length(variables)){
		vdat = dat1[,which(colnames(dat1)==variables[v])]
		if (length(which(!is.na(vdat)))>=thresh){ #remove time series shorter than 5 years
			if (length(unique(vdat))==1){ #values all equal, cannot estimate trend (=0)
				lost.sites = data.frame(rbind(lost.sites,data.frame('Variable'=variables[v],'Site'=sites[f],stringsAsFactors=F)),stringsAsFactors=F)
			} else {
				y.match = match(years,dat1$Year)
				X = y.match
				X[which(!is.na(X))] = vdat[X[which(!is.na(X))]] #fill-in values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				ph = Z
				ph[which(!is.na(Z))] = dat1$Party_Hours[which(!is.na(vdat))]
				t.scale = years - years[1]
	#			t.scale = 1:length(years); t.scale = (t.scale-min(t.scale))/max(t.scale) #original transform: scale between 0 and 1
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				trend.summary[f,v+1] = arr.Z$coef[2] #slope of trend line
				trend.summary[f,12] = length(t.scale)
				# Plot Z transform on t.scale fit
	#			png(paste0('./plots/time_trends/diversity/fitted/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.png'))
	#			plot(Z ~ t.scale, main=paste0(variables[v], ": b = ",round(arr.Z$b, digits=3)))
	#			curve(arr.Z$coef[1] + arr.Z$coef[2] * x, from=0, to=max(t.scale), add=T)
	#			dev.off()
				#Plot time series
				png(paste0('./plots/diversity_lineplots/JunAug/',variables[v],'_',f,'.png'))
				plot(t.scale,X,type='l',main=sites[f],xlab='Scaled time',ylab=variables[v],lwd=2)
	#			abline(a=arr.Z$coef[1],b=arr.Z$coef[2],lty=2,col='red',lwd=1.5)
				legend('topleft',legend=c(paste0('slope = ',round(arr.Z$coef[2],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(arr.Z$b,2)),paste0('logLik = ',round(arr.Z$logLik,2))),ncol=1,bty='n',adj=0)
				dev.off()

				# Generalized linear model, weighted by Party Hours
				if (variables[v]=='Evenness.evar.rarefied' | variables[v]=='Richness.rarefied'){
					data1 = data.frame('Z'=Z, 't.scale'=t.scale,'Party_Hours'=ph)
					gls.wt.Z <- gls(Z ~ t.scale, weights = varFixed(value = ~ exp(-Party_Hours/100)), correlation=corAR1(form = ~ 1), na.action = na.omit, data=data1) # Here is the weighting assuming the variance is proportional to Party Hours (which is like assuming Poisson sampling)
					trend.summary2[f,v+1] = gls.wt.Z$coefficients[2] #slope of trend line
					trend.summary2[f,12] = length(t.scale)
				} else {
				}
			}
		} else {
			lost.sites = data.frame(rbind(lost.sites,data.frame('Variable'=variables[v],'Site'=sites[f],stringsAsFactors=F)),stringsAsFactors=F)
		}
	}
}
write.table(trend.summary,'./diversity_trends_JunAug.txt',sep='\t',row.names=F)
write.table(trend.summary2,'./diversity_trends_wt_JunAug.txt',sep='\t',row.names=F)

length(unique(lost.sites$Site)) #21 sites with <5 years data for at least one metric

trend.summary = read.table('./diversity_trends_JunAug.txt',sep='\t',as.is=T,check.names=F,header=T)
trend.summary2 = read.table('./diversity_trends_wt_JunAug.txt',sep='\t',as.is=T,check.names=F,header=T)

length(which(trend.summary$Richness.rarefied.trend>0.05)) # 19
length(which(trend.summary$Richness.rarefied.trend>0.05)) / nrow(trend.summary) # 3.8%
length(which(trend.summary$Richness.rarefied.trend<(-0.05))) # 405
length(which(trend.summary$Richness.rarefied.trend<(-0.05))) / nrow(trend.summary) # 80.1%

length(which(trend.summary$Evenness.evar.rarefied.trend>0.05)) # 61
length(which(trend.summary$Evenness.evar.rarefied.trend>0.05)) / nrow(trend.summary) # 12.1%
length(which(trend.summary$Evenness.evar.rarefied.trend<(-0.05))) # 251
length(which(trend.summary$Evenness.evar.rarefied.trend<(-0.05))) / nrow(trend.summary) # 50.0%

library(scales)
par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(5,5,1,1))
plot(trend.summary$Richness.rarefied.trend,trend.summary2$Richness.rarefied.trend,xlab='Rarefied richness trend',ylab='Rarefied richness trend weighted',cex.lab=1.5,cex.axis=1.2,bty='n',pch=16,col=alpha('black',0.5))
abline(0,1)
plot(trend.summary$Evenness.evar.rarefied.trend,trend.summary2$Evenness.evar.rarefied.trend,xlab='Rarefied evenness trend',ylab='Rarefied evenness trend weighted',cex.lab=1.5,cex.axis=1.2,bty='n',pch=16,col=alpha('black',0.5))
abline(0,1)



# Split diversity trends by ecoregion
trend.summary = read.table('./diversity_trends_JunAug.txt',sep='\t',as.is=T,check.names=F,header=T)
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km_NoMerge.shp',layer='butterfly_sites_50km_NoMerge',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
trend.summary$EcoI = butterfly_circles@data$EcoI[match(trend.summary$Site,butterfly_circles@data$Site)]
trend.summary$EcorgnI = apply(array(trend.summary$EcoI),1,function(x){eco1@data[which(eco1@data[,1]==x)[1],2]})
ecos = unique(trend.summary$EcorgnI)
ecos = ecos[which(!is.na(ecos))]; ecos = ecos[which(ecos!='WATER')]
par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(5,5,1,1))
for (e in 1:length(ecos)){
	tdat = trend.summary[which(trend.summary$EcorgnI==ecos[e]),]
	hist(tdat$Evenness.evar.rarefied.trend,xlab='Evenness trend',breaks=seq(-0.7,0.7,0.1),ylab='No. sites',cex.lab=1.5,col='grey50',yaxt='n',xaxt='n',main=ecos[e],cex.main=1)
	axis(1,lwd=3,cex.axis=1.5,at=round(seq(-0.7,0.7,0.1),1))
	axis(2,lwd=3,cex.axis=1.5)
	abline(v=0,lwd=3,lty=1,col='pink')
}
par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(5,5,1,1))
for (e in 1:length(ecos)){
	tdat = trend.summary[which(trend.summary$EcorgnI==ecos[e]),]
	hist(tdat$Richness.rarefied.trend,xlab='Richness trend',breaks=seq(-0.7,0.7,0.1),ylab='No. sites',cex.lab=1.5,col='grey50',yaxt='n',xaxt='n',main=ecos[e],cex.main=1)
	axis(1,lwd=3,cex.axis=1.5,at=round(seq(-0.7,0.7,0.1),1))
	axis(2,lwd=3,cex.axis=1.5)
	abline(v=0,lwd=3,lty=1,col='pink')
}


##################################
# Jul set

butterfly_counts = read.table('./data/butterfly_data_NoMergeJul_w0s_m5.txt',sep='\t',as.is=T,check.names=F,header=T)

# Rarefy richness and evenness
sites = unique(butterfly_counts$Site)
species = sort(unique(butterfly_counts$Species))
years = 1993:2018

# Create list of community data matrices for each year
comm.mat.list = list()
keep.list = list()
for (y in 1:length(years)){
	print(years[y])
	mat1 = matrix(NA,nrow=length(sites),ncol=length(species))
	for (i in 1:length(sites)){
		pos1 = which(butterfly_counts$Year==years[y] & butterfly_counts$Site==sites[i])
		butterfly_counts2 = butterfly_counts[pos1,]
		mat1[i,] = butterfly_counts2$N.butterflies[match(species,butterfly_counts2$Species)]
	}
	na.count = apply(mat1,1,function(x){length(which(is.na(x)))})
	keep = which(na.count<length(species)) #position of sites that had any species counts in year y
	comm.mat.list[[y]] = mat1[keep,]
	keep.list[[y]] = keep
}

richness.rarefied = matrix(NA,nrow=length(sites),ncol=length(years))
for (y in 1:length(years)){
	print(years[y])
	mat1 = comm.mat.list[[y]]
	mat1[which(is.na(mat1))] = 0
	raremax <- min(rowSums(mat1))
	Srare <- rarefy(mat1, raremax)
#	sp.count = apply(mat1,1,function(x){length(which(x>0))})
#	plot(sp.count,Srare)
	keep = keep.list[[y]]
	richness.rarefied[keep,y] = Srare
}

library(scales)
evenness.pei.rarefied = matrix(NA,nrow=length(sites),ncol=length(years))
evenness.evar.rarefied = matrix(NA,nrow=length(sites),ncol=length(years))
for (y in 1:length(years)){
	print(years[y])
	mat1 = comm.mat.list[[y]]
	mat1[which(is.na(mat1))] = 0
	raremax <- min(rowSums(mat1))
	Erare.pei <- rarefy2(xx=mat1, sample1=raremax, divfun=get.pei)
	Erare.evar <- rarefy2(xx=mat1, sample1=raremax, divfun=get.evar)
	keep = keep.list[[y]]
	evenness.pei.rarefied[keep,y] = Erare.pei
	evenness.evar.rarefied[keep,y] = Erare.evar
}
plot(c(evenness.pei.rarefied),c(evenness.evar.rarefied),ylim=c(0,1),bty='n',xlab="Pielou's evenness rarefied",ylab="Evar rarefied",cex.lab=1.5,cex.axis=1.2,pch=16,col=alpha('black',0.5))


# Explore difference between raw and rarefied abundance
#dout = data.frame()
#for (y in 1:length(years)){ #takes a LONG time
#	print(years[y])
#	keep = keep.list[[y]]
#	sites1 = sites[keep]
#	mat1 = comm.mat.list[[y]]
#	mat1[which(is.na(mat1))] = 0
#	raremax <- min(rowSums(mat1))
#	for (i in 1:nrow(mat1)){ #for each site
#		y1 = mat1[i,][mat1[i,]>1]
#		J <- sum(y1)
#        ldiv <- lchoose(J, raremax)
#        p1 <- ifelse(J - y1 < raremax, 0, exp(lchoose(J - y1, raremax) - ldiv)) #p1 is larger for rarer individuals
#		p2 = 1-p1
#		p2sum = sum(p2)
#		p3 = mat1[i,]
#		p3[mat1[i,]>1] = p2
#		add1 = data.frame('Site'=rep(sites1[i],ncol(mat1)),'Year'=rep(years[y],ncol(mat1)),'Species'=species,'Min.sample'=rep(raremax,ncol(mat1)),'Raw.abundance'=mat1[i,],'Prob.abundance'=p3,'Sum.raw.abundance'=rep(J,ncol(mat1)),'Sum.prob.abundance'=rep(p2sum,ncol(mat1)),stringsAsFactors=F)
#		dout = data.frame(rbind(dout,add1),stringsAsFactors=F)
#	}
#}
#write.table(dout,'butterfly_abundance_probabilities.txt',sep='\t',row.names=F)
#plot(dout$Min.sample,(dout$Prob.abundance/dout$Sum.prob.abundance)-(dout$Raw.abundance/dout$Sum.raw.abundance),xlab='Minimum sample size',ylab='Difference in relative abundance',cex.lab=1.5,cex.axis=1.2,bty='n',pch=16,col=alpha('black',0.5))


#############
# Get diversity metrics

maxes = c()
div.summary = data.frame('Site'=NA,'Year'=NA,'Decay.rate'=-999,'Dominance'=-999,'Fishers.alpha'=-999,'Richness'=-999,'Richness.rarefied'=-999,'Evenness.pei'=-999,'Evenness.evar'=-999,'Evenness.pei.rarefied'=-999,'Evenness.evar.rarefied'=-999,'Dominant.species1'=NA,'Dominant.species2'=NA,'Dominant.species3'=NA,'Party_Hours'=-999,'grid_id'=NA,'Latitude'=-999,'Longitude'=-999,'Julian.date'=NA,stringsAsFactors=F)
rx = 1
for (s in 1:length(sites)){
	print(s)
	dat = butterfly_counts[which(butterfly_counts$Site==sites[s]),]
	years = sort(as.numeric(unique(dat$Year)))
	max1 = c()
	for (y in 1:length(years)){
#		print(noquote(paste0('     ',years[y])))
		# Get curve metrics
		dat2 = dat[which(dat$Year==years[y]),]
		ranks = calc.ranks(dat2$Species,dat2$N.butterflies)
		div.summary[rx,1] = sites[s]
		div.summary[rx,2] = years[y]
		div.summary[rx,3] = ranks[1] #decay rate
		div.summary[rx,4] = ranks[2] #dominance
		div.summary[rx,5] = ranks[3] #fisher's alpha
		pos1 = which(dat2$N.butterflies>0)
		sp2 = dat2$Species[pos1]
		ab2 = dat2$N.butterflies[pos1]
		order1 = order(ab2,decreasing=T)
		sp3 = sp2[order1]
		ab3 = ab2[order1]
		ypos = which(1993:2018==years[y])
		div.summary[rx,6] = length(pos1) #richness
		div.summary[rx,7] = richness.rarefied[s,ypos] #rarefied richness
		div.summary[rx,8] = get.pei(ab2) #Pielou's evenness
		div.summary[rx,9] = get.evar(ab2) #Evar evenness
		div.summary[rx,10] = evenness.pei.rarefied[s,ypos] #rarefied Peilou's evenness
		div.summary[rx,11] = evenness.evar.rarefied[s,ypos] #rarefied Evar evenness
		div.summary[rx,12] = sp3[1] #1st dominant species
		div.summary[rx,13] = sp3[2] #2nd dominant species
		div.summary[rx,14] = sp3[3] #3rd dominant species
		div.summary[rx,15] = sum(unique(dat2$Party_Hours))
		div.summary[rx,16] = unique(dat2$grid_id)
		div.summary[rx,17] = unique(dat2$Latitude)
		div.summary[rx,18] = unique(dat2$Longitude)
		div.summary[rx,19] = paste0(unique(dat2$Julian.date),collapse=';')
		rx = rx + 1
		max1 = c(max1,ab3[1]/sum(ab3))
	}
	maxes = c(maxes,max(max1,na.rm=T))
}
write.table(div.summary,'site_diversity_summary_Jul.txt',sep='\t',row.names=F,quote=F)

check = div.summary[grep('Anchorage',div.summary$Site),]

# Plot RACs
site.species.count = apply(array(sites),1,function(x){y=butterfly_counts[which(butterfly_counts$Site==x),];y2=y[which(y$N.butterflies>0),];length(unique(y2$Species))})
cols = viridis(26)
years1 = 1993:2018
for (s in 1:length(sites)){
	print(s)
	dat = butterfly_counts[which(butterfly_counts$Site==sites[s]),]
	years = sort(as.numeric(unique(dat$Year)))
	png(paste0('./plots/rank_abundance/RACs/Jul/ranks_site',s,'.png'))
	par(oma=c(0,0,0,0),mar=c(5,6,2,1))
	plot(x=1:site.species.count[s],y=c(rep(0,site.species.count[s]-1),maxes[s]),col='white',xaxt='n',yaxt='n',xlab='Species rank',ylab='Relative abundance',cex.lab=2,bty='n',main=sites[s])
	axis(1,lwd=3,cex.axis=1.5,at=1:site.species.count[s])
	axis(2,lwd=3,cex.axis=1.5)
	legend('topright',legend=years1,col=cols,ncol=4,cex=1,bty='n',lwd=3,lty=1)
	rank2 = div.summary[which(div.summary$Site==sites[s]),]
	rich1.diff = rank2$Richness[nrow(rank2)] - rank2$Richness[1]
	rich2.diff = round(rank2$Richness.rarefied[nrow(rank2)] - rank2$Richness.rarefied[1],2)
	even.diff = round(rank2$Evenness[nrow(rank2)] - rank2$Evenness[1],2)
	legend(x=site.species.count[s]/2,y=maxes[s]/2,legend=c(paste0('Richness diff = ',rich1.diff),paste0('Richness rarefied diff = ',rich2.diff),paste0('Evenness diff = ',even.diff)),cex=1,bty='n')
	for (y in 1:length(years)){
		dat2 = dat[which(dat$Year==years[y]),]
		pos1 = which(dat2$N.butterflies>0)
		sp2 = dat2$Species[pos1]
		ab2 = dat2$N.butterflies[pos1]
		order1 = order(ab2,decreasing=T)
		sp3 = sp2[order1]
		ab3 = ab2[order1]
		lines(x=1:length(pos1),y=ab3/sum(ab3),lwd=4,col=cols[which(years==years[y])])
	}
	dev.off()
}

chart.Correlation(div.summary[,2:11], histogram=TRUE, pch=19)


# Plot RACs based on species probabilities
#dout = read.table('butterfly_abundance_probabilities.txt',sep='\t',as.is=T,check.names=F,header=T)
#str(dout)
#site.species.count = apply(array(sites),1,function(x){y=dout[which(dout$Site==x),];y2=y[which(y$Prob.abundance>0),];length(unique(y2$Species))})
#cols = viridis(26)
#years1 = 1993:2018
#for (s in 1:length(sites)){
#	print(s)
#	dat = dout[which(dout$Site==sites[s]),]
#	years = sort(as.numeric(unique(dat$Year)))
#	png(paste0('./plots/rank_abundance/RACs_prob/ranks_site',s,'.png'))
#	par(oma=c(0,0,0,0),mar=c(5,6,2,1))
#	plot(x=1:site.species.count[s],y=c(rep(0,site.species.count[s]-1),maxes[s]),col='white',xaxt='n',yaxt='n',xlab='Species rank',ylab='Relative abundance',cex.lab=2,bty='n',main=sites[s])
#	axis(1,lwd=3,cex.axis=1.5,at=1:site.species.count[s])
#	axis(2,lwd=3,cex.axis=1.5)
#	legend('topright',legend=years1,col=cols,ncol=4,cex=1,bty='n',lwd=3,lty=1)
#	rank2 = div.summary[which(div.summary$Site==sites[s]),]
#	rich1.diff = rank2$Richness[nrow(rank2)] - rank2$Richness[1]
#	rich2.diff = round(rank2$Richness.rarefied[nrow(rank2)] - rank2$Richness.rarefied[1],2)
#	even.diff = round(rank2$Evenness[nrow(rank2)] - rank2$Evenness[1],2)
#	legend(x=site.species.count[s]/2,y=maxes[s]/2,legend=c(paste0('Richness diff = ',rich1.diff),paste0('Richness rarefied diff = ',rich2.diff),paste0('Evenness diff = ',even.diff)),cex=1,bty='n')
#	for (y in 1:length(years)){
#		dat2 = dat[which(dat$Year==years[y]),]
#		pos1 = which(dat2$Prob.abundance>0)
#		sp2 = dat2$Species[pos1]
#		ab2 = dat2$Prob.abundance[pos1]
#		order1 = order(ab2,decreasing=T)
#		sp3 = sp2[order1]
#		ab3 = ab2[order1]
#		lines(x=1:length(pos1),y=ab3/sum(ab3),lwd=4,col=cols[which(years==years[y])])
#	}
#	dev.off()
#}


#####################################################
# Trends in rank abundance metrics (using AR_REML)

setwd('C:/Users/mcros/Desktop/Postdoc UGA/NABA_butterflies')

source('C:/Users/mcros/Desktop/Postdoc UGA/LTER_insect_diversity/code/AR_reml.R')
library(rgdal)
library(rgeos)
library(scales)
library(nlme)

# Import diversity metrics
con = file('site_diversity_summary_Jul.txt','r')
add.lines = readLines(con)
close(con)
div.summary = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(div.summary) = strsplit(add.lines[1],'\t')[[1]]
div.summary$Decay.rate = as.numeric(div.summary$Decay.rate)
div.summary$Dominance = as.numeric(div.summary$Dominance)
div.summary$Fishers.alpha = as.numeric(div.summary$Fishers.alpha)
div.summary$Richness = as.numeric(div.summary$Richness)
div.summary$Evenness.pei = as.numeric(div.summary$Evenness.pei)
div.summary$Evenness.pei.rarefied = as.numeric(div.summary$Evenness.pei.rarefied)
div.summary$Evenness.evar = as.numeric(div.summary$Evenness.evar)
div.summary$Evenness.evar.rarefied = as.numeric(div.summary$Evenness.evar.rarefied)
div.summary$Richness.rarefied = as.numeric(div.summary$Richness.rarefied)
div.summary$Latitude = as.numeric(div.summary$Latitude)
div.summary$Longitude = as.numeric(div.summary$Longitude)
div.summary$Year = as.numeric(div.summary$Year)
div.summary$grid_id = as.integer(div.summary$grid_id)
div.summary$Party_Hours = as.numeric(div.summary$Party_Hours)

check = div.summary[grep('Anchorage',div.summary$Site),]

sites = unique(div.summary$Site)
years = 1993:2018
trend.summary = data.frame('Site'=NA,'Richness.trend'=-999,'Richness.rarefied.trend'=-999,'Evenness.pei.trend'=-999,'Evenness.pei.rarefied.trend'=-999,'Evenness.evar.trend'=-999,'Evenness.evar.rarefied.trend'=-999,'Dominance.trend'=-999,'Fishers.alpha.trend'=-999,'Decay.rate.trend'=-999,'Party_Hours.trend'=-999,'Series.length'=-999)
trend.summary2 = data.frame('Site'=NA,'Richness.trend'=-999,'Richness.rarefied.trend'=-999,'Evenness.pei.trend'=-999,'Evenness.pei.rarefied.trend'=-999,'Evenness.evar.trend'=-999,'Evenness.evar.rarefied.trend'=-999,'Dominance.trend'=-999,'Fishers.alpha.trend'=-999,'Decay.rate.trend'=-999,'Party_Hours.trend'=-999,'Series.length'=-999)
variables = c('Richness','Richness.rarefied','Evenness.pei','Evenness.pei.rarefied','Evenness.evar','Evenness.evar.rarefied','Dominance','Fishers.alpha','Decay.rate','Party_Hours')
thresh = 5
lost.sites = data.frame()
for (f in 1:length(sites)){
	print(f)
	dat1 = div.summary[which(div.summary$Site==sites[f]),]
	trend.summary[f,1] = trend.summary2[f,1] = sites[f]
	for (v in 1:length(variables)){
		vdat = dat1[,which(colnames(dat1)==variables[v])]
		if (length(which(!is.na(vdat)))>=thresh){ #remove time series shorter than 5 years
			if (length(unique(vdat))==1){ #values all equal, cannot estimate trend (=0)
				lost.sites = data.frame(rbind(lost.sites,data.frame('Variable'=variables[v],'Site'=sites[f],stringsAsFactors=F)),stringsAsFactors=F)
			} else {
				y.match = match(years,dat1$Year)
				X = y.match
				X[which(!is.na(X))] = vdat[X[which(!is.na(X))]] #fill-in values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				ph = Z
				ph[which(!is.na(Z))] = dat1$Party_Hours[which(!is.na(vdat))]
				t.scale = years - years[1]
	#			t.scale = 1:length(years); t.scale = (t.scale-min(t.scale))/max(t.scale) #original transform: scale between 0 and 1
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				trend.summary[f,v+1] = arr.Z$coef[2] #slope of trend line
				trend.summary[f,12] = length(t.scale)
				# Plot Z transform on t.scale fit
	#			png(paste0('./plots/time_trends/diversity/fitted/',variables[v],'_',unique(dat1$LTER.site),'_',uL[l],'.png'))
	#			plot(Z ~ t.scale, main=paste0(variables[v], ": b = ",round(arr.Z$b, digits=3)))
	#			curve(arr.Z$coef[1] + arr.Z$coef[2] * x, from=0, to=max(t.scale), add=T)
	#			dev.off()
				#Plot time series
				png(paste0('./plots/diversity_lineplots/JunAug/',variables[v],'_',f,'.png'))
				plot(t.scale,X,type='l',main=sites[f],xlab='Scaled time',ylab=variables[v],lwd=2)
	#			abline(a=arr.Z$coef[1],b=arr.Z$coef[2],lty=2,col='red',lwd=1.5)
				legend('topleft',legend=c(paste0('slope = ',round(arr.Z$coef[2],2)),paste0('length = ',length(which(!is.na(Z)))),paste0('autocor = ',round(arr.Z$b,2)),paste0('logLik = ',round(arr.Z$logLik,2))),ncol=1,bty='n',adj=0)
				dev.off()

				# Generalized linear model, weighted by Party Hours
				if (variables[v]=='Evenness.evar.rarefied'){
					data1 = data.frame('Z'=Z, 't.scale'=t.scale,'Party_Hours'=ph)
					gls.wt.Z <- gls(Z ~ t.scale, weights = varFixed(value = ~ exp(-Party_Hours/100)), correlation=corAR1(form = ~ 1), na.action = na.omit, data=data1) # Here is the weighting assuming the variance is proportional to Party Hours (which is like assuming Poisson sampling)
					trend.summary2[f,v+1] = gls.wt.Z$coefficients[2] #slope of trend line
					trend.summary2[f,13] = length(t.scale)
				} else {
				}
			}
		} else {
			lost.sites = data.frame(rbind(lost.sites,data.frame('Variable'=variables[v],'Site'=sites[f],stringsAsFactors=F)),stringsAsFactors=F)
		}
	}
}
write.table(trend.summary,'./diversity_trends_Jul.txt',sep='\t',row.names=F)
write.table(trend.summary2,'./diversity_trends_wt_Jul.txt',sep='\t',row.names=F)

length(unique(lost.sites$Site)) #21 sites with <5 years data for at least one metric

trend.summary = read.table('./diversity_trends_Jul.txt',sep='\t',as.is=T,check.names=F,header=T)
trend.summary2 = read.table('./diversity_trends_wt_Jul.txt',sep='\t',as.is=T,check.names=F,header=T)

length(which(trend.summary$Richness.rarefied.trend>0.05)) # 36
length(which(trend.summary$Richness.rarefied.trend>0.05)) / nrow(trend.summary) # 9.8%
length(which(trend.summary$Richness.rarefied.trend<(-0.05))) # 215
length(which(trend.summary$Richness.rarefied.trend<(-0.05))) / nrow(trend.summary) # 58.4%

length(which(trend.summary$Evenness.evar.rarefied.trend>0.05)) # 55
length(which(trend.summary$Evenness.evar.rarefied.trend>0.05)) / nrow(trend.summary) # 14.9%
length(which(trend.summary$Evenness.evar.rarefied.trend<(-0.05))) # 162
length(which(trend.summary$Evenness.evar.rarefied.trend<(-0.05))) / nrow(trend.summary) # 44.0%


# Create shapefile with diversity trends
library(classInt)
library(sf)
library(rgdal)
library(rgeos)
library(ggplot2)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

trend.summary = read.table('./diversity_trends_JunAug.txt',sep='\t',as.is=T,check.names=F,header=T)
proj1 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km_NoMerge.shp',layer='butterfly_sites_50km_NoMerge',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with grid
match1 = match(butterfly_circles@data$Site,trend.summary$Site)
trend.summary$Lat = butterfly_circles@coords[match1[which(!is.na(match1))],2]
trend.summary$Lon = butterfly_circles@coords[match1[which(!is.na(match1))],1]

dat = spTransform(SpatialPointsDataFrame(coords=cbind(trend.summary$Lon,trend.summary$Lat),data=trend.summary,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')),CRS(proj1))
writeOGR(dat,dsn='./shapefiles/butterfly_trends_JunAug.shp',layer='butterfly_trends_JunAug',driver='ESRI Shapefile',overwrite=T)


