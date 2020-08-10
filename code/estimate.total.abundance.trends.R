

setwd('your/working/directory')

library(rgdal)
library(sp)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(rgeos)
library(maptools)
library(spdep)
library(vegan)
library(inlabru)
library(INLA)
library(brinla)

options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

# plot theme
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# map theme
theme_map <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.background = element_blank(),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text.align=1,
          legend.text = element_text(size=9),
          legend.title=element_text(hjust=0, size=11),
          legend.justification=c(0, 0.5),
          plot.title = element_text(size=14, hjust = 0.7))
}

# Define function for calculating Pielou's evenness
calc.evenness = function(sp1,ab1){
	#sp1 = a vector of species names (species records do not have to be unique)
	#ab1 = a vector of corresponding abundances of species in sp1
	usp = unique(sp1) #unique species names
	spc = apply(array(usp),1,function(x){sum(ab1[which(sp1==x)],na.rm=T)}) #species counts
	H = diversity(spc)
	pei = H / log(specnumber(spc))
	return(pei)
}


################################################################################
# 50km scale (no pseudo-absences)

# Import butterfly data
butterfly_counts1 = read.table('./data/butterfly_data_gridded_50km_1993-2018_wNAs.txt',sep='\t',as.is=T,check.names=F,header=T)
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km.shp',layer='butterfly_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
butterfly_grid = readOGR(dsn='./shapefiles/butterfly_grid_50km.shp',layer='butterfly_grid_50km',verbose=F,stringsAsFactors=F) #shapefile with grid
grid50 <- as(butterfly_grid, "sf")
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

# Aggregate butterfly abundance across species per site
butterfly_counts = data.frame()
sites = unique(butterfly_counts1$Site)
for (s in 1:length(sites)){
	print(s)
	sdat = butterfly_counts1[which(butterfly_counts1$Site==sites[s]),]
	years = unique(sdat$Year)
	for (y in 1:length(years)){
		sdat2 = sdat[which(sdat$Year==years[y]),]
		evenness = calc.evenness(sdat2$Species,sdat2$N.butterflies)
		N.species = length(unique(sdat2$Species[which(sdat2$N.butterflies>0 & !is.na(sdat2$N.butterflies))]))
		N.butterflies = sum(sdat2$N.butterflies,na.rm=T)
		Party_Hours = mean(sdat2$Party_Hours)
		Longitude = unique(sdat2$Longitude)[1]
		Latitude = unique(sdat2$Latitude)[1]
		grid_id = unique(sdat2$grid_id)
		EcoI = unique(sdat2$EcoI)
		EcoII = unique(sdat2$EcoII)
		EcoregionI = unique(sdat2$EcoregionI)
		EcoregionII = unique(sdat2$EcoregionII)
		add.dat = data.frame('Site'=sites[s],'Year'=years[y],'N.species'=N.species,'Evenness'=evenness,'N.butterflies'=N.butterflies,'Party_Hours'=Party_Hours,'Longitude'=Longitude,
			'Latitude'=Latitude,'grid_id'=grid_id,'EcoI'=EcoI,'EcoII'=EcoII,'EcoregionI'=EcoregionI,'EcoregionII'=EcoregionII,stringsAsFactors=F)
		butterfly_counts = data.frame(rbind(butterfly_counts,add.dat),stringsAsFactors=F)
	}
}
write.table(butterfly_counts,'butterfly_total_abundance_50km_wNAs.txt',sep='\t',row.names=F)


# INLA model
sp_counts = butterfly_counts
sp_sites = unique(sp_counts$Site)
sp.years = sort(as.numeric(unique(sp_counts$Year)))

# Create spatial dataframe
sp_counts.shp = spTransform(SpatialPointsDataFrame(coords=cbind(sp_counts$Longitude,sp_counts$Latitude),data=sp_counts,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')),CRS(proj1))
sp_sites = unique(sp_counts$Site)
sp_circles = butterfly_circles[which(!is.na(match(butterfly_circles$Site,sp_sites))),]

# Summarize circles per cell (used for mapping later)
gids = sort(unique(sp_counts.shp@data$grid_id))
gids_count = apply(array(gids),1,function(x){gd = sp_counts.shp[which(sp_counts.shp$grid_id==x),];length(unique(gd$Site))})
circles_per_cell = data.frame('grid_id'=gids,'number_circles'=gids_count)

# Transform effort & standardize years
sp_counts$ln.Party_Hours = log(as.numeric(sp_counts$Party_Hours))
sp_counts$std_Year = sp_counts$Year - max(sp_counts$Year)

# Index and sort
sp_counts$eps_i <- sp_counts$alpha_i <- sp_counts$grid_id
sp_counts$tau_i <- sp_counts$eps_i
sp_counts$kappa_k <- as.integer(factor(sp_counts$Site))
sp_counts <- arrange(sp_counts, grid_id, std_Year)
n_circs <- max(sp_counts$kappa_k, na.rm=T)
n_cells <- max(sp_counts$alpha_i, na.rm=T)

# Make negative binomial model with raw counts
form1 <- N.butterflies ~ -1 + # remove grand mean
  # cell ICAR random intercepts
  f(alpha_i, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # cell ICAR random effort slopes
  f(eps_i, ln.Party_Hours, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # cell ICAR random year slopes
  f(tau_i, std_Year, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # random circle intercepts
  f(kappa_k, model="iid", constr=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

# Run model
print(noquote('---running model---'))
out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3)

# Overdispersion parameter
sm1 = summary(out1)[[3]]

# Vector of grid IDs of cells with counts
cells_with_counts = unique(sp_counts$grid_id)

# get alpha summaries
alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
alph_iw <- alph_ul - alph_ll
#	par(mfrow=c(1,3))
#	hist(alph); summary(alph)
#	hist(alph_ll); summary(alph_ll)
#	hist(alph_ul); summary(alph_ul)

# get epsilon summaries
eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
eps_iw <- eps_ul - eps_ll
#	par(mfrow=c(1,3))
#	hist(eps); summary(eps); round(sum(eps<1)/length(eps), 2)
#	hist(eps_ll); summary(eps_ll)
#	hist(eps_ul); summary(eps_ul)

c1 = cor.test(eps, alph, method="spearman")

# get tau summaries
tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
tau_iw <- tau_ul - tau_ll
#	par(mfrow=c(2,2))
#	hist(tau); summary(tau); round(sum(tau>=0)/length(tau), 2)
#	hist(tau_ll); summary(tau_ll)
#	hist(tau_ul); summary(tau_ul)
#	hist(tau_iw); summary(tau_iw)

# Correlation between epsilon and alpha
c2 = cor.test(alph, tau, method="spearman")

# Visualize goodness of fit with Probability Integral Transform
#	sum(out1$cpo$failure, na.rm=T) -2 * sum(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
pit1 <- data.frame(PIT=out1$cpo$pit) %>%
  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
pit2 <- ggplot(data=pit1, aes(x=PIT)) +
  geom_histogram(col="white") +
  xlab("Probability integral transform (PIT)") +
  ylab("Count"); pit2; summary(pit1$PIT)
#	ggsave(paste0(species[i],'_pit.png'),plot=pit2,device='png',path='./plots/PIT/',width=8,height=8,units="in",dpi=300)
ggsave('Allspecies_pit_1993-2018_50km_wNAs.png',plot=pit2,device='png',path='./plots/PIT/',width=8,height=8,units="in",dpi=300)

# collect posterior summaries into one dataframe
grid_key <- unique(sp_counts[, c("grid_id", "EcoI", "EcoII")])	
post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,eps, eps_ll, eps_ul, eps_iw, eps_sig=NA,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
post_sum$eps_sig <- ifelse((post_sum$eps_ll < 0 & post_sum$eps_ul > 0),post_sum$eps_sig <- NA,post_sum$eps_sig <- post_sum$eps)
post_sum$tau_sig <- ifelse((post_sum$tau_ll < 0 & post_sum$tau_ul > 0),post_sum$tau_sig <- NA,post_sum$tau_sig <- post_sum$tau)
post_sum$EcoI <- grid_key$EcoI[match(cells_with_counts,grid_key$grid_id)]
post_sum$EcoII <- grid_key$EcoII[match(cells_with_counts,grid_key$grid_id)]
# Add medians of slope estimates to data frame
sm2 = summary(post_sum)
sm2[3,] = trimws(gsub('Median :','',sm2[3,]),which='both')

# Make cell level maps
results_cells <- merge(butterfly_grid, post_sum, by="grid_id", all=F)
write.table(results_cells@data,'Allspecies_trends_1993-2018_50km_wNAs.txt',sep='\t',quote=F,row.names=F) # Write model output to file
res_sf <- as(results_cells, "sf")

# grid
grid2 <- merge(grid50, circles_per_cell, all=T)
grid2$number_circles[is.na(grid2$number_circles)] <- 0
grid_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=grid2, aes(fill=number_circles), col="gray40", size=0.3) +
  scale_fill_gradient("Circles\nper cell  ", low = ("white"),
					  high = ("red4"), space = "Lab", trans="log1p",
					  na.value = "grey40", guide = "colourbar",
					  breaks=c(0,2,7,20)) +
  theme_map() +
  theme(panel.grid.major=element_line(colour="transparent"))

# map tau
tau_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
					   mid = "white",
					   high = ("royalblue4"), midpoint = 0, space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p2 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_iw), col="gray40", size=0.3) +
  scale_fill_gradient2("Abund.\ntrend\ncredible\ninterval\nwidth\n(%/year)",
					   low = ("purple4"), mid = "white",
					   high = ("green4"), midpoint = median(res_sf$tau_iw), space = "Lab",
					   na.value = "grey40", guide = "colourbar",
					   breaks=c(3,6,9,12)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p3 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_sig), col="gray40", size=0.3) +
  scale_fill_gradient2("Sig.\nabund.\ntrend\n(%/year)",
					   low = muted("red4"), mid = "gray95",
					   high = muted("royalblue4"), midpoint = 0, space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# map epsilon
eps_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
  scale_fill_gradient2("Sampling\neffort", low = muted("purple4"), mid = "white",
					   high = muted("green4"), midpoint = median(res_sf$eps), space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# map alpha
alph_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=alph), col="gray40", size=0.3) +
  scale_fill_gradient2("Relative\nabund.", low = "tan4", mid = "white",
					   high = "green4", midpoint = (max(res_sf$alph) + min(res_sf$alph)) / 2, space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# print cell maps
print(noquote('---drawing maps---'))
ggsave('Allspecies_effort_1993-2018_50km_wNAs.png',plot =eps_p1,device='png',path='./',width=8,height=8,units="in",dpi=300)
ggsave('Allspecies_change_per_year_1993-2018_50km_wNAs.png',plot =tau_p1,device='png',path='./',width=8,height=8,units="in",dpi=300)
ggsave('Allspecies_relative_abundance_1993-2018_50km_wNAs.png',plot =alph_p1,device='png',path='./',width=8,height=8,units="in",dpi=300)

res_sf = read.table('Allspecies_trends_1993-2018_50km_wNAs.txt',sep='\t',as.is=T,check.names=F,header=T)
par(oma=c(0,0,0,0),mar=c(0,5,0,0))
boxplot(res_sf$tau,ylab='Total abundance trend',ylim=c(-10,10),cex.lab=2,boxlwd=2,staplelwd=2,whisklwd=2,whisklty=1,col='grey40',outpch=16,frame=F,yaxt='n')
axis(2,lwd=3,at=seq(-10,10,2),cex.axis=1.5)
abline(a=0,b=0,lwd=2,col='red',lty=2)
legend('topright',legend=c('No. sites = 497','No. cells = 414'),bty='n',cex=2,ncol=1)

dev.off()
hist(res_sf$tau,main='',ylim=c(0,200),breaks=seq(-12,6,2),ylab='No. species',xlab='Median total abundance trend',cex.lab=1.5,cex.axis=1.2)
abline(v=0,lwd=2,col='red')


############################
# Repeat INLA for total abundance (no pseudoabsences),
# but this time excluding the four most abundant species

# Import butterfly data
butterfly_counts1 = read.table('./data/butterfly_data_gridded_50km_1993-2018_wNAs.txt',sep='\t',as.is=T,check.names=F,header=T)
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km.shp',layer='butterfly_sites_50km',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
butterfly_grid = readOGR(dsn='./shapefiles/butterfly_grid_50km.shp',layer='butterfly_grid_50km',verbose=F,stringsAsFactors=F) #shapefile with grid
grid50 <- as(butterfly_grid, "sf")
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

species = unique(butterfly_counts1$Species)
scount = apply(array(species),1,function(x){sum(butterfly_counts1$N.butterflies[which(butterfly_counts1$Species==x)],na.rm=T)})
stab = cbind(species,scount)
stab[order(scount),]
butterfly_counts1 = butterfly_counts1[which(butterfly_counts1$Species!='Thymelicus lineola' & butterfly_counts1$Species!='Pieris rapae' & butterfly_counts1$Species!='Phyciodes tharos' & butterfly_counts1$Species!='Colias philodice'),]
species = unique(butterfly_counts1$Species)
scount = apply(array(species),1,function(x){sum(butterfly_counts1$N.butterflies[which(butterfly_counts1$Species==x)],na.rm=T)})
stab = cbind(species,scount)
stab[order(scount),]

# Aggregate butterfly abundance across species per site
butterfly_counts = data.frame()
sites = unique(butterfly_counts1$Site)
for (s in 1:length(sites)){
	print(s)
	sdat = butterfly_counts1[which(butterfly_counts1$Site==sites[s]),]
	years = unique(sdat$Year)
	for (y in 1:length(years)){
		sdat2 = sdat[which(sdat$Year==years[y]),]
		evenness = calc.evenness(sdat2$Species,sdat2$N.butterflies)
		N.species = length(unique(sdat2$Species[which(sdat2$N.butterflies>0 & !is.na(sdat2$N.butterflies))]))
		N.butterflies = sum(sdat2$N.butterflies,na.rm=T)
		Party_Hours = mean(sdat2$Party_Hours)
		Longitude = unique(sdat2$Longitude)[1]
		Latitude = unique(sdat2$Latitude)[1]
		grid_id = unique(sdat2$grid_id)
		EcoI = unique(sdat2$EcoI)
		EcoII = unique(sdat2$EcoII)
		EcoregionI = unique(sdat2$EcoregionI)
		EcoregionII = unique(sdat2$EcoregionII)
		add.dat = data.frame('Site'=sites[s],'Year'=years[y],'N.species'=N.species,'Evenness'=evenness,'N.butterflies'=N.butterflies,'Party_Hours'=Party_Hours,'Longitude'=Longitude,
			'Latitude'=Latitude,'grid_id'=grid_id,'EcoI'=EcoI,'EcoII'=EcoII,'EcoregionI'=EcoregionI,'EcoregionII'=EcoregionII,stringsAsFactors=F)
		butterfly_counts = data.frame(rbind(butterfly_counts,add.dat),stringsAsFactors=F)
	}
}
write.table(butterfly_counts,'butterfly_total_abundance_50km_wNAs_noTop4.txt',sep='\t',row.names=F)

# INLA model
sp_counts = butterfly_counts
sp_sites = unique(sp_counts$Site)
sp.years = sort(as.numeric(unique(sp_counts$Year)))

# Create spatial dataframe
sp_counts.shp = spTransform(SpatialPointsDataFrame(coords=cbind(sp_counts$Longitude,sp_counts$Latitude),data=sp_counts,proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')),CRS(proj1))
sp_sites = unique(sp_counts$Site)
sp_circles = butterfly_circles[which(!is.na(match(butterfly_circles$Site,sp_sites))),]

# Summarize circles per cell (used for mapping later)
gids = sort(unique(sp_counts.shp@data$grid_id))
gids_count = apply(array(gids),1,function(x){gd = sp_counts.shp[which(sp_counts.shp$grid_id==x),];length(unique(gd$Site))})
circles_per_cell = data.frame('grid_id'=gids,'number_circles'=gids_count)

# Transform effort & standardize years
sp_counts$ln.Party_Hours = log(as.numeric(sp_counts$Party_Hours))
sp_counts$std_Year = sp_counts$Year - max(sp_counts$Year)

# Index and sort
sp_counts$eps_i <- sp_counts$alpha_i <- sp_counts$grid_id
sp_counts$tau_i <- sp_counts$eps_i
sp_counts$kappa_k <- as.integer(factor(sp_counts$Site))
sp_counts <- arrange(sp_counts, grid_id, std_Year)
n_circs <- max(sp_counts$kappa_k, na.rm=T)
n_cells <- max(sp_counts$alpha_i, na.rm=T)

# Make negative binomial model with raw counts
form1 <- N.butterflies ~ -1 + # remove grand mean
  # cell ICAR random intercepts
  f(alpha_i, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # cell ICAR random effort slopes
  f(eps_i, ln.Party_Hours, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # cell ICAR random year slopes
  f(tau_i, std_Year, model="besag", graph=g50, constr=FALSE, scale.model=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # random circle intercepts
  f(kappa_k, model="iid", constr=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

# Run model
print(noquote('---running model---'))
out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3)

# Overdispersion parameter
sm1 = summary(out1)[[3]]

# Vector of grid IDs of cells with counts
cells_with_counts = unique(sp_counts$grid_id)

# get alpha summaries
alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
alph_iw <- alph_ul - alph_ll
#	par(mfrow=c(1,3))
#	hist(alph); summary(alph)
#	hist(alph_ll); summary(alph_ll)
#	hist(alph_ul); summary(alph_ul)

# get epsilon summaries
eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
eps_iw <- eps_ul - eps_ll
#	par(mfrow=c(1,3))
#	hist(eps); summary(eps); round(sum(eps<1)/length(eps), 2)
#	hist(eps_ll); summary(eps_ll)
#	hist(eps_ul); summary(eps_ul)

c1 = cor.test(eps, alph, method="spearman")

# get tau summaries
tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
tau_iw <- tau_ul - tau_ll
#	par(mfrow=c(2,2))
#	hist(tau); summary(tau); round(sum(tau>=0)/length(tau), 2)
#	hist(tau_ll); summary(tau_ll)
#	hist(tau_ul); summary(tau_ul)
#	hist(tau_iw); summary(tau_iw)

# Correlation between epsilon and alpha
c2 = cor.test(alph, tau, method="spearman")

# Visualize goodness of fit with Probability Integral Transform
#	sum(out1$cpo$failure, na.rm=T) -2 * sum(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
pit1 <- data.frame(PIT=out1$cpo$pit) %>%
  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
pit2 <- ggplot(data=pit1, aes(x=PIT)) +
  geom_histogram(col="white") +
  xlab("Probability integral transform (PIT)") +
  ylab("Count"); pit2; summary(pit1$PIT)
#	ggsave(paste0(species[i],'_pit.png'),plot=pit2,device='png',path='./plots/PIT/',width=8,height=8,units="in",dpi=300)
ggsave('Allspecies_pit_1993-2018_50km_wNAs_noTop4.png',plot=pit2,device='png',path='./plots/PIT/',width=8,height=8,units="in",dpi=300)

# collect posterior summaries into one dataframe
grid_key <- unique(sp_counts[, c("grid_id", "EcoI", "EcoII")])	
post_sum <- data.frame(grid_id=cells_with_counts,alph, alph_ll, alph_ul, alph_iw,eps, eps_ll, eps_ul, eps_iw, eps_sig=NA,tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
post_sum$eps_sig <- ifelse((post_sum$eps_ll < 0 & post_sum$eps_ul > 0),post_sum$eps_sig <- NA,post_sum$eps_sig <- post_sum$eps)
post_sum$tau_sig <- ifelse((post_sum$tau_ll < 0 & post_sum$tau_ul > 0),post_sum$tau_sig <- NA,post_sum$tau_sig <- post_sum$tau)
post_sum$EcoI <- grid_key$EcoI[match(cells_with_counts,grid_key$grid_id)]
post_sum$EcoII <- grid_key$EcoII[match(cells_with_counts,grid_key$grid_id)]
# Add medians of slope estimates to data frame
sm2 = summary(post_sum)
sm2[3,] = trimws(gsub('Median :','',sm2[3,]),which='both')

# Make cell level maps
results_cells <- merge(butterfly_grid, post_sum, by="grid_id", all=F)
write.table(results_cells@data,'Allspecies_trends_1993-2018_50km_wNAs_noTop4.txt',sep='\t',quote=F,row.names=F) # Write model output to file
res_sf <- as(results_cells, "sf")

# grid
grid2 <- merge(grid50, circles_per_cell, all=T)
grid2$number_circles[is.na(grid2$number_circles)] <- 0
grid_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=grid2, aes(fill=number_circles), col="gray40", size=0.3) +
  scale_fill_gradient("Circles\nper cell  ", low = ("white"),
					  high = ("red4"), space = "Lab", trans="log1p",
					  na.value = "grey40", guide = "colourbar",
					  breaks=c(0,2,7,20)) +
  theme_map() +
  theme(panel.grid.major=element_line(colour="transparent"))

# map tau
tau_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
					   mid = "white",
					   high = ("royalblue4"), midpoint = 0, space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p2 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_iw), col="gray40", size=0.3) +
  scale_fill_gradient2("Abund.\ntrend\ncredible\ninterval\nwidth\n(%/year)",
					   low = ("purple4"), mid = "white",
					   high = ("green4"), midpoint = median(res_sf$tau_iw), space = "Lab",
					   na.value = "grey40", guide = "colourbar",
					   breaks=c(3,6,9,12)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p3 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_sig), col="gray40", size=0.3) +
  scale_fill_gradient2("Sig.\nabund.\ntrend\n(%/year)",
					   low = muted("red4"), mid = "gray95",
					   high = muted("royalblue4"), midpoint = 0, space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# map epsilon
eps_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
  scale_fill_gradient2("Sampling\neffort", low = muted("purple4"), mid = "white",
					   high = muted("green4"), midpoint = median(res_sf$eps), space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# map alpha
alph_p1 <- ggplot() +
  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=alph), col="gray40", size=0.3) +
  scale_fill_gradient2("Relative\nabund.", low = "tan4", mid = "white",
					   high = "green4", midpoint = (max(res_sf$alph) + min(res_sf$alph)) / 2, space = "Lab",
					   na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# print cell maps
print(noquote('---drawing maps---'))
ggsave('Allspecies_effort_1993-2018_50km_wNAs_noTop4.png',plot =eps_p1,device='png',path='./',width=8,height=8,units="in",dpi=300)
ggsave('Allspecies_change_per_year_1993-2018_50km_wNAs_noTop4.png',plot =tau_p1,device='png',path='./',width=8,height=8,units="in",dpi=300)
ggsave('Allspecies_relative_abundance_1993-2018_50km_wNAs_noTop4.png',plot =alph_p1,device='png',path='./',width=8,height=8,units="in",dpi=300)

res_sf = read.table('Allspecies_trends_1993-2018_50km_wNAs_noTop4.txt',sep='\t',as.is=T,check.names=F,header=T)
par(oma=c(0,0,0,0),mar=c(0,5,0,0))
boxplot(res_sf$tau,ylab='Total abundance trend',ylim=c(-10,10),cex.lab=2,boxlwd=2,staplelwd=2,whisklwd=2,whisklty=1,col='grey40',outpch=16,frame=F,yaxt='n')
axis(2,lwd=3,at=seq(-10,10,2),cex.axis=1.5)
abline(a=0,b=0,lwd=2,col='red',lty=2)
legend('topright',legend=c('No. sites = 497','No. cells = 414'),bty='n',cex=2,ncol=1)

dev.off()
hist(res_sf$tau,main='',ylim=c(0,200),breaks=seq(-12,6,2),ylab='No. species',xlab='Median total abundance trend',cex.lab=1.5,cex.axis=1.2)
abline(v=0,lwd=2,col='red')

res_sf1 = read.table('Allspecies_trends_1993-2018_50km_wNAs.txt',sep='\t',as.is=T,check.names=F,header=T)
res_sf2 = read.table('Allspecies_trends_1993-2018_50km_wNAs_noTop4.txt',sep='\t',as.is=T,check.names=F,header=T)
quantile(res_sf1$tau)
quantile(res_sf2$tau)

res_sf1[which(res_sf1$tau<(-10)),] #grid_id 3255 trend = -10.8%/year
sp_counts[which(sp_counts$grid_id==3255),]

res_sf1[which(res_sf1$tau>5),] #grid_id 688 trend = +5.96%/year
sp_counts[which(sp_counts$grid_id==688),]


