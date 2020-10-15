
# Code used to estimate trends in butterfly abundance, per species per grid cell

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


##########################################################################################################################################################
# 50 km scale, pseudoabsences reported as NA, JunAug set

# Import butterfly data
butterfly_counts = read.table('./data/butterfly_data_gridded_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km_NoMergeJunAug_m5.shp',layer='butterfly_sites_50km_NoMergeJunAug_m5',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
butterfly_grid = readOGR(dsn='./shapefiles/butterfly_grid_50km_NoMergeJunAug_m5.shp',layer='butterfly_grid_50km_NoMergeJunAug_m5',verbose=F,stringsAsFactors=F) #shapefile with grid
grid50 <- as(butterfly_grid, "sf")
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

species = sort(unique(butterfly_counts$Species)) #397 species

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.butterflies'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
for (i in 1:length(species)){
	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]
	sp_counts = butterfly_counts[which(butterfly_counts$Species==species[i]),] 	# Subset data for species i
	sp_sites = unique(sp_counts$Site)
	sp.summary$Total.N.butterflies[i] = sum(sp_counts$N.butterflies,na.rm=T)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]

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
	out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3)

	# Overdispersion parameter
	sm1 = summary(out1)[[3]]
	sp.summary$Overdispersion[i] = sm1[1,4] # > 1 implies overdispersion relative to a poisson distribution. Justifies use of the negative binomial (formula: exp(-log(1/sm1[1,4]))

	# Vector of grid IDs of cells with counts
	cells_with_counts = unique(sp_counts$grid_id)

	# get alpha summaries
	alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
	alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
	alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
	alph_iw <- alph_ul - alph_ll

	# get epsilon summaries
	eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
	eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
	eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
	eps_iw <- eps_ul - eps_ll

	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c1 = cor.test(eps, alph, method="spearman")
		sp.summary$Epsilon_Alpha_rho[i] = c1[[4]]
		sp.summary$Epsilon_Alpha_rho_p[i] = c1[[3]]
	}
	
	# get tau summaries
	tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
	tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
	tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
	tau_iw <- tau_ul - tau_ll
	
	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c2 = cor.test(alph, tau, method="spearman")
		sp.summary$Alpha_Tau_rho[i] = c2[[4]]
		sp.summary$Alpha_Tau_rho_p[i] = c2[[3]]
	}

	# Visualize goodness of fit with Probability Integral Transform
	pit1 <- data.frame(PIT=out1$cpo$pit) %>%
	  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
	pit2 <- ggplot(data=pit1, aes(x=PIT)) +
	  geom_histogram(col="white") +
	  xlab("Probability integral transform (PIT)") +
	  ylab("Count"); pit2; summary(pit1$PIT)
	ggsave(paste0(species[i],'_pit_1993-2018_50km_noMergeJunAug_m5.png'),plot=pit2,device='png',path='./plots/PIT/',width=8,height=8,units="in",dpi=300)

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
	sp.summary$mdn_Alph[i] = as.numeric(sm2[3,2])
	sp.summary$mdn_Alph_ll[i] = as.numeric(sm2[3,3])
	sp.summary$mdn_Alph_ul[i] = as.numeric(sm2[3,4])
	sp.summary$mdn_Alph_iw[i] = as.numeric(sm2[3,5])
	sp.summary$mdn_Eps[i] = as.numeric(sm2[3,6])
	sp.summary$mdn_Eps_ll[i] = as.numeric(sm2[3,7])
	sp.summary$mdn_Eps_ul[i] = as.numeric(sm2[3,8])
	sp.summary$mdn_Eps_iw[i] = as.numeric(sm2[3,9])
	sp.summary$mdn_Tau[i] = as.numeric(sm2[3,11])
	sp.summary$mdn_Tau_ll[i] = as.numeric(sm2[3,12])
	sp.summary$mdn_Tau_ul[i] = as.numeric(sm2[3,13])
	sp.summary$mdn_Tau_iw[i] = as.numeric(sm2[3,14])

	# Make cell level maps
	results_cells <- merge(butterfly_grid, post_sum, by="grid_id", all=F)
	write.table(results_cells@data,paste0('./inla_model_output/50km_wNAs_noMerge/',species[i],'_trends_1993-201_50km_noMergeJunAug.txt'),sep='\t',quote=F,row.names=F) # Write model output to file
	res_sf <- as(results_cells, "sf")

	# map tau
	tau_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
	  geom_sf(data=res_sf, aes(fill=tau), col="gray40", size=0.3) +
	  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
						   mid = "white",
						   high = ("royalblue4"), midpoint = 0, space = "Lab",
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
	ggsave(paste0(species[i],'_effort_1993-2018_50km_noMergeJunAug.png'),plot =eps_p1,device='png',path='./plots/model_maps/50km_wNAs_noMerge/effort/',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_change_per_year_1993-2018_50km_noMergeJunAug.png'),plot =tau_p1,device='png',path='./plots/model_maps/50km_wNAs_noMerge/change_per_year/',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_relative_abundance_1993-2018_50km_noMergeJunAug.png'),plot =alph_p1,device='png',path='./plots/model_maps/50km_wNAs_noMerge/relative_abundance/',width=8,height=8,units="in",dpi=300)
	
}
write.table(sp.summary,'./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5.txt',sep='\t',quote=F,row.names=F)

trends = read.table('./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
N = nrow(trends); N
up = length(which(trends$mdn_Tau>1)); up
100*up/N
down = length(which(trends$mdn_Tau<(-1))); down
100*down/N


##########################################################################################################################################################
# 50 km scale, pseudoabsences reported as NA, Jul set

# Import butterfly data
butterfly_counts = read.table('./data/butterfly_data_gridded_50km_noMergeJul_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km_NoMergeJul_m5.shp',layer='butterfly_sites_50km_NoMergeJul_m5',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
butterfly_grid = readOGR(dsn='./shapefiles/butterfly_grid_50km_NoMergeJul_m5.shp',layer='butterfly_grid_50km_NoMergeJul_m5',verbose=F,stringsAsFactors=F) #shapefile with grid
grid50 <- as(butterfly_grid, "sf")
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors
eco1 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l1/NA_CEC_Eco_Level1.shp',layer='NA_CEC_Eco_Level1',verbose=F,stringsAsFactors=F),CRS(proj1)) #ecoregion shapefiles - used for mapping later
eco_sf = as(eco1, 'sf')

species = sort(unique(butterfly_counts$Species)) #397 species

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.butterflies'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
for (i in 1:length(species)){
	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]
	sp_counts = butterfly_counts[which(butterfly_counts$Species==species[i]),] 	# Subset data for species i
	sp_sites = unique(sp_counts$Site)
	sp.summary$Total.N.butterflies[i] = sum(sp_counts$N.butterflies,na.rm=T)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]

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
	out1 <- inla(form1, family="nbinomial", data=sp_counts, control.compute=list(cpo=T, config=T), control.inla=list(strategy="adaptive", int.strategy="auto"), num.threads=3)

	# Overdispersion parameter
	sm1 = summary(out1)[[3]]
	sp.summary$Overdispersion[i] = sm1[1,4] # > 1 implies overdispersion relative to a poisson distribution. Justifies use of the negative binomial (formula: exp(-log(1/sm1[1,4]))

	# Vector of grid IDs of cells with counts
	cells_with_counts = unique(sp_counts$grid_id)

	# get alpha summaries
	alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
	alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
	alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
	alph_iw <- alph_ul - alph_ll

	# get epsilon summaries
	eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
	eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
	eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
	eps_iw <- eps_ul - eps_ll

	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c1 = cor.test(eps, alph, method="spearman")
		sp.summary$Epsilon_Alpha_rho[i] = c1[[4]]
		sp.summary$Epsilon_Alpha_rho_p[i] = c1[[3]]
	}
	
	# get tau summaries
	tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts]) - 1) * 100
	tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts]) - 1) * 100
	tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts]) - 1) * 100
	tau_iw <- tau_ul - tau_ll
	
	# Correlation between epsilon and alpha
	if (length(cells_with_counts)<2){
	} else {
		c2 = cor.test(alph, tau, method="spearman")
		sp.summary$Alpha_Tau_rho[i] = c2[[4]]
		sp.summary$Alpha_Tau_rho_p[i] = c2[[3]]
	}

	# Visualize goodness of fit with Probability Integral Transform
	pit1 <- data.frame(PIT=out1$cpo$pit) %>%
	  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
	pit2 <- ggplot(data=pit1, aes(x=PIT)) +
	  geom_histogram(col="white") +
	  xlab("Probability integral transform (PIT)") +
	  ylab("Count"); pit2; summary(pit1$PIT)
	ggsave(paste0(species[i],'_pit_1993-2018_50km_noMergeJul_m5.png'),plot=pit2,device='png',path='./plots/PIT/',width=8,height=8,units="in",dpi=300)

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
	sp.summary$mdn_Alph[i] = as.numeric(sm2[3,2])
	sp.summary$mdn_Alph_ll[i] = as.numeric(sm2[3,3])
	sp.summary$mdn_Alph_ul[i] = as.numeric(sm2[3,4])
	sp.summary$mdn_Alph_iw[i] = as.numeric(sm2[3,5])
	sp.summary$mdn_Eps[i] = as.numeric(sm2[3,6])
	sp.summary$mdn_Eps_ll[i] = as.numeric(sm2[3,7])
	sp.summary$mdn_Eps_ul[i] = as.numeric(sm2[3,8])
	sp.summary$mdn_Eps_iw[i] = as.numeric(sm2[3,9])
	sp.summary$mdn_Tau[i] = as.numeric(sm2[3,11])
	sp.summary$mdn_Tau_ll[i] = as.numeric(sm2[3,12])
	sp.summary$mdn_Tau_ul[i] = as.numeric(sm2[3,13])
	sp.summary$mdn_Tau_iw[i] = as.numeric(sm2[3,14])

	# Make cell level maps
	results_cells <- merge(butterfly_grid, post_sum, by="grid_id", all=F)
	write.table(results_cells@data,paste0('./inla_model_output/50km_wNAs_noMerge/',species[i],'_trends_1993-201_50km_noMergeJul.txt'),sep='\t',quote=F,row.names=F) # Write model output to file
	res_sf <- as(results_cells, "sf")

	# map tau
	tau_p1 <- ggplot() +
	  geom_sf(data=eco_sf, fill="gray40", col="gray40") +
	  geom_sf(data=res_sf, aes(fill=tau), col="gray40", size=0.3) +
	  scale_fill_gradient2("Abund.\ntrend\n(%/year)", low = ("red4"),
						   mid = "white",
						   high = ("royalblue4"), midpoint = 0, space = "Lab",
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
	ggsave(paste0(species[i],'_effort_1993-2018_50km_noMergeJul.png'),plot =eps_p1,device='png',path='./plots/model_maps/50km_wNAs_noMerge/effort/',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_change_per_year_1993-2018_50km_noMergeJul.png'),plot =tau_p1,device='png',path='./plots/model_maps/50km_wNAs_noMerge/change_per_year/',width=8,height=8,units="in",dpi=300)
	ggsave(paste0(species[i],'_relative_abundance_1993-2018_50km_noMergeJul.png'),plot =alph_p1,device='png',path='./plots/model_maps/50km_wNAs_noMerge/relative_abundance/',width=8,height=8,units="in",dpi=300)
	
}
write.table(sp.summary,'./butterfly_trends_summary_1993-2018_50km_noMergeJul_m5.txt',sep='\t',quote=F,row.names=F)


