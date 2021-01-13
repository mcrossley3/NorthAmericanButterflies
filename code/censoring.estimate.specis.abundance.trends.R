

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

##########################################################################################################################################################
# 50 km scale, pseudoabsences reported as NA, JunAug set


#####
# Left-censored (trim first 5 years before estimating trends)

# Import butterfly data
butterfly_counts = read.table('./data/butterfly_data_gridded_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors

butterfly_counts2 = butterfly_counts[which(butterfly_counts$Year>1997),]
species = sort(unique(butterfly_counts2$Species))
out = data.frame('Species'=NA,'grid_id'=NA,'Nyears'=NA,'Nobs'=NA)
ox = 1
for (s in 1:length(species)){
	dat = butterfly_counts2[which(butterfly_counts2$Species==species[s]),]
	gids = unique(dat$grid_id)
	for (g in 1:length(gids)){
		dat2 = dat[which(dat$grid_id==gids[g]),]
		out[ox,1] = species[s]
		out[ox,2] = gids[g]
		out[ox,3] = length(unique(dat2$Year))
		out[ox,4] = nrow(dat2)
		ox = ox + 1
	}
}; dim(out)
out2 = out[which(out$Nyears>=10 & out$Nobs>=5),]; dim(out2)
keep.pos = c()
for (i in 1:nrow(out2)){
	print(noquote(i))
	keep.pos = c(keep.pos,which(butterfly_counts2$Species==out2$Species[i] & butterfly_counts2$grid_id==out2$grid_id[i]))
}
butterfly_counts3 = butterfly_counts2[keep.pos,]
species = sort(unique(butterfly_counts3$Species)) #down to 373 species (from 456)

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.butterflies'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
for (i in 1:length(species)){
	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]
	sp_counts = butterfly_counts3[which(butterfly_counts3$Species==species[i]),] 	# Subset data for species i
	sp.summary$Total.N.butterflies[i] = sum(sp_counts$N.butterflies,na.rm=T)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]
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
	write.table(results_cells@data,paste0('./inla_model_output/50km_wNAs_noMerge/leftcensored/',species[i],'_trends_1993-201_50km_noMergeJunAug_leftcensored.txt'),sep='\t',quote=F,row.names=F) # Write model output to file	
}
write.table(sp.summary,'./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5_leftcensored.txt',sep='\t',quote=F,row.names=F)


#####
# Right-censoring

# Import butterfly data
butterfly_counts = read.table('./data/butterfly_data_gridded_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
g50 <- inla.read.graph('./shapefiles/nb50.graph') #graph network of grid neighbors

butterfly_counts2 = butterfly_counts[which(butterfly_counts$Year<2014),]
species = sort(unique(butterfly_counts2$Species))
out = data.frame('Species'=NA,'grid_id'=NA,'Nyears'=NA,'Nobs'=NA)
ox = 1
for (s in 1:length(species)){
	dat = butterfly_counts2[which(butterfly_counts2$Species==species[s]),]
	gids = unique(dat$grid_id)
	for (g in 1:length(gids)){
		dat2 = dat[which(dat$grid_id==gids[g]),]
		out[ox,1] = species[s]
		out[ox,2] = gids[g]
		out[ox,3] = length(unique(dat2$Year))
		out[ox,4] = nrow(dat2)
		ox = ox + 1
	}
}; dim(out)
out2 = out[which(out$Nyears>=10 & out$Nobs>=5),]; dim(out2)
keep.pos = c()
for (i in 1:nrow(out2)){
	print(noquote(i))
	keep.pos = c(keep.pos,which(butterfly_counts2$Species==out2$Species[i] & butterfly_counts2$grid_id==out2$grid_id[i]))
}
butterfly_counts3 = butterfly_counts2[keep.pos,]
species = sort(unique(butterfly_counts3$Species)) #down to 380 species (from 456)

f1 = rep(NA,length(species))
sp.summary = data.frame('Species'=f1,'Overdispersion'=f1,'Epsilon_Alpha_rho'=f1,'Epsilon_Alpha_rho_p'=f1,'Alpha_Tau_rho'=f1,'Alpha_Tau_rho_p'=f1,'mdn_Alph'=f1,
	'mdn_Alph_ll'=f1,'mdn_Alph_ul'=f1,'mdn_Alph_iw'=f1,'mdn_Eps'=f1,'mdn_Eps_ll'=f1,'mdn_Eps_ul'=f1,'mdn_Eps_iw'=f1,'mdn_Tau'=f1,'mdn_Tau_ll'=f1,'mdn_Tau_ul'=f1,
	'mdn_Tau_iw'=f1,'Total.N.butterflies'=f1,'N.sites'=f1,'N.years'=f1,'First.year'=f1,'Last.year'=f1)
for (i in 1:length(species)){
	print(noquote(paste0(i,'. ',species[i])))
	sp.summary$Species[i] = species[i]
	sp_counts = butterfly_counts3[which(butterfly_counts3$Species==species[i]),] 	# Subset data for species i
	sp.summary$Total.N.butterflies[i] = sum(sp_counts$N.butterflies,na.rm=T)
	sp.summary$N.sites[i] = length(unique(sp_counts$Site))
	sp.years = sort(as.numeric(unique(sp_counts$Year)))
	sp.summary$N.years[i] = length(sp.years)
	sp.summary$First.year[i] = sp.years[1]
	sp.summary$Last.year[i] = rev(sp.years)[1]
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
	write.table(results_cells@data,paste0('./inla_model_output/50km_wNAs_noMerge/rightcensored/',species[i],'_trends_1993-201_50km_noMergeJunAug_rightcensored.txt'),sep='\t',quote=F,row.names=F) # Write model output to file
}
write.table(sp.summary,'./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5_rightcensored.txt',sep='\t',quote=F,row.names=F)


######
# Compare trends with and without left or right censoring

trends1 = read.table('./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T)
trends2 = read.table('./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5_leftcensored.txt',sep='\t',as.is=T,check.names=F,header=T)
trends3 = read.table('./butterfly_trends_summary_1993-2018_50km_noMergeJunAug_m5_rightcensored.txt',sep='\t',as.is=T,check.names=F,header=T)

merge1 = merge(trends1,trends2,by='Species')
merge2 = merge(merge1,trends3,by='Species')

par(mfrow=c(1,3))
plot(merge2$mdn_Tau.x,merge2$mdn_Tau.y,xlab='Median abund. trend uncensored',ylab='Median abund. trend left-censored'); abline(h=0); abline(v=0)
plot(merge2$mdn_Tau.x,merge2$mdn_Tau,xlab='Median abund. trend uncensored',ylab='Median abund. trend right-censored'); abline(h=0); abline(v=0)
plot(merge2$mdn_Tau.y,merge2$mdn_Tau,xlab='Median abund. trend left-censored',ylab='Median abund. trend right-censored'); abline(h=0); abline(v=0)


