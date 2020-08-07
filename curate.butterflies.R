
setwd('<your/working/directory>')


#################################################################################
# Check & Merge data (merge sub-species & identical sites)

library(lubridate)

# Define function that ensures an R-readable date format
convert.date = function(mdy){
	if (is.na(mdy)){
		return(NA)
	} else {
		mdy.split = strsplit(mdy,'/')[[1]]
		YYYY = mdy.split[3]
		if (as.numeric(YYYY)<70){YYYY=paste0('20',YYYY)}else{YYYY=paste0('19',YYYY)}
		MM = apply(array(mdy.split[1]),1,function(x){if(nchar(x)<2){paste0(0,x)}else{x}})
		DD = apply(array(mdy.split[2]),1,function(x){if(nchar(x)<2){paste0(0,x)}else{x}})
		return(paste(YYYY,MM,DD,sep='-'))
	}
}


# Import data
dat = read.csv('./data/NABA_Count_Data_through_2018.csv',as.is=T,check.names=F,header=T)
colnames(dat)[1] = 'CountName' #overwrite original column name, which contains invisible characters
lonlat = paste(dat$Lng,dat$Lat,sep='_')
dat$SiteLonLat = paste(dat$CountName,lonlat,sep='_')
str(dat)
sum(dat$NumSeen,na.rm=T) #11,950,146 butterflies
length(unique(dat$CountName)) #1,221 sites


# Remove genus-level records
species.split = apply(array(dat$ScientificName),1,function(x){strsplit(x,' ')[[1]]})
sslen = unlist(lapply(species.split,function(x){length(x)}))
keep1 = which(sslen>1)
dat = dat[keep1,]; str(dat) #388,675 records with species-level records (lost 15,737)


# Get R-readable dates
dat$Date2 = apply(array(dat$Date),1,function(x){convert.date(x)}); str(dat$Date2) #R-readable date format
dat$Julian = yday(dat$Date2); str(dat$Julian)
hist(dat$Julian) #most observations occur in June-July


# Merge subspecies-level records to species-level
species.split = species.split[keep1]
sslen = sslen[keep1]
dat$GenusSpecies = unlist(lapply(species.split,function(x){paste(x[1],x[2],sep=' ')}))
species = sort(unique(dat$GenusSpecies))
dat$GenusSpecies[which(dat$GenusSpecies=='Colias philodiceXeurytheme')] = 'Colias philodice'
dat$GenusSpecies[which(dat$GenusSpecies=='Colias philodice/eurytheme')] = 'Colias philodice'
dat$GenusSpecies[which(dat$GenusSpecies=='Limenitis archippusXweidemeyerii')] = 'Limenitis archippus'
dat$GenusSpecies[which(dat$GenusSpecies=='Limenitis arthemisXarchippus')] = 'Limenitis arthemis'
dat$GenusSpecies[which(dat$GenusSpecies=='Limenitis arthemisXweidemeyerii')] = 'Limenitis arthemis'
dat$GenusSpecies[which(dat$GenusSpecies=='Limenitis lorquiniXarthemis')] = 'Limenitis lorquini'
dat$GenusSpecies[which(dat$GenusSpecies=='Boloria montinus')] = 'Boloria chariclea' #merge synonym
dat$GenusSpecies[which(dat$GenusSpecies=='Hemiargus ammon')] = 'Cyclargus ammon' #might be synonymous with C. thomasi, but both records appear in the same site*year, so we'll assume observers were referring to two different species
dat$GenusSpecies[which(dat$GenusSpecies=='Hemiargus thomasi')] = 'Cyclargus thomasi'
uspecies = sort(unique(dat$GenusSpecies)) #644 unique species
species.merge = sort(unique(c(dat$GenusSpecies[which(sslen>2)],'Colias philodice','Limenitis archippus','Limenitis lorquini','Satyrium titus','Callophrys augustinus')))
site.remove.pos = unlist(apply(array(species.merge),1,function(x){which(dat$GenusSpecies==x)}))
dat2 = dat[-site.remove.pos,] #331,274 records left; 57,401 records will be merged
add.dat = c()
for (m in 1:length(species.merge)){ #WARNING: this loop takes several minutes
	print(paste(species.merge[m],m,'of',length(species.merge)))
	mdat = dat[which(dat$GenusSpecies==species.merge[m]),]
	msites = unique(mdat$SiteLonLat)
	for (s in 1:length(msites)){
		smdat = mdat[which(mdat$SiteLonLat==msites[s]),]
		sdates = unique(smdat$Date2)
		for (d in 1:length(sdates)){
			dsmdat = smdat[which(smdat$Date2==sdates[d]),]
			add.line = dsmdat[1,]
			add.line$NumSeen = sum(dsmdat$NumSeen,na.rm=T)
			add.line$NumObservers = sum(dsmdat$NumObservers,na.rm=T)
			add.line$Party_Hours = sum(dsmdat$Party_Hours,na.rm=T)
			add.dat = data.frame(rbind(add.dat,add.line))
		}
	}
}
dat2 = data.frame(rbind(dat2,add.dat)); str(dat2) #dataframe contains species-level records merged by site*date; 387,916 records left (54,401 records merged into 759).


# Merge sites that share coordinates but are inconsistently named 
con = file('./data/sites_to_be_merged2.txt','r') #manipulated in excel to identify sites that should be merged. 90 merges are proposed.
site.lines = readLines(con); close(con)
site.merge.key = t(apply(array(site.lines),1,function(x){strsplit(x,'\t')[[1]]})); colnames(site.merge.key) = site.merge.key[1,]; site.merge.key = site.merge.key[-1,]
key.sitelonlat = paste(site.merge.key[,1],site.merge.key[,3],site.merge.key[,2],sep='_')
site.remove.pos = unlist(apply(array(key.sitelonlat[which(site.merge.key[,4]!='NA')]),1,function(x){which(dat2$SiteLonLat==x)}))
dat3 = dat2[-site.remove.pos,]; str(dat3) #326,665 records retained (62,010 records will be merged)
mergers = unique(site.merge.key[,4]); mergers = mergers[which(mergers!='NA')]
add.dat = c()
for (m in 1:length(mergers)){
	print(m)
	sites2merge = key.sitelonlat[which(site.merge.key[,4]==mergers[m])]
	schars = nchar(sites2merge)
	merge.name = sites2merge[order(schars)][1] #use shorter name as the name of merged records
	mdat = c()
	for (i in 1:length(sites2merge)){
		mdat = data.frame(rbind(mdat,dat2[which(dat2$Site==sites2merge[i]),]))
	}
	umdates = unique(mdat$Date2)
	for (d in 1:length(umdates)){
		dmdat = mdat[which(mdat$Date2==umdates[d]),]
		umspecies = unique(dmdat$GenusSpecies)
		for (s in 1:length(umspecies)){
			sdmdat = dmdat[which(dmdat$GenusSpecies==umspecies[s]),]
			if (nrow(sdmdat)>1){print(umspecies[s])}else{}
			sdmdat$SiteLonLat = paste(merge.name,'merged',sep='_')
			add.dat = data.frame(rbind(add.dat,sdmdat))
		}
	}
}
dat3 = data.frame(rbind(dat3,add.dat)); str(dat3) # 387,916 records left; 62,010 records merged into 61,251 


# Compile into a trim table
butterflies = data.frame(
	'Site'=dat3$SiteLonLat,
	'Year'=apply(array(dat3$Date2),1,function(x){strsplit(x,'-')[[1]][1]}),
	'Species'=dat3$GenusSpecies,
	'N.butterflies'=dat3$NumSeen,
	'NumObservers'=dat3$NumObservers,
	'Party_Hours'=dat3$Party_Hours,
	'Latitude'=dat3$Lat,
	'Longitude'=dat3$Lng,
	'Julian.date'=dat3$Julian,
	stringsAsFactors=F)


# Fix typoes in species names
butterflies$Species[which(butterflies$Species=='Colias scudderi')] = 'Colias scudderii'
butterflies$Species[which(butterflies$Species=='Neophasia terlootii')] = 'Neophasia terlooii'


#########
# Merge butterfly counts at site*years that had more than one count event (e.g. counts on two different dates)

# Check species duplicates (occurs when site was sampled multiple times within a year
sites = unique(butterflies$Site)
site.year.check = data.frame('Site'=NA,'Year'=NA,'Check'=NA)
sx = 1
for (s in 1:length(sites)){
	print(s)
	dat = butterflies[which(butterflies$Site==sites[s]),]
	years = sort(as.numeric(unique(dat$Year)))
	for (y in 1:length(years)){
		dat2 = dat[which(dat$Year==years[y]),]
		species = unique(dat2$Species)
		check.counter = 0
		for (i in 1:length(species)){
			check = which(dat2$Species==species[i])
			if (length(check)>1){
				check.counter = check.counter + 1
			} else {
			}
		}
		site.year.check[sx,1] = sites[s]
		site.year.check[sx,2] = years[y]
		site.year.check[sx,3] = check.counter
		sx = sx + 1
	}
}


# Generate merged counts
merge.site.years = site.year.check[which(site.year.check[,3]>0),1:2]
add.back = data.frame()
for (i in 1:nrow(merge.site.years)){
	print(i)
	dat = butterflies[which(butterflies$Site==merge.site.years[i,1] & butterflies$Year==merge.site.years[i,2]),]
	juliandate = paste0(unique(dat$Julian.date),collapse=';')
	latitude = unique(dat$Latitude)
	if (length(latitude)>1){
		latitude = latitude[2]
	} else {
	}
	longitude = unique(dat$Longitude)
	if (length(longitude)>1){
		longitude = longitude[2]
	} else {
	}
	site = merge.site.years[i,1]
	year = merge.site.years[i,2]
	numobservers = sum(unique(dat$NumObservers))
	partyhours = sum(unique(dat$Party_Hours))
	species = unique(dat$Species)
	for (j in 1:length(species)){
		spos = which(dat$Species==species[j])
		dat2 = dat[spos,]
		species1 = species[j]
		nbutterflies = sum(dat2$N.butterflies)
		add1 = data.frame('Site'=site,'Year'=year,'Species'=species1,'N.butterflies'=nbutterflies,'NumObservers'=numobservers,'Party_Hours'=partyhours,'Latitude'=latitude,
			'Longitude'=longitude,'Julian.date'=juliandate,stringsAsFactors=F)
		add.back = data.frame(rbind(add.back,add1),stringsAsFactors=F)
	}
}


# Replace duplicated species records with merged counts
remove.pos = unlist(apply(merge.site.years,1,function(x){which(butterflies$Site==x[1] & butterflies$Year==x[2])}))
butterflies2 = butterflies[-remove.pos,]
butterflies2 = data.frame(rbind(butterflies2,add.back),stringsAsFactors=F)

# Double check
sites = unique(butterflies2$Site)
for (s in 1:length(sites)){
	print(s)
	dat = butterflies2[which(butterflies2$Site==sites[s]),]
	years = sort(as.numeric(unique(dat$Year)))
	for (y in 1:length(years)){
		dat2 = dat[which(dat$Year==years[y]),]
		species = unique(dat2$Species)
		for (i in 1:length(species)){
			check = which(dat2$Species==species[i])
			if (length(check)>1){
				print(noquote('............................flag'))
			} else {
			}
		}
	}
}

str(butterflies2)
length(unique(butterflies2$Site))
length(unique(butterflies2$Species))
sort(unique(butterflies2$Year))

# Write data to file
write.table(butterflies2,'./data/butterfly_data_merged.txt',sep='\t',row.names=F)


################################################################################
# Explore data

library(rgdal)
library(sp)

# Import tidy merged dataset
butterflies = read.table('./data/butterfly_data_merged.txt',sep='\t',as.is=T,check.names=F,header=T)
sum(butterflies$N.butterflies,na.rm=T) #11,949,896

# Check that the number of species never exceeds number of butterflies
which(butterflies$N.species > butterflies$N.butterflies)

sites = sort(unique(butterflies$Site))
years = sort(unique(butterflies$Year))

# Site summary
out = data.frame('Site'=NA,'N.species'=-999,'Total.n.butterflies'=-999,'Total.party.hrs'=-999,'Lon'=-999,'Lat'=-999,'First.year'=NA,'Last.year'=NA,'N.years'=-999)
for (s in 1:length(sites)){
	s.dat = butterflies[which(butterflies$Site==sites[s]),]
	n.species = length(unique(s.dat$Species))
	sum.count = sum(s.dat$N.butterflies,na.rm=T)
	sum.partyhours = sum(s.dat$Party_Hours,na.rm=T)
	lon = unique(s.dat$Longitude)
	lat = unique(s.dat$Latitude)
	if (length(lon)>1 | length(lat)>1){ #deal with case where sites were merged, but noisy differences in coordinates remain
		lonlat = strsplit(sites[s],'_')[[1]][2:3]
		lon = lonlat[1]
		lat = lonlat[2]
	} else {
	}
	syears = sort(unique(s.dat$Year))
	out[s,1] = sites[s]
	out[s,2] = n.species
	out[s,3] = sum.count
	out[s,4] = sum.partyhours
	out[s,5] = lon
	out[s,6] = lat
	out[s,7] = syears[1]
	out[s,8] = rev(syears)[1]
	out[s,9] = length(syears)
}
write.table(out,'site_summary.txt',sep='\t',row.names=F)


# Write site summary to shapefile
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
shp = spTransform(SpatialPointsDataFrame(coords=cbind(as.numeric(out$Lon),as.numeric(out$Lat)),data=out,proj4string=CRS(wgs84)),CRS(proj1))
# Obtained ecoregions shapefile from https://www.epa.gov/eco-research/ecoregions-north-america
eco2 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l2/NA_CEC_Eco_Level2.shp',layer='NA_CEC_Eco_Level2',verbose=F,stringsAsFactors=F),CRS(proj1))
over2 = over(shp,eco2)
shp@data$EcoregionI = over2$NA_L1CODE
shp@data$EcoregionII = over2$NA_L2CODE
writeOGR(shp,dsn='./shapefiles/butterfly_sites.shp',layer='butterfly_sites',driver='ESRI Shapefile',overwrite=T)
over2$Site = shp@data$Site
write.table(over2,'Ecoregion_key.txt',sep='\t',row.names=F)


# Site*Year summary
out = data.frame('Site'=NA,'Year'=NA,'N.species'=-999,'Total.n.butterflies'=-999,'Total.party.hrs'=-999,'Lon'=-999,'Lat'=-999,'Julian.date'=-999,'Max.n.butterflies'=-999)
ox = 1
for (s in 1:length(sites)){
	s.dat = butterflies[which(butterflies$Site==sites[s]),]
	syears = sort(unique(s.dat$Year))
	for (y in 1:length(syears)){
		y.dat = s.dat[which(s.dat$Year==syears[y]),]
		n.species = length(unique(y.dat$Species))
		sum.count = sum(y.dat$N.butterflies,na.rm=T)
		sum.partyhours = sum(y.dat$Party_Hours,na.rm=T)
		lon = unique(y.dat$Longitude)
		lat = unique(y.dat$Latitude)
		if (length(lon)>1 | length(lat)>1){ #deal with case where sites were merged, but noisy differences in coordinates remain
			lonlat = strsplit(sites[s],'_')[[1]][2:3]
			lon = lonlat[1]
			lat = lonlat[2]
		} else {
		}
		julian1 = unique(y.dat$Julian.date)
		out[ox,1] = sites[s]
		out[ox,2] = syears[y]
		out[ox,3] = n.species
		out[ox,4] = sum.count
		out[ox,5] = sum.partyhours
		out[ox,6] = lon
		out[ox,7] = lat
		out[ox,8] = julian1
		out[ox,9] = max(y.dat$N.butterflies,na.rm=T)
		ox = ox + 1
	}
}
write.table(out,'site_year_summary.txt',sep='\t',row.names=F)


# Examine coverage among sites & years
out = read.table('site_year_summary.txt',sep='\t',as.is=T,check.names=F,header=T); str(out)
length(unique(out$Site))
out2 = out[which(out$Year>1992),]
sites2 = unique(out2$Site)
length(sites2) #1,049 of 1,173 sites remain
hist(as.numeric(out$Julian.date))


#################################################################################################
# Prepare data for modeling

library(rgdal)
library(scales)
library(sp)
library(rgeos)
library(maptools)
library(ggplot2)
library(spdep)
library(raster)


# Import butterfly shapefile & transform CRS to match spatial grids
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
shp = spTransform(readOGR(dsn='./shapefiles/butterfly_sites.shp',layer='butterfly_sites',verbose=F,stringsAsFactors=F),CRS(proj1))
dp = which(shp@data$N_years>9)
shp2 = shp[dp,] #498 sites left
length(unique(shp@data$Site)); length(unique(shp2@data$Site))
writeOGR(shp2,dsn='./shapefiles/butterfly_sites_10yr.shp',layer='butterfly_sites_10yr',driver='ESRI Shapefile',overwrite=T)


# Import butterfly count data & prune sites with too few years
butterflies = read.table('./data/butterfly_data_merged.txt',sep='\t',as.is=T,check.names=F,header=T)
deep_series = shp@data$Site[dp]
butterflies2 = butterflies[which(!is.na(match(butterflies$Site,deep_series))),] #318,219 out of 387,916 records remain


# Add grid_id & ecoregions to site shapefile
eco2 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l2/NA_CEC_Eco_Level2.shp',layer='NA_CEC_Eco_Level2',verbose=F,stringsAsFactors=F),CRS(proj1))
over2 = over(shp2,eco2)
shp2@data$EcoI = over2$NA_L1CODE
shp2@data$EcoII = over2$NA_L2CODE
shp2@data$EcoregionI = over2$NA_L1NAME
shp2@data$EcoregionII = over2$NA_L2NAME
# 3 sites in ocean not assigned an ecoregion - fix that by expanding a buffer around them until they touch land
shp3 = shp2[which(is.na(shp2@data$EcoI)),]
buff1 = buffer(shp3[1,],width=50000)
buff2 = buffer(shp3[2,],width=50000)
buff3 = buffer(shp3[3,],width=50000)
ov1 = over(buff1,eco2)
ov2 = over(buff2,eco2)
ov3 = over(buff3,eco2)
ovs = rbind(ov1,ov2,ov3)
shp2@data$EcoI[match(shp3@data$Site,shp2@data$Site)] = ovs$NA_L1CODE
shp2@data$EcoII[match(shp3@data$Site,shp2@data$Site)] = ovs$NA_L2CODE
shp2@data$EcoregionI[match(shp3@data$Site,shp2@data$Site)] = ovs$NA_L1NAME
shp2@data$EcoregionII[match(shp3@data$Site,shp2@data$Site)] = ovs$NA_L2NAME
na_cells2 = readOGR(na_cells2,dsn='./shapefiles/butterfly_grid_100km.shp',layer='butterfly_grid_100km',driver='ESRI Shapefile',overwrite=T)
over3 = over(shp2,na_cells2)
shp2$grid_id = over3$grid_id #grid id added to site data shapefile
writeOGR(shp2,dsn='./shapefiles/butterfly_sites_10yr_wEco.shp',layer='butterfly_sites_10yr_wEco',driver='ESRI Shapefile',overwrite=T)


# Add grid_id & ecoregions to butterfly count data
over3 = over(shp2,na_cells2)
over3$Site = shp2@data$Site
over3$EcoI = shp2@data$EcoI
over3$EcoII = shp2@data$EcoII
over3$EcoregionI = shp2@data$EcoregionI
over3$EcoregionII = shp2@data$EcoregionII
sites = unique(butterflies2$Site)
butterflies2$grid_id = NA
butterflies2$EcoregionII = butterflies2$EcoregionI = butterflies2$EcoII = butterflies2$EcoI = NA
for (i in 1:length(sites)){
	print(i)
	pos1 = which(butterflies2$Site==sites[i])
	pos2 = which(over3$Site==sites[i])
	butterflies2$grid_id[pos1] = over3$grid_id[pos2]
	butterflies2$EcoI[pos1] = over3$EcoI[pos2]
	butterflies2$EcoII[pos1] = over3$EcoII[pos2]
	butterflies2$EcoregionI[pos1] = over3$EcoregionI[pos2]
	butterflies2$EcoregionII[pos1] = over3$EcoregionII[pos2]
}
summary(butterflies2)
write.table(butterflies2,'./data/butterfly_data_gridded_100km.txt',sep='\t',row.names=F)


##############################################################
# Create 50 km grid

butterfly_counts = read.table('./data/butterfly_data_gridded_100km.txt',sep='\t',as.is=T,check.names=F,header=T)
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_10yr.shp',layer='butterfly_sites_10yr',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
eco2 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l2/NA_CEC_Eco_Level2.shp',layer='NA_CEC_Eco_Level2',verbose=F,stringsAsFactors=F),CRS(proj1))
grid50 = readOGR(dsn='./shapefiles/50km_grid_clip.shp',layer='50km_grid_clip',verbose=F,stringsAsFactors=F)
grid50@data$grid_id = seq(1,length(grid50),1)
centroids = gCentroid(grid50,byid=T)
grid50@data$centroid_x = centroids@coords[,1]
grid50@data$centroid_y = centroids@coords[,2]
grid50@data = grid50@data[,-1]
writeOGR(grid50,dsn='./shapefiles/butterfly_grid_50km.shp',layer='butterfly_grid_50km',driver='ESRI Shapefile',overwrite=T)
bover = over(butterfly_circles,grid50)
butterfly_circles$grid_id = bover$grid_id
writeOGR(butterfly_circles,dsn='./shapefiles/butterfly_sites_50km.shp',layer='butterfly_sites_50km',driver='ESRI Shapefile',overwrite=T)
gids = sort(unique(butterfly_circles@data$grid_id))
butterfly_counts2 = butterfly_counts #replace 100km grid ids with 50 km grid ids
butterfly_counts2$grid_id = NA
for (g in 1:length(gids)){
	site1 = butterfly_circles@data$Site[which(butterfly_circles@data$grid_id==gids[g])]
	for (i in 1:length(site1)){
		row.pos = which(butterfly_counts2$Site==site1[i])
		butterfly_counts2$grid_id[row.pos] = gids[g]	
	}
}
unique(butterfly_counts2$Species)
write.table(butterfly_counts2,'./data/butterfly_data_gridded_50km.txt',sep='\t',row.names=F)
nb50 <- poly2nb(grid50, row.names=grid50$grid_id); nb50
is.symmetric.nb(nb50, verbose = FALSE, force = TRUE)
nb2INLA("nb50.graph", nb50)


#######################################################################
# Prune years<1993 & add implicit zeroes to data

butterflies100 = read.table('./data/butterfly_data_gridded_100km.txt',sep='\t',as.is=T,check.names=F,header=T); str(butterflies100)
butterflies50 = read.table('./data/butterfly_data_gridded_50km.txt',sep='\t',as.is=T,check.names=F,header=T); str(butterflies50) #310,329 obs

butterflies100.2 = butterflies100[which(butterflies100$Year>1992),]; dim(butterflies100.2)
butterflies50.2 = butterflies50[which(butterflies50$Year>1992),]; dim(butterflies50.2) #287,707 obs left

site.years = unique(butterflies100.2[,1:2])
species = unique(butterflies100.2$Species)
con = file('./data/add_100km_1993-2018_w0s.txt','w')
writeLines(paste0(c('Site','Year','Species','N.butterflies','NumObservers','Party_Hours','Latitude','Longitude','Julian.date','grid_id','EcoI','EcoII','EcoregionI','EcoregionII'),collapse='\t'),con)
for (i in 1:nrow(site.years)){
	print(i)
	s1 = site.years[i,1]
	y1 = site.years[i,2]
	sydata = butterflies100.2[which(butterflies100.2[,1]==s1 & butterflies100.2[,2]==y1),]
	sydata$N.butterflies[which(is.na(sydata$N.butterflies))] = 0
	c12 = sydata[1,1:2]
	c515 = sydata[1,5:14]
	missing.species = species[which(is.na(match(species,unique(sydata$Species))))]
	nm = length(missing.species)
	add.data = data.frame('Site'=rep(c12[1,1],nm),'Year'=rep(c12[1,2],nm),'Species'=missing.species,'N.butterflies'=rep(0,nm),'NumObservers'=rep(c515[1,1],nm),
		'Party_Hours'=rep(c515[1,2],nm),'Latitude'=rep(c515[1,3],nm),'Longitude'=rep(c515[1,4],nm),'Julian.date'=rep(c515[1,5],nm),
		'grid_id'=rep(c515[1,6],nm),'EcoI'=rep(c515[1,7],nm),'EcoII'=rep(c515[1,8],nm),'EcoregionI'=rep(c515[1,9],nm),'EcoregionII'=rep(c515[1,10],nm),stringsAsFactors=F)
	for (j in 1:nrow(add.data)){
		writeLines(paste0(add.data[j,],collapse='\t'),con)
	}
}
close(con)
con = file('./data/add_100km_1993-2018_w0s.txt','r')
add.lines = readLines(con)
close(con)
add.data = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(add.data) = strsplit(add.lines[1],'\t')[[1]]; str(add.data)
butterflies100.2 = data.frame(rbind(butterflies100.2,add.data),stringsAsFactors=F)
butterflies100.2$N.butterflies = as.numeric(butterflies100.2$N.butterflies)
butterflies100.2$NumObservers = as.numeric(butterflies100.2$NumObservers)
butterflies100.2$Party_Hours = as.numeric(butterflies100.2$Party_Hours)
butterflies100.2$Latitude = as.numeric(butterflies100.2$Latitude)
butterflies100.2$Longitude = as.numeric(butterflies100.2$Longitude)
write.table(butterflies100.2,'./data/butterfly_data_gridded_100km_1993-2018_w0s.txt',sep='\t',row.names=F)


# 50km
site.years = unique(butterflies50.2[,1:2])
species = unique(butterflies50.2$Species)
con = file('./data/add_50km_1993-2018_w0s.txt','w')
writeLines(paste0(c('Site','Year','Species','N.butterflies','NumObservers','Party_Hours','Latitude','Longitude','Julian.date','grid_id','EcoI','EcoII','EcoregionI','EcoregionII'),collapse='\t'),con)
for (i in 1:nrow(site.years)){
	print(i)
	s1 = site.years[i,1]
	y1 = site.years[i,2]
	sydata = butterflies50.2[which(butterflies50.2[,1]==s1 & butterflies50.2[,2]==y1),]
	sydata$N.butterflies[which(is.na(sydata$N.butterflies))] = 0
	c12 = sydata[1,1:2]
	c515 = sydata[1,5:14]
	missing.species = species[which(is.na(match(species,unique(sydata$Species))))]
	nm = length(missing.species)
	add.data = data.frame('Site'=rep(c12[1,1],nm),'Year'=rep(c12[1,2],nm),'Species'=missing.species,'N.butterflies'=rep(0,nm),'NumObservers'=rep(c515[1,1],nm),
		'Party_Hours'=rep(c515[1,2],nm),'Latitude'=rep(c515[1,3],nm),'Longitude'=rep(c515[1,4],nm),'Julian.date'=rep(c515[1,5],nm),
		'grid_id'=rep(c515[1,6],nm),'EcoI'=rep(c515[1,7],nm),'EcoII'=rep(c515[1,8],nm),'EcoregionI'=rep(c515[1,9],nm),'EcoregionII'=rep(c515[1,10],nm),stringsAsFactors=F)
	for (j in 1:nrow(add.data)){
		writeLines(paste0(add.data[j,],collapse='\t'),con)
	}
}
close(con)
con = file('./data/add_50km_1993-2018_w0s.txt','r')
add.lines = readLines(con)
close(con)
add.data = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(add.data) = strsplit(add.lines[1],'\t')[[1]]
butterflies50.2 = data.frame(rbind(butterflies50.2,add.data),stringsAsFactors=F)
butterflies50.2$N.butterflies = as.numeric(butterflies50.2$N.butterflies)
butterflies50.2$NumObservers = as.numeric(butterflies50.2$NumObservers)
butterflies50.2$Party_Hours = as.numeric(butterflies50.2$Party_Hours)
butterflies50.2$Latitude = as.numeric(butterflies50.2$Latitude)
butterflies50.2$Longitude = as.numeric(butterflies50.2$Longitude)
write.table(butterflies50.2,'./data/butterfly_data_gridded_50km_1993-2018_w0s.txt',sep='\t',row.names=F)


# 50 km, implicit zeroes recorded as NA (consider them pseudoabsences)
butterflies50 = read.table('./data/butterfly_data_gridded_50km.txt',sep='\t',as.is=T,check.names=F,header=T); str(butterflies50)
butterflies50.2 = butterflies50[which(butterflies50$Year>1992),]; dim(butterflies50.2)
site.years = unique(butterflies50.2[,1:2])
species = unique(butterflies50.2$Species)
con = file('./data/add_50km_1993-2018_wNAs.txt','w')
writeLines(paste0(c('Site','Year','Species','N.butterflies','NumObservers','Party_Hours','Latitude','Longitude','Julian.date','grid_id','EcoI','EcoII','EcoregionI','EcoregionII'),collapse='\t'),con)
for (i in 1:nrow(site.years)){
	print(i)
	s1 = site.years[i,1]
	y1 = site.years[i,2]
	sydata = butterflies50.2[which(butterflies50.2[,1]==s1 & butterflies50.2[,2]==y1),]
	c12 = sydata[1,1:2]
	c515 = sydata[1,5:14]
	missing.species = species[which(is.na(match(species,unique(sydata$Species))))]
	nm = length(missing.species)
	add.data = data.frame('Site'=rep(c12[1,1],nm),'Year'=rep(c12[1,2],nm),'Species'=missing.species,'N.butterflies'=rep(NA,nm),'NumObservers'=rep(c515[1,1],nm),
		'Party_Hours'=rep(c515[1,2],nm),'Latitude'=rep(c515[1,3],nm),'Longitude'=rep(c515[1,4],nm),'Julian.date'=rep(c515[1,5],nm),
		'grid_id'=rep(c515[1,6],nm),'EcoI'=rep(c515[1,7],nm),'EcoII'=rep(c515[1,8],nm),'EcoregionI'=rep(c515[1,9],nm),'EcoregionII'=rep(c515[1,10],nm),stringsAsFactors=F)
	for (j in 1:nrow(add.data)){
		writeLines(paste0(add.data[j,],collapse='\t'),con)
	}
}
close(con)
con = file('./data/add_50km_1993-2018_wNAs.txt','r')
add.lines = readLines(con)
close(con)
add.data = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(add.data) = strsplit(add.lines[1],'\t')[[1]]
butterflies50.2 = data.frame(rbind(butterflies50.2,add.data),stringsAsFactors=F)
butterflies50.2$N.butterflies = as.numeric(butterflies50.2$N.butterflies)
butterflies50.2$NumObservers = as.numeric(butterflies50.2$NumObservers)
butterflies50.2$Party_Hours = as.numeric(butterflies50.2$Party_Hours)
butterflies50.2$Latitude = as.numeric(butterflies50.2$Latitude)
butterflies50.2$Longitude = as.numeric(butterflies50.2$Longitude)
write.table(butterflies50.2,'./data/butterfly_data_gridded_50km_1993-2018_wNAs.txt',sep='\t',row.names=F)
