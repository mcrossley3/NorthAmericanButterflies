

#################################################################################
# Check & Merge data (merge sub-species)

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
dat = read.csv('./data/NABA_Count_Data_through_2018.csv',as.is=T,check.names=F,header=T) #raw data available from NABA
dat = dat[which(!is.na(dat$NumSeen)),]
colnames(dat)[1] = 'CountName' #overwrite original column name, which contains invisible characters
#na.count = apply(dat,1,function(x){length(which(is.na(x)))})
lonlat = paste(dat$Lng,dat$Lat,sep='_')
dat$SiteLonLat = paste(dat$CountName,lonlat,sep='_')
str(dat)
dim(dat) #404,331 records
sum(dat$NumSeen,na.rm=T) #12,148,588 butterflies
length(unique(dat$CountName)) #1,221 sites

# Remove genus-level records
species.split = apply(array(dat$ScientificName),1,function(x){strsplit(x,' ')[[1]]})
sslen = unlist(lapply(species.split,function(x){length(x)}))
keep1 = which(sslen>1)
dat = dat[keep1,]; str(dat) #388,601 species-level records

# Get R-readable dates
dat$Date2 = apply(array(dat$Date),1,function(x){convert.date(x)}); str(dat$Date2) #R-readable date format
dat$Julian = yday(dat$Date2); str(dat$Julian)

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
#dat$GenusSpecies[which(dat$GenusSpecies=='Boloria montinus')] = 'Boloria chariclea' #NABA considers separate species
#dat$GenusSpecies[which(dat$GenusSpecies=='Hemiargus ammon')] = 'Cyclargus ammon' #NABA considers this Hemiargus ammon
#dat$GenusSpecies[which(dat$GenusSpecies=='Hemiargus thomasi')] = 'Cyclargus thomasi' #NABA considers this Hemiargus thomasi
dat$GenusSpecies[which(dat$GenusSpecies=='Colias scudderi')] = 'Colias scudderii' #fix typo
dat$GenusSpecies[which(dat$GenusSpecies=='Neophasia terlootii')] = 'Neophasia terlooii' #fix typo

uspecies = sort(unique(dat$GenusSpecies)) #644 unique species
species.merge = sort(unique(c(dat$GenusSpecies[which(sslen>2)],'Colias philodice','Limenitis archippus','Limenitis lorquini','Satyrium titus','Callophrys augustinus')))
length(species.merge) #45 species with records of separate subspecies to be merged
site.remove.pos = unlist(apply(array(species.merge),1,function(x){which(dat$GenusSpecies==x)}))
dat2 = dat[-site.remove.pos,]; dim(dat2)
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
dat2 = data.frame(rbind(dat2,add.dat)); dim(dat2) #387,842 records left


# Create a tidy table for exploration
butterflies = data.frame(
	'Site'=dat2$SiteLonLat,
	'Year'=apply(array(dat2$Date2),1,function(x){strsplit(x,'-')[[1]][1]}),
	'Species'=dat2$GenusSpecies,
	'N.butterflies'=dat2$NumSeen,
	'NumObservers'=dat2$NumObservers,
	'Party_Hours'=dat2$Party_Hours,
	'Latitude'=dat2$Lat,
	'Longitude'=dat2$Lng,
	'Julian.date'=dat2$Julian,
	stringsAsFactors=F)

write.table(butterflies,'./data/butterfly_data_NoMerge.txt',sep='\t',row.names=F)


################################################################################
# Explore data

library(rgdal)
library(sp)
library(raster)

# Import tidy merged dataset
butterflies = read.table('./data/butterfly_data_NoMerge.txt',sep='\t',as.is=T,check.names=F,header=T) #387,842
sum(butterflies$N.butterflies,na.rm=T) #11,950,146

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
	syears = sort(as.numeric(unique(s.dat$Year)))
	out[s,1] = sites[s]
	out[s,2] = n.species
	out[s,3] = sum.count
	out[s,4] = sum.partyhours
	out[s,5] = lon
	out[s,6] = lat
	out[s,7] = syears[1]
	out[s,8] = rev(syears)[1]
	out[s,9] = (rev(syears)[1] - syears[1]) + 1
}
str(out)
write.table(out,'site_summary_NoMerge.txt',sep='\t',row.names=F)

#out = read.table('site_summary_NoMerge.txt',sep='\t',as.is=T,check.names=F,header=T)

# Write site summary to shapefile
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
shp = spTransform(SpatialPointsDataFrame(coords=cbind(as.numeric(out$Lon),as.numeric(out$Lat)),data=out,proj4string=CRS(wgs84)),CRS(proj1))
eco2 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l2/NA_CEC_Eco_Level2.shp',layer='NA_CEC_Eco_Level2',verbose=F,stringsAsFactors=F),CRS(proj1))
#over1 = over(shp,eco1)
over2 = over(shp,eco2)
shp@data$EcoI = over2$NA_L1CODE
shp@data$EcoII = over2$NA_L2CODE
shp@data$EcoregionI = over2$NA_L1NAME
shp@data$EcoregionII = over2$NA_L2NAME
# Fix NA and "WATER" ecoregion assignments
shp2 = shp[which(is.na(over2[,2])),]
shp2.5km = buffer(shp2,width=50000,dissolve=F)
ov1 = over(shp2.5km,eco2)
shp@data$EcoI[match(shp2@data$Site,shp@data$Site)] = ov1$NA_L1CODE
shp@data$EcoII[match(shp2@data$Site,shp@data$Site)] = ov1$NA_L2CODE
shp@data$EcoregionI[match(shp2@data$Site,shp@data$Site)] = ov1$NA_L1NAME
shp@data$EcoregionII[match(shp2@data$Site,shp@data$Site)] = ov1$NA_L2NAME
writeOGR(shp,dsn='./shapefiles/butterfly_sites_NoMerge.shp',layer='butterfly_sites_NoMerge',driver='ESRI Shapefile',overwrite=T)

##############
# Filter by >1992
butterflies = read.table('./data/butterfly_data_NoMerge.txt',sep='\t',as.is=T,check.names=F,header=T); dim(butterflies) #387,842 records
butterflies2 = butterflies[which(butterflies$Year>1992),]; dim(butterflies2) #353,933 records

##############
# Filter by Julian date
# June-August 152:243
# July 182:212
hist(butterflies2$Julian.date)
100 * length(which(butterflies2$Julian.date >= 152 & butterflies2$Julian.date <= 243)) / nrow(butterflies2) #91% kept
100 * length(which(butterflies2$Julian.date >= 182 & butterflies2$Julian.date <= 212)) / nrow(butterflies2) #59% kept

butterflies.JunAug = butterflies2[which(butterflies2$Julian.date >= 152 & butterflies2$Julian.date <= 243),]
butterflies.Jul = butterflies2[which(butterflies2$Julian.date >= 182 & butterflies2$Julian.date <= 212),]

###############
# Filter by series length
shp = readOGR(dsn='./shapefiles/butterfly_sites_NoMerge.shp',layer='butterfly_sites_NoMerge',verbose=F,stringsAsFactors=F)
dp = which(shp@data$N_years>9)
shp2 = shp[dp,]; dim(shp2@data)
length(unique(shp@data$Site)); length(unique(shp2@data$Site)) #539 out of 1,270 sites left
writeOGR(shp2,dsn='./shapefiles/butterfly_sites_NoMerge_10yr.shp',layer='butterfly_sites_NoMerge_10yr',driver='ESRI Shapefile',overwrite=T)
deep_series = shp@data$Site[dp]
butterflies.JunAug2 = butterflies.JunAug[which(!is.na(match(butterflies.JunAug$Site,deep_series))),]; dim(butterflies.JunAug2) #273,099 records remain
butterflies.Jul2 = butterflies.Jul[which(!is.na(match(butterflies.Jul$Site,deep_series))),]; dim(butterflies.Jul2) #177,129 records remain

###############
# Assign to grid cell ids

library(rgdal)
library(scales)
library(sp)
library(rgeos)
library(maptools)
library(ggplot2)
library(spdep)
library(raster)

# Overlay 50x50 km grid
butterfly_circles = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_NoMerge_10yr.shp',layer='butterfly_sites_NoMerge_10yr',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile with butterfly count circles
eco2 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l2/NA_CEC_Eco_Level2.shp',layer='NA_CEC_Eco_Level2',verbose=F,stringsAsFactors=F),CRS(proj1)) #shapefile available at https://www.epa.gov/eco-research/ecoregions-north-america
grid50 = readOGR(dsn='./shapefiles/50km_grid_clip.shp',layer='50km_grid_clip',verbose=F,stringsAsFactors=F)
grid50@data$grid_id = seq(1,length(grid50),1)
centroids = gCentroid(grid50,byid=T)
grid50@data$centroid_x = centroids@coords[,1]
grid50@data$centroid_y = centroids@coords[,2]
grid50@data = grid50@data[,-1]
writeOGR(grid50,dsn='./shapefiles/butterfly_grid_50km_NoMerge.shp',layer='butterfly_grid_50km_NoMerge',driver='ESRI Shapefile',overwrite=T)
bover = over(butterfly_circles,grid50)
butterfly_circles$grid_id = bover$grid_id
writeOGR(butterfly_circles,dsn='./shapefiles/butterfly_sites_50km_NoMerge.shp',layer='butterfly_sites_50km_NoMerge',driver='ESRI Shapefile',overwrite=T)
gids = sort(unique(butterfly_circles@data$grid_id))

butterflies.JunAug2$grid_id = butterflies.JunAug2$EcoI = butterflies.JunAug2$EcoII = butterflies.JunAug2$EcoregionI = butterflies.JunAug2$EcoregionII = NA
for (g in 1:length(gids)){
	print(g)
	site1 = butterfly_circles@data$Site[which(butterfly_circles@data$grid_id==gids[g])]
	eco1 = butterfly_circles@data$EcoI[which(butterfly_circles@data$grid_id==gids[g])]
	eco2 = butterfly_circles@data$EcoII[which(butterfly_circles@data$grid_id==gids[g])]
	ecoregion1 = butterfly_circles@data$EcorgnI[which(butterfly_circles@data$grid_id==gids[g])]
	ecoregion2 = butterfly_circles@data$EcrgnII[which(butterfly_circles@data$grid_id==gids[g])]
	for (i in 1:length(site1)){
		row.pos = which(butterflies.JunAug2$Site==site1[i])
		butterflies.JunAug2$grid_id[row.pos] = gids[g]	
		butterflies.JunAug2$EcoI[row.pos] = eco1
		butterflies.JunAug2$EcoII[row.pos] = eco2	
		butterflies.JunAug2$EcoregionI[row.pos] = ecoregion1	
		butterflies.JunAug2$EcoregionII[row.pos] = ecoregion2	
	}
}
sort(unique(butterflies.JunAug2$Species)) #529 species
write.table(butterflies.JunAug2,'./data/butterfly_data_gridded_50km_NoMergeJunAug.txt',sep='\t',row.names=F)

butterflies.Jul2$grid_id = butterflies.Jul2$EcoI = butterflies.Jul2$EcoII = butterflies.Jul2$EcoregionI = butterflies.Jul2$EcoregionII = NA
for (g in 1:length(gids)){
	site1 = butterfly_circles@data$Site[which(butterfly_circles@data$grid_id==gids[g])]
	eco1 = butterfly_circles@data$EcoI[which(butterfly_circles@data$grid_id==gids[g])]
	eco2 = butterfly_circles@data$EcoII[which(butterfly_circles@data$grid_id==gids[g])]
	ecoregion1 = butterfly_circles@data$EcorgnI[which(butterfly_circles@data$grid_id==gids[g])]
	ecoregion2 = butterfly_circles@data$EcrgnII[which(butterfly_circles@data$grid_id==gids[g])]
	for (i in 1:length(site1)){
		row.pos = which(butterflies.Jul2$Site==site1[i])
		butterflies.Jul2$grid_id[row.pos] = gids[g]	
		butterflies.Jul2$EcoI[row.pos] = eco1
		butterflies.Jul2$EcoII[row.pos] = eco2	
		butterflies.Jul2$EcoregionI[row.pos] = ecoregion1	
		butterflies.Jul2$EcoregionII[row.pos] = ecoregion2	
	}
}
sort(unique(butterflies.Jul2$Species)) #499 species
write.table(butterflies.Jul2,'./data/butterfly_data_gridded_50km_NoMergeJul.txt',sep='\t',row.names=F)

# Prepare neighbor graph for trend model input
nb50 <- poly2nb(grid50, row.names=grid50$grid_id); nb50
is.symmetric.nb(nb50, verbose = FALSE, force = TRUE)
nb2INLA("nb50.graph", nb50)

# Prepare grid cells with sites
gid.junaug = unique(butterflies.JunAug2$grid_id)
grid50.junaug = grid50[match(gid.junaug,grid50@data$grid_id),]
writeOGR(grid50.junaug,dsn='./shapefiles/butterfly_sites_50km_NoMergeJunAug.shp',layer='butterfly_sites_50km_NoMergeJunAug',driver='ESRI Shapefile',overwrite=T)
butterfly_circles.junaug = butterfly_circles[match(gid.junaug,butterfly_circles@data$grid_id),]
writeOGR(butterfly_circles.junaug,dsn='./shapefiles/butterfly_sites_50km_NoMergeJunAug.shp',layer='butterfly_sites_50km_NoMergeJunAug',driver='ESRI Shapefile',overwrite=T)

gid.jul = unique(butterflies.Jul2$grid_id)
grid50.jul = grid50[match(gid.jul,grid50@data$grid_id),]
writeOGR(grid50.jul,dsn='./shapefiles/butterfly_sites_50km_NoMergeJul.shp',layer='butterfly_sites_50km_NoMergeJul',driver='ESRI Shapefile',overwrite=T)
butterfly_circles.jul = butterfly_circles[match(gid.jul,butterfly_circles@data$grid_id),]
writeOGR(butterfly_circles.jul,dsn='./shapefiles/butterfly_sites_50km_NoMergeJul.shp',layer='butterfly_sites_50km_NoMergeJul',driver='ESRI Shapefile',overwrite=T)


####################################
# Remove species*sites that have < 5 data points

# JunAug
species = sort(unique(butterflies.JunAug2$Species))
ssjunaug = data.frame('Species'=NA,'Site'=NA,'N'=-999)
sx = 1
for (i in 1:length(species)){
	noquote(print(species[i]))
	sites = sort(unique(butterflies.JunAug2$Site[which(butterflies.JunAug2$Species==species[i])]))
	for (j in 1:length(sites)){
		ssjunaug[sx,1] = species[i]
		ssjunaug[sx,2] = sites[j]
		ssjunaug[sx,3] = length(which(butterflies.JunAug2$Species==species[i] & butterflies.JunAug2$Site==sites[j]))
		sx = sx + 1
	}
}
nrow(ssjunaug)
length(which(ssjunaug$N<5)); 100 * length(which(ssjunaug$N<5)) / nrow(ssjunaug) #39% of data lost

# Jul
species = sort(unique(butterflies.Jul2$Species))
ssjul = data.frame('Species'=NA,'Site'=NA,'N'=-999)
sx = 1
for (i in 1:length(species)){
	noquote(print(species[i]))
	sites = sort(unique(butterflies.Jul2$Site[which(butterflies.Jul2$Species==species[i])]))
	for (j in 1:length(sites)){
		ssjul[sx,1] = species[i]
		ssjul[sx,2] = sites[j]
		ssjul[sx,3] = length(which(butterflies.Jul2$Species==species[i] & butterflies.Jul2$Site==sites[j]))
		sx = sx + 1
	}
}
nrow(ssjul)
length(which(ssjul$N<5)); 100 * length(which(ssjul$N<5)) / nrow(ssjul) #50% of data lost


keep.junaug = ssjunaug[which(ssjunaug$N>4),]
keep.junaug.row = c()
for (i in 1:nrow(keep.junaug)){
	pos = which(butterflies.JunAug2$Site==keep.junaug$Site[i] & butterflies.JunAug2$Species==keep.junaug$Species[i])
	keep.junaug.row = c(keep.junaug.row,pos)
}

keep.jul = ssjul[which(ssjul$N>4),]
keep.jul.row = c()
for (i in 1:nrow(keep.jul)){
	pos = which(butterflies.Jul2$Site==keep.jul$Site[i] & butterflies.Jul2$Species==keep.jul$Species[i])
	keep.jul.row = c(keep.jul.row,pos)
}

butterflies.JunAug3 = butterflies.JunAug2[keep.junaug.row,]; dim(butterflies.JunAug3) #245,971 records
butterflies.Jul3 = butterflies.Jul2[keep.jul.row,]; dim(butterflies.Jul3) #149,954 records
write.table(butterflies.JunAug3,'./data/butterfly_data_gridded_50km_NoMergeJunAug_m5.txt',sep='\t',row.names=F)
write.table(butterflies.Jul3,'./data/butterfly_data_gridded_50km_NoMergeJul_m5.txt',sep='\t',row.names=F)

# Prepare grid cells with sites
gid.junaug = unique(butterflies.JunAug3$grid_id)
grid50.junaug = grid50[match(gid.junaug,grid50@data$grid_id),]
writeOGR(grid50.junaug,dsn='./shapefiles/butterfly_grid_50km_NoMergeJunAug_m5.shp',layer='butterfly_grid_50km_NoMergeJunAug_m5',driver='ESRI Shapefile',overwrite=T)
butterfly_circles.junaug = butterfly_circles[match(gid.junaug,butterfly_circles@data$grid_id),]
writeOGR(butterfly_circles.junaug,dsn='./shapefiles/butterfly_sites_50km_NoMergeJunAug_m5.shp',layer='butterfly_sites_50km_NoMergeJunAug_m5',driver='ESRI Shapefile',overwrite=T)

gid.jul = unique(butterflies.Jul3$grid_id)
grid50.jul = grid50[match(gid.jul,grid50@data$grid_id),]
writeOGR(grid50.jul,dsn='./shapefiles/butterfly_grid_50km_NoMergeJul_m5.shp',layer='butterfly_grid_50km_NoMergeJul_m5',driver='ESRI Shapefile',overwrite=T)
butterfly_circles.jul = butterfly_circles[match(gid.jul,butterfly_circles@data$grid_id),]
writeOGR(butterfly_circles.jul,dsn='./shapefiles/butterfly_sites_50km_NoMergeJul_m5.shp',layer='butterfly_sites_50km_NoMergeJul_m5',driver='ESRI Shapefile',overwrite=T)

check = butterflies.Jul3[which(butterflies.Jul3$Species=='Achlyodes thraso'),]
check2 = check[which(check$Site=='Bentsen-Rio Grande Valley State Park, TX._-98.3135_26.1915'),]
plot(check2$Year,check2$N.butterflies,type='l')


# Map sample sites in full and JunAug set
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj1 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
out = read.table('site_year_summary.txt',sep='\t',as.is=T,check.names=F,header=T)
shp = spTransform(SpatialPointsDataFrame(coords=cbind(as.numeric(out$Lon),as.numeric(out$Lat)),data=out,proj4string=CRS(wgs84)),CRS(proj1))
eco2 = spTransform(readOGR(dsn='./shapefiles/na_cec_eco_l2/NA_CEC_Eco_Level2.shp',layer='NA_CEC_Eco_Level2',verbose=F,stringsAsFactors=F),CRS(proj1))
png('./plots/map of all sites.png',res=300,width=480*6,height=480*6)
plot(shp,pch=16,col='white')
plot(eco2,col='grey70',border='grey70',add=T)
plot(shp,pch=16,add=T,col='black')
dev.off()
# Filtered sites, highlighting places with two dates of sampling within a year
shp = spTransform(readOGR(dsn='./shapefiles/butterfly_sites_50km_NoMergeJunAug.shp',layer='butterfly_sites_50km_NoMergeJunAug',stringsAsFactors=F),CRS(proj1))
png('./plots/map of JunAug sites.png',res=300,width=480*6,height=480*6)
plot(eco2,col='grey70',border='grey70')
plot(shp,pch=16,add=T,col='black')
dev.off()


########################################
# Create implicit zeroes to data set for richness & evenness trend estimation

setwd('C:/Users/mcros/Desktop/Postdoc UGA/NABA_butterflies')

# Jun-Aug set
butterfliesJunAug = read.table('./data/butterfly_data_gridded_50km_NoMergeJunAug_m5.txt',sep='\t',as.is=T,check.names=F,header=T); str(butterfliesJunAug) #310,329 obs
sum(butterfliesJunAug$N.butterflies) #8,448,945
length(unique(butterfliesJunAug$Species)) #456
length(unique(butterfliesJunAug$Site)) #503

site.years = unique(butterfliesJunAug[,1:2])
species = unique(butterfliesJunAug$Species)
con = file('./data/add_NoMergeJunAug_w0s_m5.txt','w')
writeLines(paste0(c('Site','Year','Species','N.butterflies','NumObservers','Party_Hours','Latitude','Longitude','Julian.date','grid_id','EcoI','EcoII','EcoregionI','EcoregionII'),collapse='\t'),con)
for (i in 1:nrow(site.years)){
	print(i)
	s1 = site.years[i,1]
	y1 = site.years[i,2]
	sydata = butterfliesJunAug[which(butterfliesJunAug[,1]==s1 & butterfliesJunAug[,2]==y1),]
	sydata$N.butterflies[which(is.na(sydata$N.butterflies))] = 0
	missing.species = species[which(is.na(match(species,unique(sydata$Species))))]
	nm = length(missing.species)
	add.data = data.frame('Site'=rep(sydata$Site[1],nm),'Year'=rep(sydata$Year[1],nm),'Species'=missing.species,'N.butterflies'=rep(0,nm),'NumObservers'=rep(sydata$NumObservers[1],nm),
		'Party_Hours'=rep(sydata$Party_Hours[1],nm),'Latitude'=rep(sydata$Latitude[1],nm),'Longitude'=rep(sydata$Longitude[1],nm),'Julian.date'=rep(sydata$Julian.date[1],nm),
		'grid_id'=rep(sydata$grid_id[1],nm),'EcoI'=rep(sydata$EcoI[1],nm),'EcoII'=rep(sydata$EcoII[1],nm),'EcoregionI'=rep(sydata$EcoregionI[1],nm),'EcoregionII'=rep(sydata$EcoregionII[1],nm),stringsAsFactors=F)
	for (j in 1:nrow(add.data)){
		writeLines(paste0(add.data[j,],collapse='\t'),con)
	}
}
close(con)
con = file('./data/add_NoMergeJunAug_w0s_m5.txt','r')
add.lines = readLines(con)
close(con)
add.data = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(add.data) = strsplit(add.lines[1],'\t')[[1]]
butterfliesJunAug = data.frame(rbind(butterfliesJunAug,add.data),stringsAsFactors=F)
butterfliesJunAug$N.butterflies = as.numeric(butterfliesJunAug$N.butterflies)
butterfliesJunAug$NumObservers = as.numeric(butterfliesJunAug$NumObservers)
butterfliesJunAug$Party_Hours = as.numeric(butterfliesJunAug$Party_Hours)
butterfliesJunAug$Latitude = as.numeric(butterfliesJunAug$Latitude)
butterfliesJunAug$Longitude = as.numeric(butterfliesJunAug$Longitude)
str(butterfliesJunAug)
write.table(butterfliesJunAug,'./data/butterfly_data_NoMergeJunAug_w0s_m5.txt',sep='\t',row.names=F)

# Jul set
butterfliesJul = read.table('./data/butterfly_data_gridded_50km_NoMergeJul_m5.txt',sep='\t',as.is=T,check.names=F,header=T); str(butterfliesJul) #310,329 obs
site.years = unique(butterfliesJul[,1:2])
species = unique(butterfliesJul$Species)
con = file('./data/add_NoMergeJul_w0s_m5.txt','w')
writeLines(paste0(c('Site','Year','Species','N.butterflies','NumObservers','Party_Hours','Latitude','Longitude','Julian.date','grid_id','EcoI','EcoII','EcoregionI','EcoregionII'),collapse='\t'),con)
for (i in 1:nrow(site.years)){
	print(i)
	s1 = site.years[i,1]
	y1 = site.years[i,2]
	sydata = butterfliesJul[which(butterfliesJul[,1]==s1 & butterfliesJul[,2]==y1),]
	sydata$N.butterflies[which(is.na(sydata$N.butterflies))] = 0
	missing.species = species[which(is.na(match(species,unique(sydata$Species))))]
	nm = length(missing.species)
	add.data = data.frame('Site'=rep(sydata$Site[1],nm),'Year'=rep(sydata$Year[1],nm),'Species'=missing.species,'N.butterflies'=rep(0,nm),'NumObservers'=rep(sydata$NumObservers[1],nm),
		'Party_Hours'=rep(sydata$Party_Hours[1],nm),'Latitude'=rep(sydata$Latitude[1],nm),'Longitude'=rep(sydata$Longitude[1],nm),'Julian.date'=rep(sydata$Julian.date[1],nm),
		'grid_id'=rep(sydata$grid_id[1],nm),'EcoI'=rep(sydata$EcoI[1],nm),'EcoII'=rep(sydata$EcoII[1],nm),'EcoregionI'=rep(sydata$EcoregionI[1],nm),'EcoregionII'=rep(sydata$EcoregionII[1],nm),stringsAsFactors=F)
	for (j in 1:nrow(add.data)){
		writeLines(paste0(add.data[j,],collapse='\t'),con)
	}
}
close(con)
con = file('./data/add_NoMergeJul_w0s_m5.txt','r')
add.lines = readLines(con)
close(con)
add.data = data.frame(t(apply(array(add.lines[-1]),1,function(x){strsplit(x,'\t')[[1]]})),stringsAsFactors=F)
colnames(add.data) = strsplit(add.lines[1],'\t')[[1]]
butterfliesJul = data.frame(rbind(butterfliesJul,add.data),stringsAsFactors=F)
butterfliesJul$N.butterflies = as.numeric(butterfliesJul$N.butterflies)
butterfliesJul$NumObservers = as.numeric(butterfliesJul$NumObservers)
butterfliesJul$Party_Hours = as.numeric(butterfliesJul$Party_Hours)
butterfliesJul$Latitude = as.numeric(butterfliesJul$Latitude)
butterfliesJul$Longitude = as.numeric(butterfliesJul$Longitude)
str(butterfliesJul)
write.table(butterfliesJul,'./data/butterfly_data_NoMergeJul_w0s_m5.txt',sep='\t',row.names=F)


#####################################################
# What proportion of counts do Thymelicus lineola and Pieris rapae account for?
butterfly_counts = read.table('./data/butterfly_data_gridded_50km_NoMergeJunAug.txt',sep='\t',as.is=T,check.names=F,header=T)

# Proportion of counts that are top 2 invasives
N = sum(butterfly_counts$N.butterflies,na.rm=T)
Nabund = sum(butterfly_counts$N.butterflies[which(butterfly_counts$Species=='Thymelicus lineola' | butterfly_counts$Species=='Pieris rapae' | butterfly_counts$Species=='Colias philodice' | butterfly_counts$Species=='Phycoides tharos')],na.rm=T)
Nabund / N #27.5% of JunAug; 30.6% of Jul



