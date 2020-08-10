
setwd('<your/working/directory>')

library(mvabund)
library(lattice)

data1 = read.table('./data/butterfly_traits_envars_trends_50km_wNAs_trim_forModelling.txt',sep='\t',as.is=T,check.names=F,header=T)
gids = sort(unique(data1$grid_id))
species = sort(unique(data1$Species))
traits = c('AdultSize','Diet.breadth.families')
envars = c('Cropland.2015','Cropland.diff','Built.2015','Built.diff','Precip.2018','Precip.diff','Temp.2018','Temp.diff')

# R matrix: rows=sites, columns=envars
R.mat = matrix(NA,nrow=length(gids),ncol=length(envars))
rownames(R.mat) = gids; colnames(R.mat) = envars
for (g in 1:length(gids)){
	for (e in 1:length(envars)){
		R.mat[g,e] = unique(data1[which(data1$grid_id==gids[g]),which(colnames(data1)==envars[e])])
	}
}

# L matrix: rows=sites, columns=species abundances
L.mat = matrix(0,nrow=length(gids),ncol=length(species)) #species absences coded as 0
rownames(L.mat) = gids; colnames(L.mat) = species
for (g in 1:length(gids)){
	for (s in 1:length(species)){
		val1 = data1$Abundance.trend[which(data1$grid_id==gids[g] & data1$Species==species[s])]
		if (length(val1)>0){
			L.mat[g,s] = val1
		} else {
		}
	}
}

# Q matrix: rows=species, columns=traits
Q.mat = matrix(NA,nrow=length(species),ncol=length(traits))
rownames(Q.mat) = species; colnames(Q.mat) = traits
for (s in 1:length(species)){
	for (r in 1:length(traits)){
		Q.mat[s,r] = unique(data1[which(data1$Species==species[s]),which(colnames(data1)==traits[r])])
	}
}


fit1 = traitglm(R=data.frame(R.mat),L=data.frame(L.mat),Q=data.frame(Q.mat),method='manylm')

fit1$fourth


