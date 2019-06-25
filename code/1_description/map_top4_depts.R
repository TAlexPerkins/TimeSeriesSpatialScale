# load necessary R packages
if(!require(maptools)){install.packages('maptools');library(maptools)}
if(!require(RColorBrewer)){install.packages('RColorBrewer');library(RColorBrewer)}

# load maps and process cartograms
map.muni = readShapePoly('../../data/2_classification/COL_adm2_RAPID.shp')
depts = unique(as.character(map.muni@data$COD_DEPTO))
map.dept = unionSpatialPolygons(map.muni,map.muni@data$COD_DEPTO)

# specify top four departments that we want to plot
depts.top4 = rbind(c('76','Valle_Del_Cauca'),c('54','Norte_Santander'),c('68','Santander'),c('73','Tolima'))
cols.top4 = brewer.pal(n = 6, "Set2")[1:4]

# assign colors to each department
cols.all = rep('white',length(map.dept))
for(ii in 1:4){
  cols.all[which(names(map.dept)==depts.top4[ii,1])] = cols.top4[ii]
}

# plot a map with the top 4 departments colored appropriately
jpeg(filename='../../output/1_description/map_top4_depts.jpeg',width=5,height=5,units='in',res=500)
plot(map.dept,col=cols.all)
dev.off()
