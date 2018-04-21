# load stuff
if(!require(cluster)){install.packages('cluster');library(cluster)}
if(!require(fpc)){install.packages('fpc');library(fpc)}
if(!require(maptools)){install.packages('maptools');library(maptools)}
if(!require(RColorBrewer)){install.packages('RColorBrewer');library(RColorBrewer)}
if(!require(tidyverse)){install.packages('tidyverse')};library(tidyverse)
load('../../data/2_classification/data_municipalities_new_nov5.Rdata')
map.muni = readShapePoly('../../data/2_classification/COL_adm2_RAPID.shp')
depts = unique(as.character(map.muni@data$COD_DEPTO))
map.dept = unionSpatialPolygons(map.muni,map.muni@data$COD_DEPTO)
case.cums.muni = numeric(length=length(map.muni@data$ID_ESPACIA))
names(case.cums.muni) = as.character(map.muni@data$ID_ESPACIA)
load('../../output/2_classification/muni/curve_data.RData')


names_muni = names(munip.list)[-which.zero]

# load co-variates
load('../data/3_simulation/R0.RData')

muni_R0 = as.tibble(map.muni@data) %>% select(c(ID_ESPACIA,sum)) %>% rename(Pop = sum)  %>% 
  mutate(ID_ESPACIA = as.integer(as.character(ID_ESPACIA))) %>%
  left_join(as.tibble(x), by = "ID_ESPACIA")

muni_aggregate_non_ts = read_csv('../../data/3_simulation/municip_aggregate_non_ts.csv') %>%
  select(c(ID_ESPACIA,Wpop2015, UrbanPop, MeanGCP_2005USD)) %>% 
  mutate(PerUrban = UrbanPop / Wpop2015)

muni_aegypti = read_csv('../../data/3_simulation/municip_Ae_aegypti_weeks_weighted.csv') %>% 
  mutate(Aegypti = rowMeans(.[grep("week",names(.))])) %>%
  select(c(ID_ESPACIA, Aegypti))

muni_ndvi_terra = read_csv('../../data/3_simulation/municip_ndvi_terra_weekly_weighted.csv') %>%
  mutate(Ndvi_terra = rowMeans(.[grep('[0-9]+', names(.))])) %>%
  select(ID_ESPACIA,Ndvi_terra)

muni_ndvi_aqua = read_csv('../../data/3_simulation/municip_ndvi_aqua_weekly_weighted.csv') %>%
  mutate(Ndvi_aqua = rowMeans(.[grep('[0-9]+', names(.))])) %>%
  select(ID_ESPACIA,Ndvi_aqua)

muni_tmean = read_csv('../../data/3_simulation/municip_tmean_weekly_weighted.csv') %>%
  mutate(Tmean = rowMeans(.[grep('[0-9]+', names(.))])) %>%
  select(ID_ESPACIA,Tmean)

anova_results_table = data.frame(
  row.names =  c('Aegypti', 'Ndvi_terra', 'Ndvi_aqua', 'Tmean', 'PerUrban', 'Wpop2015', 'MeanGCP_2005USD', 'R0'),
  K2_F = rep(0,8), K2_P = 0, K3_F = 0, K3_P = 0, K4_F = 0, K4_P = 0
)

# assign groupings when there are k groups
for (k in 2:4){
  clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
  
  clusters_df = data.frame(ID_ESPACIA = as.integer(names_muni), Cluster = clusters) %>%
    left_join(muni_aggregate_non_ts, by = 'ID_ESPACIA') %>%
    left_join(muni_aegypti, by = 'ID_ESPACIA') %>%
    left_join(muni_ndvi_terra, by = 'ID_ESPACIA') %>%
    left_join(muni_ndvi_aqua, by = 'ID_ESPACIA') %>%
    left_join(muni_tmean, by = 'ID_ESPACIA') %>%
    left_join(muni_R0, by = "ID_ESPACIA") %>% drop_na()
  for(rr in rownames(anova_results_table)){
    anova_results_table[rr, sprintf('K%.0f_F',k)] = anova(aov(
      as.formula(sprintf('Cluster ~ %s',rr)), data = clusters_df))$`F value`[1]
    anova_results_table[rr, sprintf('K%.0f_P',k)] = anova(aov(
      as.formula(sprintf('Cluster ~ %s',rr)), data = clusters_df))$`Pr(>F)`[1]
  }
}

write.csv(anova_results_table, '../../output/3_simulation/anova_analysis_muni.csv')

#Figures K = 2------
k = 2
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering

clusters_df = data.frame(ID_ESPACIA = as.integer(names_muni), Cluster = clusters) %>%
  left_join(muni_aggregate_non_ts, by = 'ID_ESPACIA') %>%
  left_join(muni_aegypti, by = 'ID_ESPACIA') %>%
  left_join(muni_ndvi_terra, by = 'ID_ESPACIA') %>%
  left_join(muni_ndvi_aqua, by = 'ID_ESPACIA') %>%
  left_join(muni_tmean, by = 'ID_ESPACIA') %>%
  left_join(muni_R0, by = "ID_ESPACIA") %>% drop_na()


layout(
  matrix(1:2,1,2)
)
par(mar=c(1,0,1,0), oma = c(2,4,2,2))

clusters_covariates_df = clusters_df %>% group_by(Cluster) %>% summarize(Aegypti = mean(Aegypti)) %>% ungroup()
 
plot(0,0,ylim = c(min(clusters_df$Aegypti),max(clusters_df$Aegypti)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(Aegypti), names = 'Aegypti', axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("Aegypti", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$Aegypti),max(clusters_df$Aegypti)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(Aegypti), names = 'Aegypti', axes = F, col = "lightgray", add = T)
mtext("Cluser 2", side = 1, line = 2)

# Terra
clusters_covariates_df = clusters_df %>% group_by(Cluster) %>% summarize(Ndvi_terra = mean(Ndvi_terra))
plot(0,0,ylim = c(min(clusters_df$Ndvi_terra),max(clusters_df$Ndvi_terra)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(Ndvi_terra), names = 'Ndvi_terra', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("Ndvi_terra", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$Ndvi_terra),max(clusters_df$Ndvi_terra)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(Ndvi_terra), names = 'Ndvi_terra',
        axes = F, col = "lightgray", add = T)
mtext("Cluser 2", side = 1, line = 2)

#Aqua
clusters_covariates_df = left_join(clusters_covariates_df, 
                                   clusters_df %>% group_by(Cluster) %>% summarize(Ndvi_aqua = mean(Ndvi_aqua)) %>% ungroup(), by = 'Cluster')
plot(0,0,ylim = c(min(clusters_df$Ndvi_aqua),max(clusters_df$Ndvi_aqua)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(Ndvi_aqua), names = 'Ndvi_aqua',
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("Ndvi_aqua", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$Ndvi_aqua),max(clusters_df$Ndvi_aqua)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(Ndvi_aqua), names = 'Ndvi_aqua', 
        axes = F, col = "lightgray", add = T)
mtext("Cluser 2", side = 1, line = 2)

#Tmean
clusters_covariates_df = left_join(clusters_covariates_df, 
  clusters_df %>% group_by(Cluster) %>% summarize(Tmean = mean(Tmean)) %>% ungroup(), by = 'Cluster')
plot(0,0,ylim = c(min(clusters_df$Tmean),max(clusters_df$Tmean)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(Tmean), names = 'Tmean', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("Tmean", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$Tmean),max(clusters_df$Tmean)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(Tmean), names = 'Tmean', 
        axes = F, col = "lightgray", add = T)
mtext("Cluser 2", side = 1, line = 2)

# Urban
clusters_covariates_df = left_join(clusters_covariates_df, 
  clusters_df %>% group_by(Cluster) %>% summarize(PerUrban = mean(PerUrban)) %>% ungroup(), by = 'Cluster')
plot(0,0,ylim = c(min(clusters_df$PerUrban),max(clusters_df$PerUrban)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(PerUrban), names = 'PerUrban', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("PerUrban", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$PerUrban),max(clusters_df$PerUrban)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(PerUrban), names = 'PerUrban', 
        axes = F, col = "lightgray", add = T)
mtext("Cluser 2", side = 1, line = 2)

#Pop
clusters_covariates_df = left_join(clusters_covariates_df, 
  clusters_df %>% group_by(Cluster) %>% summarize(Wpop2015 = mean(Wpop2015)) %>% ungroup(), by = 'Cluster')
plot(0,0,ylim = c(min(clusters_df$Wpop2015),max(clusters_df$Wpop2015)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(Wpop2015), names = 'Wpop2015', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("Wpop2015", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$Wpop2015),max(clusters_df$Wpop2015)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(Wpop2015), names = 'Wpop2015',
         col = "lightgray", add = T,axes = F)
mtext("Cluser 2", side = 1, line = 2)

#GCP
clusters_covariates_df = left_join(clusters_covariates_df, 
  clusters_df %>% group_by(Cluster) %>% summarize(MeanGCP_2005USD = mean(MeanGCP_2005USD)) %>% ungroup(), by = 'Cluster')
plot(0,0,ylim = c(min(clusters_df$MeanGCP_2005USD),max(clusters_df$MeanGCP_2005USD)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(MeanGCP_2005USD), names = 'MeanGCP_2005USD', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("MeanGCP_2005USD", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$MeanGCP_2005USD),max(clusters_df$MeanGCP_2005USD)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(MeanGCP_2005USD), names = 'MeanGCP_2005USD',
        col = "lightgray", add = T,axes = F)
mtext("Cluser 2", side = 1, line = 2)

#R0
clusters_covariates_df = left_join(clusters_covariates_df, 
  clusters_df %>% group_by(Cluster) %>% summarize(R0 = mean(R0)) %>% ungroup(), by = 'Cluster')
plot(0,0,ylim = c(min(clusters_df$R0),max(clusters_df$R0)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(R0), names = 'R0', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("R0", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$R0),max(clusters_df$R0)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(R0), names = 'R0',
        col = "lightgray", add = T,axes = F)
mtext("Cluser 2", side = 1, line = 2)

write.csv( clusters_covariates_df, '../../output/3_elucidation/clusters_k2_mean_difference_muni.csv', row.names = F)
#Figures K = 3------
k = 3
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering

clusters_df = data.frame(ID_ESPACIA = as.integer(names_muni), Cluster = clusters) %>%
  left_join(muni_aggregate_non_ts, by = 'ID_ESPACIA') %>%
  left_join(muni_aegypti, by = 'ID_ESPACIA') %>%
  left_join(muni_ndvi_terra, by = 'ID_ESPACIA') %>%
  left_join(muni_ndvi_aqua, by = 'ID_ESPACIA') %>%
  left_join(muni_tmean, by = 'ID_ESPACIA') %>%
  left_join(muni_R0, by = "ID_ESPACIA") %>% drop_na()


layout(
  matrix(1:3,1,3)
)
par(mar=c(1,0,1,0), oma = c(2,4,2,2))

# Urban
clusters_covariates_df = clusters_df %>% group_by(Cluster) %>% summarize(PerUrban = mean(PerUrban)) %>% ungroup()
plot(0,0,ylim = c(min(clusters_df$PerUrban),max(clusters_df$PerUrban)), xlim = c(0,2),axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 1) %>% pull(PerUrban), names = 'PerUrban', 
        axes = F, col = "lightgray", add = T)
axis(2)
mtext("Cluser 1", side = 1, line = 2)
mtext("PerUrban", side = 2, line = 2)
plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$PerUrban),max(clusters_df$PerUrban)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 2) %>% pull(PerUrban), names = 'PerUrban', 
        axes = F, col = "lightgray", add = T)
mtext("Cluser 2", side = 1, line = 2)

plot(0,0,xlim = c(0,2),ylim = c(min(clusters_df$PerUrban),max(clusters_df$PerUrban)), axes = F, lty = 0, type = "l")
boxplot(filter(clusters_df, Cluster == 3) %>% pull(PerUrban), names = 'PerUrban', 
        axes = F, col = "lightgray", add = T)
mtext("Cluser 3", side = 1, line = 2)

write.csv( clusters_covariates_df, '../../output/3_elucidation/clusters_k3_mean_difference_muni.csv', row.names = F)
