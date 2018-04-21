# load stuff -----
if(!require(cluster)){install.packages('cluster');library(cluster)}
if(!require(fpc)){install.packages('fpc');library(fpc)}
if(!require(maptools)){install.packages('maptools');library(maptools)}
if(!require(RColorBrewer)){install.packages('RColorBrewer');library(RColorBrewer)}
if(!require(tidyverse)){install.packages('tidyverse')};library(tidyverse)
load('../../data/2_classification/departamental_data_ZIKV_apr28_2017.RData')
map.muni = readShapePoly('../../data/2_classification/COL_adm2_RAPID.shp')
depts = unique(as.character(map.muni@data$COD_DEPTO))
map.dept = unionSpatialPolygons(map.muni,map.muni@data$COD_DEPTO)
map.natl = unionSpatialPolygons(map.muni,rep(1,nrow(map.muni@data)))
pop.dept = as.numeric(by(map.muni@data$sum,sapply(map.muni@data$COD_DEPTO,function(cd)which(unique(map.muni@data$COD_DEPTO)==cd)),sum))
names(pop.dept) = as.character(unique(map.muni@data$COD_DEPTO))
pop.dept = pop.dept[order(names(pop.dept))]
load('../../output/2_classification/dept/curve_data.RData')

names_dept = names(pop.dept[-c(which.zero,34)])

# load co-variates
load('../../data/3_elucidation/R0.RData')

muni_data = as.tibble(map.muni@data) %>% select(c(ID_ESPACIA,sum)) %>% rename(Pop = sum)  %>% 
  mutate(ID_ESPACIA = as.integer(as.character(ID_ESPACIA))) %>%
  left_join(as.tibble(x), by = "ID_ESPACIA")

muni_R0 =  left_join(muni_data, 
                     group_by(muni_data, COD_DEPTO) %>% summarize(POP_DEPTO = sum(Pop)),
                     by = "COD_DEPTO")

dept_R0 = group_by(muni_R0, COD_DEPTO) %>% summarize( R0 = sum(R0 * Pop / POP_DEPTO))


dept_aggregate_non_ts = read_csv('../../data/3_elucidation/dept_aggregate_non_ts.csv') %>%
  select(c(DEPTID,Wpop2015, UrbanPop, MeanGCP_2005USD)) %>% 
  mutate(PerUrban = UrbanPop / Wpop2015)

dept_aegypti = read_csv('../../data/3_elucidation/dept_Ae_aegypti_weeks_weighted.csv') %>% 
  mutate(Aegypti = rowMeans(.[grep("week",names(.))])) %>%
  select(c(DEPTID, Aegypti))

dept_ndvi_terra = read_csv('../../data/3_elucidation/dept_ndvi_terra_weekly_weighted.csv') %>%
  mutate(Ndvi_terra = rowMeans(.[grep('[0-9]+', names(.))])) %>%
  select(DEPTID,Ndvi_terra)

dept_ndvi_aqua = read_csv('../../data/3_elucidation/dept_ndvi_aqua_weekly_weighted.csv') %>%
  mutate(Ndvi_aqua = rowMeans(.[grep('[0-9]+', names(.))])) %>%
  select(DEPTID,Ndvi_aqua)

dept_tmean = read_csv('../../data/3_elucidation/dept_tmean_weekly_weighted.csv') %>%
  mutate(Tmean = rowMeans(.[grep('[0-9]+', names(.))])) %>%
  select(DEPTID,Tmean)

anova_results_table = data.frame(
  row.names =  c('Aegypti', 'Ndvi_terra', 'Ndvi_aqua', 'Tmean', 'PerUrban', 'Wpop2015', 'MeanGCP_2005USD', 'R0'),
  K2_F = rep(0,8), K2_P = 0, K3_F = 0, K3_P = 0, K4_F = 0, K4_P = 0
)

# assign groupings when there are k groups-----
for (k in 2:4){
  clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
  
  clusters_df = data.frame(DEPTID = as.integer(names_dept), Cluster = clusters) %>%
    left_join(dept_aggregate_non_ts, by = 'DEPTID') %>%
    left_join(dept_aegypti, by = 'DEPTID') %>%
    left_join(dept_ndvi_terra, by = 'DEPTID') %>%
    left_join(dept_ndvi_aqua, by = 'DEPTID') %>%
    left_join(dept_tmean, by = 'DEPTID') %>%
    left_join(dept_R0, by = c("DEPTID" = "COD_DEPTO"))
  print(dim(clusters_df))
  print(length(which.zero))
  for(rr in rownames(anova_results_table)){
    anova_results_table[rr, sprintf('K%.0f_F',k)] = anova(aov(
      as.formula(sprintf('Cluster ~ %s',rr)), data = clusters_df))$`F value`[1]
    anova_results_table[rr, sprintf('K%.0f_P',k)] = anova(aov(
      as.formula(sprintf('Cluster ~ %s',rr)), data = clusters_df))$`Pr(>F)`[1]
  }
}



write.csv( anova_results_table, '../../output/3_elucidation/anova_analysis_dept.csv')

#Figures K = 2------
k = 2
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering

clusters_df = data.frame(DEPTID = as.integer(names_dept), Cluster = clusters) %>%
  left_join(dept_aggregate_non_ts, by = 'DEPTID') %>%
  left_join(dept_aegypti, by = 'DEPTID') %>%
  left_join(dept_ndvi_terra, by = 'DEPTID') %>%
  left_join(dept_ndvi_aqua, by = 'DEPTID') %>%
  left_join(dept_tmean, by = 'DEPTID') %>%
  left_join(dept_R0, by = c("DEPTID" = "COD_DEPTO"))
layout(
  matrix(1:2,1,2)
)
par(mar=c(1,0,1,0), oma = c(2,4,2,2))

clusters_covariates_df = clusters_df %>% group_by(Cluster) %>% summarize(Tmean = mean(Tmean)) %>% ungroup()
boxplot(filter(clusters_df, Cluster == 1) %>% select(Tmean), 
        ylim = c(15,30), axes = F)
axis(2)
mtext("Cluster 1", side = 1, line = 2)
mtext("Temperature", side = 2, line = 2)
boxplot(filter(clusters_df, Cluster == 2) %>% select(Tmean),
        ylim = c(15,30), yaxt = 'n', axes =F)
mtext("Cluster 2", side = 1, line = 2)

write.csv( clusters_covariates_df, '../../output/3_elucidation/clusters_k2_mean_difference_dept.csv', row.names = F)
