# load necessary R packages
if(!require(cluster)){install.packages('cluster');library(cluster)}
if(!require(fpc)){install.packages('fpc');library(fpc)}
if(!require(ggfortify)){install.packages('ggfortify');library(ggfortify)}
if(!require(maptools)){install.packages('maptools');library(maptools)}
if(!require(RColorBrewer)){install.packages('RColorBrewer');library(RColorBrewer)}
library(getcartr)



# load departmental level case report data
load('../../data/2_classification/departamental_data_ZIKV_apr28_2017.RData')



# load simulated data
load('../../output/3_elucidation/simulations.RData')



# load municipality level case report data and shape file
load('../../data/2_classification/data_municipalities_new_nov5.Rdata')
map.muni = readShapePoly('../../data/2_classification/COL_adm2_RAPID.shp')
depts = unique(as.character(map.muni@data$COD_DEPTO))
map.dept = unionSpatialPolygons(map.muni,map.muni@data$COD_DEPTO)
carto.pop.muni = quick.carto(map.muni,map.muni@data$sum)
carto.pop.muni.tmp = gBuffer(carto.pop.muni, byid=TRUE, width=0)
carto.pop.dept = unionSpatialPolygons(carto.pop.muni.tmp,carto.pop.muni.tmp@data$COD_DEPTO)
case.cums.muni = numeric(length=length(map.muni@data$ID_ESPACIA))
names(case.cums.muni) = as.character(map.muni@data$ID_ESPACIA)
for(ml in names(case.cums.muni)){
  if(ml %in% names(munip.list)){
    case.cums.muni[ml] = sum(sapply(
      which(names(munip.list)==ml),function(ii)sum(munip.list[[ii]]$cases_combo)))
  }
}
case.cums.muni[is.na(case.cums.muni)] = 0
case.cums.muni[case.cums.muni==0] = 1
carto.cases.muni = quick.carto(map.muni,case.cums.muni)
carto.cases.muni.tmp = gBuffer(carto.cases.muni, byid=TRUE, width=0)
carto.cases.dept = unionSpatialPolygons(carto.cases.muni.tmp,carto.cases.muni.tmp@data$COD_DEPTO)



# function for ggplot colors
gg_color_hue <- function(n, k, alpha.in) {
  hues = seq(15, 375, length = n + 1)
  sapply(1:length(k), function(ii)
    hcl(h = hues, l = 65, c = 100, alpha = alpha.in[ii])[1:n][k[ii]])
}



# record data from simulations in an array
sims.mat = array(0,dim=c(length(sims.list),121,100))
for(ss in 1:100){
  for(ii in 1:length(sims.list)){
    sims.mat[ii,,ss] = subset(sims.list[[ii]],sim==ss)$C
  }
}

# fit S curves to simulated data
curve.data = array(0,dim=c(dim(sims.mat)[1],7,100))
which.zero = matrix(0,dim(sims.mat)[1],100)
for(ss in 1:100){
  for(ii in 1:dim(sims.mat)[1]){
    if(sum(sims.mat[ii,,ss]) > 0){
      cases.norm = sims.mat[ii,,ss] / sum(sims.mat[ii,,ss])
      cases.norm.cum = cumsum(sims.mat[ii,,ss]) / sum(sims.mat[ii,,ss])
      tmp = optim(c(30,1),function(par){
        obs = cases.norm.cum
        pred = pnorm(1:length(cases.norm.cum),par[1],par[2])
        R2 = 1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2)
        return(-R2)})
      curve.data[ii,1,ss] = tmp$par[1]
      curve.data[ii,2,ss] = tmp$par[2]
      curve.data[ii,3,ss] = -tmp$value
      curve.data[ii,4,ss] =
        max(which(0.05 >
                    cases.norm.cum)) -
        qnorm(0.05,tmp$par[1],tmp$par[2])
      curve.data[ii,5,ss] =
        min(which(0.95 <
                    cases.norm.cum)) -
        qnorm(0.95,tmp$par[1],tmp$par[2])
      curve.data[ii,6,ss] =
        length(min(which(sims.mat[ii,,ss] > 0)):max(which(sims.mat[ii,,ss] > 0)))
      curve.data[ii,7,ss] =
        1 - sum(cases.norm[
          min(which(sims.mat[ii,,ss] > 0)):max(which(sims.mat[ii,,ss] > 0))] > 0) / curve.data[ii,6,ss]
    }
    if(sum(sims.mat[ii,,ss]) == 0){
      which.zero[ii,ss] = 1
    }
  }
  curve.data[is.infinite(curve.data)] = 0
  colnames(curve.data) = c(
    'Mean of fitted normal cdf',
    'Standard deviation of fitted normal cdf',
    'Coefficient of determination of fitted normal cdf',
    'Difference between 5% quantiles of empirical cdf and fitted normal cdf',
    'Difference between 95% quantiles of empirical cdf and fitted normal cdf',
    'Weeks between first and last case reports',
    'Proportion of weeks between first and last case reports with no case reports')
}
save(curve.data,which.zero,file='../../output/3_elucidation/curve_data_sim.RData')



# compute distributions of various statistics about the simulations and the classification algorithm applied to them
Pr.correct.source = Pr.correct.sink = Pr.correct.Rgt1 = Pr.correct.Rlt1 = rep(0,100)
Pr.source = Pr.sink = Pr.Rgt1 = Pr.Rlt1 = pct.zero = pct.zero.Rgt1 = rep(0,100)
table.output.2 = table.output.3 = list()
sil.wid.2 = sil.wid.3 = rep(0,100)
k.highest = rep(0,100)
for(ss in 1:100){
  k.highest[ss] = pamk(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])))$nc
  clusters = pam(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),2)$clustering
  Pr.correct.source[ss] = sum(x$R0[which(which.zero[,ss]==0)][which(clusters==1)] >= 1) / sum(clusters==1)
  Pr.correct.sink[ss] = sum(x$R0[which(which.zero[,ss]==0)][which(clusters==2)] < 1) / sum(clusters==2)
  Pr.correct.Rgt1[ss] = sum(x$R0[which(which.zero[,ss]==0)][which(clusters==1)] >= 1) / sum(x$R0[which(which.zero[,ss]==0)] >= 1)
  Pr.correct.Rlt1[ss] = sum(x$R0[which(which.zero[,ss]==0)][which(clusters==2)] < 1) / sum(x$R0[which(which.zero[,ss]==0)] < 1)
  Pr.source[ss] = sum(clusters == 1) / length(clusters)
  Pr.sink[ss] = sum(clusters == 2) / length(clusters)
  Pr.Rgt1[ss] = sum(x$R0[which(which.zero[,ss]==0)] >= 1) / length(clusters)
  Pr.Rlt1[ss] = sum(x$R0[which(which.zero[,ss]==0)] < 1) / length(clusters)
  pct.zero[ss] = sum(which.zero[,ss]==1) / nrow(which.zero)
  pct.zero.Rgt1[ss] = sum(which.zero[,ss]==1 & x$R0 >= 1) / sum(which.zero[,ss]==1)
  
  clusters.2 = pam(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),2)$clustering
  table.output.2[[ss]] = table(clusters.2,x$R0[which(which.zero[,ss]==0)] >= 1)

  clusters.3 = pam(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),3)$clustering
  table.output.3[[ss]] = table(clusters.3,x$R0[which(which.zero[,ss]==0)] >= 1)
  
  sil.wid.2[ss] = pamk(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),2)$crit[2]
  sil.wid.3[ss] = pamk(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),3)$crit[3]
}

source.2.1.1 = sapply(1:100,function(ss)table.output.2[[ss]][1,1]/sum(table.output.2[[ss]][,1]))
source.2.2.1 = sapply(1:100,function(ss)table.output.2[[ss]][2,1]/sum(table.output.2[[ss]][,1]))
source.3.1.1 = sapply(1:100,function(ss)table.output.3[[ss]][1,1]/sum(table.output.3[[ss]][,1]))
source.3.2.1 = sapply(1:100,function(ss)table.output.3[[ss]][2,1]/sum(table.output.3[[ss]][,1]))
source.3.3.1 = sapply(1:100,function(ss)table.output.3[[ss]][3,1]/sum(table.output.3[[ss]][,1]))
source.2.1.2 = sapply(1:100,function(ss)table.output.2[[ss]][1,2]/sum(table.output.2[[ss]][,2]))
source.2.2.2 = sapply(1:100,function(ss)table.output.2[[ss]][2,2]/sum(table.output.2[[ss]][,2]))
source.3.1.2 = sapply(1:100,function(ss)table.output.3[[ss]][1,2]/sum(table.output.3[[ss]][,2]))
source.3.2.2 = sapply(1:100,function(ss)table.output.3[[ss]][2,2]/sum(table.output.3[[ss]][,2]))
source.3.3.2 = sapply(1:100,function(ss)table.output.3[[ss]][3,2]/sum(table.output.3[[ss]][,2]))



# plot pairwise combinations of curve features stratified by grouping for one simulated dataset
ss = 1
for(k in 2:6){
  clusters = pam(data.frame(scale(curve.data[-which(which.zero[,ss]==1),2:7,ss])),k)$clustering
  sil.data = silhouette(pam(scale(curve.data[-which(which.zero[,ss]==1),2:7,ss]),k))
  sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
  row.names(sil.data) = sort(row.names(sil.data))
  sil.max = max(sil.data[,3])
  alpha.sil = pmax(0, sil.data[,3] / sil.max)
  
  jpeg(filename=paste('../../output/3_elucidation/pairwise_vars_muni_k0',k,'_sim.jpeg',sep=''),width=6.5,height=6,units='in',res=500)
  
  layout(matrix(1:49,7,7,byrow=T))
  par(oma=c(8.5,10.5,0.5,0.5),mar=rep(0,4))
  
  hist(curve.data[-which(which.zero[,ss]==1),1,ss],10,xaxt='n',yaxt='n',xlab='',ylab='',main=''); box()
  mtext('Mean of fitted\nnormal cdf',2,3.3,cex=0.7,las=1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),2,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),3,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),4,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),5,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),6,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),7,ss])),cex=1.1)
  
  plot(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),2,ss],xaxt='n',xlab='',ylab='',las=1,
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('S.D. of fitted\nnormal cdf',2,3.3,cex=0.7,las=1)
  hist(curve.data[-which(which.zero[,ss]==1),2,ss],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),3,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),4,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),5,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),6,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),7,ss])),cex=1.1)
  
  plot(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),3,ss],xaxt='n',xlab='',ylab='',las=1,
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('R^2 of fitted\nnormal cdf',2,3.3,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),3,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  hist(curve.data[-which(which.zero[,ss]==1),3,ss],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),4,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),5,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),6,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),7,ss])),cex=1.1)
  
  plot(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),4,ss],xaxt='n',xlab='',ylab='',las=1,
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Weeks btwn.\n5% quantile\nof empirical\nand fitted cdf',2,3.3,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),4,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  plot(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),4,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  hist(curve.data[-which(which.zero[,ss]==1),4,ss],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),4,ss],curve.data[-which(which.zero[,ss]==1),5,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),4,ss],curve.data[-which(which.zero[,ss]==1),6,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),4,ss],curve.data[-which(which.zero[,ss]==1),7,ss])),cex=1.1)
  
  plot(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),5,ss],xaxt='n',xlab='',ylab='',las=1,
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Weeks btwn.\n95% quantile\nof empirical\nand fitted cdf',2,3.3,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),5,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  plot(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),5,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  plot(curve.data[-which(which.zero[,ss]==1),4,ss],curve.data[-which(which.zero[,ss]==1),5,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  hist(curve.data[-which(which.zero[,ss]==1),5,ss],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),5,ss],curve.data[-which(which.zero[,ss]==1),6,ss])),cex=1.1)
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),5,ss],curve.data[-which(which.zero[,ss]==1),7,ss])),cex=1.1)
  
  plot(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),6,ss],xaxt='n',xlab='',ylab='',las=1,
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Weeks btwn.\nfirst and last\ncase report',2,3.3,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),6,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  plot(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),6,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  plot(curve.data[-which(which.zero[,ss]==1),4,ss],curve.data[-which(which.zero[,ss]==1),6,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  plot(curve.data[-which(which.zero[,ss]==1),5,ss],curve.data[-which(which.zero[,ss]==1),6,ss],xaxt='n',yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  hist(curve.data[-which(which.zero[,ss]==1),6,ss],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
  plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
  text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which(which.zero[,ss]==1),6,ss],curve.data[-which(which.zero[,ss]==1),7,ss])),cex=1.1)
  
  plot(curve.data[-which(which.zero[,ss]==1),1,ss],curve.data[-which(which.zero[,ss]==1),7,ss],xlab='',ylab='',las=1,
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Propn. of weeks\nbetween first and\nlast case report\nw/ no case report',2,3.3,cex=0.7,las=1)
  mtext('Mean of\nfitted\nnormal cdf',1,4.2,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),2,ss],curve.data[-which(which.zero[,ss]==1),7,ss],yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('S.D. of\nfitted\nnormal cdf',1,4.2,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),3,ss],curve.data[-which(which.zero[,ss]==1),7,ss],yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('R^2 of\nfitted\nnormal cdf',1,4.2,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),4,ss],curve.data[-which(which.zero[,ss]==1),7,ss],yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Weeks btwn.\n5% quantile\nof empirical\nand fitted cdf',1,5.25,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),5,ss],curve.data[-which(which.zero[,ss]==1),7,ss],yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Weeks btwn.\n95% quantile\nof empirical\nand fitted cdf',1,5.25,cex=0.7,las=1)
  plot(curve.data[-which(which.zero[,ss]==1),6,ss],curve.data[-which(which.zero[,ss]==1),7,ss],yaxt='n',xlab='',ylab='',
       col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
  mtext('Weeks btwn.\nfirst and last\ncase report',1,4.2,cex=0.7,las=1)
  hist(curve.data[-which(which.zero[,ss]==1),7,ss],10,yaxt='n',xlab='',ylab='',main=''); box()
  mtext('Propn. of\nweeks btwn.\nfirst and\nlast case\nreport w/ no\ncase report',1,7.4,cex=0.7,las=1)
  
  dev.off()
}



# silhouette plots for each set of groupings
ss = 1
asw = pamk(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),2:10)$crit

jpeg('../../output/3_elucidation/silhouette_muni_sim.jpeg',width=6.5,height=6.5,units='in',res=400)

layout(matrix(1:9,3,3,byrow=T))
par(mar=c(0.5,2.5,3,0.5),oma=c(0.25,2,0.25,0.25))

for(k in 2:10){
  pam.data = pam(data.frame(scale(curve.data[-which(which.zero[,ss]==1),2:7,ss])),k)
  sil.data = silhouette(pam.data)
  
  barplot(sil.data[,3],col=gg_color_hue(k,sil.data[,1],rep(1,nrow(sil.data))),ylim=c(-0.3,1),space=0,border=NA,xaxt='n',las=1)
  abline(h=0)
  
  mtext(paste('Number of groups =',k),3,line=1.3,cex=0.7)
  mtext(paste('Avg. silhouette width =',sprintf("%0.3f",asw[k])),3,line=0.1,cex=0.7)
  
  if(k == 5){
    mtext('Silhouette width',2,line=3,cex=0.7)
  }
}

dev.off()



# plot centered, normalized cumulative incidence curves for each set of groupings
jpeg(paste('../../output/3_elucidation/pca_curves_weighted_muni_sim.jpeg',sep=''),width=6.5,height=5.5,units='in',res=500)

layout(matrix(1:25,5,5,byrow=T),widths=c(1,1,0.1,1,1))
par(mar=rep(0.2,4),oma=c(3.5,3.6,2.5,1.5))

k = 2
for(ss in 1:5)
{
  clusters = pam(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss]),2)$clustering
  sil.data = silhouette(pam(curve.data[-which(which.zero[,ss]==1),2:7,ss],2))
  sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
  row.names(sil.data) = sort(row.names(sil.data))
  sil.max = max(sil.data[,3])
  alpha.sil = pmax(0, sil.data[,3] / sil.max)

  plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
  polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
  abline(v=seq(-52,52,13),col='gray')
  abline(h=seq(0,1,0.2),col='gray')
  for(jj in (1:nrow(x))[-which(which.zero[,ss]==1)][which(clusters==2)]){
    cases.norm.cum = cumsum(sims.mat[jj,,ss])/sum(sims.mat[jj,,ss])
    lines(
      -curve.data[jj,1,ss]+(1:121),
      cases.norm.cum,
      type='l',col=gg_color_hue(k,2,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
  }
  if(ss == 5){
    axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
  }
  axis(2,cex.axis=0.6,at=seq(0,1,0.2),las=1)
  if(ss == 3){
    mtext('Cumulative proportion of infections',2,line=2.5,cex=0.7)
  }
  if(ss == 1){
    mtext(expression('Small '*F[Delta*t]*' and '*F[SD]),3,line=0.1,cex=0.7)
    mtext('Classified by curve features',3,line=1.45,cex=0.7,at=58)
  }
  
  plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
  polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
  abline(v=seq(-52,52,13),col='gray')
  abline(h=seq(0,1,0.2),col='gray')
  for(jj in (1:nrow(x))[-which(which.zero[,ss]==1)][which(clusters==1)]){
    cases.norm.cum = cumsum(sims.mat[jj,,ss])/sum(sims.mat[jj,,ss])
    lines(
      -curve.data[jj,1,ss]+(1:121),
      cases.norm.cum,
      type='l',col=gg_color_hue(k,1,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
  }
  if(ss == 5){
    axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
  }
  if(ss == 1){
    mtext(expression('Large '*F[Delta*t]*' and '*F[SD]),3,line=0.1,cex=0.7)
  }
  
  plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i',bty='n')
  if(ss == 5){
    mtext('Week relative to mean of fitted normal',1,line=2.25,cex=0.7)
  }
  
  plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
  polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
  abline(v=seq(-52,52,13),col='gray')
  abline(h=seq(0,1,0.2),col='gray')
  for(jj in (1:nrow(x))[which.zero[,ss]==0 & x$R0<1]){
    cases.norm.cum = cumsum(sims.mat[jj,,ss])/sum(sims.mat[jj,,ss])
    lines(
      -curve.data[jj,1,ss]+(1:121),
      cases.norm.cum,
      type='l',col=gg_color_hue(k,2,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
  }
  if(ss == 5){
    axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
  }
  if(ss == 1){
    mtext(expression(R[0]<1),3,line=0.1,cex=0.7)
    mtext(expression('Classified by '*R[0]),3,line=1.45,cex=0.7,at=58)
  }
  
  plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
  polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
  abline(v=seq(-52,52,13),col='gray')
  abline(h=seq(0,1,0.2),col='gray')
  for(jj in (1:nrow(x))[which.zero[,ss]==0 & x$R0>=1]){
    cases.norm.cum = cumsum(sims.mat[jj,,ss])/sum(sims.mat[jj,,ss])
    lines(
      -curve.data[jj,1,ss]+(1:121),
      cases.norm.cum,
      type='l',col=gg_color_hue(k,1,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
  }
  mtext(paste('Simulated dataset ',ss,sep=''),4,line=0.25,cex=0.5)
  if(ss == 5){
    axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
  }
  if(ss == 1){
    mtext(expression(R[0]>=1),3,line=0.25,cex=0.7)
  }
}
dev.off()



# make maps of R0 above or below 1 and of the two groups that municipalities were assigned to based on simulated data
colors = rep('',length(map.muni))

jpeg(filename=paste('../../output/3_elucidation/cluster_map_muni_sim.jpeg',sep=''),width=6.5,height=6.5,units='in',res=500)
  layout(matrix(1:6,2,3,byrow=T))
  par(oma=rep(0,4),mar=rep(0,4))

  colors = ifelse(x$R0 >= 1, gg_color_hue(2,1,alpha=1), gg_color_hue(2,2,alpha=1))
  plot(map.muni,col=colors,border=NA,xlim=c(-79.4,-66.7),ylim=c(-4.5,12.6))
  plot(map.dept,border=1,add=T)
  mtext(expression({R[0]>=1} * ' (red) or ' * {R[0]<1} * ' (blue)'),3,cex=0.7,line=-1.5)
  
  for(ss in 1:5){
    clusters = pam(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),k)$clustering
    sil.data = silhouette(pam(scale(data.frame(curve.data[-which(which.zero[,ss]==1),2:7,ss])),k))
    sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
    row.names(sil.data) = sort(row.names(sil.data))
    sil.max = max(sil.data[,3])
    alpha.sil = pmax(0, sil.data[,3] / sil.max)
    colors[which(which.zero[,ss]==0)] = gg_color_hue(k,clusters,alpha=alpha.sil)
    colors[which(which.zero[,ss]==1)] = '#000000'
    plot(map.muni,col=colors,border=NA,xlim=c(-79.4,-66.7),ylim=c(-4.5,12.6))
    plot(map.dept,border=1,add=T)
    mtext(paste('Simulated dataset ',ss,sep=''),3,cex=0.7,line=-1.5)
  }
dev.off()
