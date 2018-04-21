# load necessary R packages
if(!require(cluster)){install.packages('cluster');library(cluster)}
if(!require(fpc)){install.packages('fpc');library(fpc)}
if(!require(ggfortify)){install.packages('ggfortify');library(ggfortify)}
if(!require(maptools)){install.packages('maptools');library(maptools)}
if(!require(RColorBrewer)){install.packages('RColorBrewer');library(RColorBrewer)}
library(getcartr)



# function for ggplot colors
gg_color_hue <- function(n, k, alpha.in) {
  hues = seq(15, 375, length = n + 1)
  sapply(1:length(k), function(ii)
    hcl(h = hues, l = 65, c = 100, alpha = alpha.in[ii])[1:n][k[ii]])
}



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
case.cums.muni = 5 * case.cums.muni
case.cums.muni[case.cums.muni==0] = 1
carto.cases.muni = quick.carto(map.muni,case.cums.muni)
carto.cases.muni.tmp = gBuffer(carto.cases.muni, byid=TRUE, width=0)
carto.cases.dept = unionSpatialPolygons(carto.cases.muni.tmp,carto.cases.muni.tmp@data$COD_DEPTO)



# fit S curves to case report data
curve.data = matrix(0,length(munip.list),7)
which.zero = rep(0,length(munip.list))
for(ii in 1:length(munip.list)){
  if(sum(munip.list[[ii]]$cases_combo) > 0){
    cases.norm = munip.list[[ii]]$cases_combo / sum(munip.list[[ii]]$cases_combo)
    cases.norm.cum = cumsum(munip.list[[ii]]$cases_combo) / sum(munip.list[[ii]]$cases_combo)
    tmp = optim(c(30,1),function(par){
      obs = cases.norm.cum
      pred = pnorm(1:length(cases.norm.cum),par[1],par[2])
      R2 = 1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2)
      return(-R2)})
    curve.data[ii,1] = tmp$par[1]
    curve.data[ii,2] = tmp$par[2]
    curve.data[ii,3] = -tmp$value
    curve.data[ii,4] =
      max(which(0.05 >
        cases.norm.cum)) -
      qnorm(0.05,tmp$par[1],tmp$par[2])
    curve.data[ii,5] =
      min(which(0.95 <
        cases.norm.cum)) -
      qnorm(0.95,tmp$par[1],tmp$par[2])
    curve.data[ii,6] =
      length(min(which(munip.list[[ii]]$cases_combo > 0)):max(which(munip.list[[ii]]$cases_combo > 0)))
    curve.data[ii,7] =
      1 - sum(cases.norm[
        min(which(munip.list[[ii]]$cases_combo > 0)):max(which(munip.list[[ii]]$cases_combo > 0))] > 0) / curve.data[ii,6]
  }
  if(sum(munip.list[[ii]]$cases_combo) == 0){
    which.zero[ii] = 1
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
which.zero = which(which.zero==1)
save(curve.data,which.zero,file='../../output/2_classification/muni/curve_data.RData')



# plot pairwise combinations of curve features stratified by grouping
for(k in 2:6){
  clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
  sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
  sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
  row.names(sil.data) = sort(row.names(sil.data))
  sil.max = max(sil.data[,3])
  alpha.sil = pmax(0, sil.data[,3] / sil.max)

  jpeg(filename=paste('../../output/2_classification/muni/pairwise_vars_muni_k0',k,'.jpeg',sep=''),width=6.5,height=6,units='in',res=500)
  
    layout(matrix(1:49,7,7,byrow=T))
    par(oma=c(8.5,10.5,0.5,0.5),mar=rep(0,4))
    
    hist(curve.data[-which.zero,1],10,xaxt='n',yaxt='n',xlab='',ylab='',main=''); box()
    mtext('Mean of fitted\nnormal cdf',2,3.3,cex=0.7,las=1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,1],curve.data[-which.zero,2])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,1],curve.data[-which.zero,3])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,1],curve.data[-which.zero,4])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,1],curve.data[-which.zero,5])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,1],curve.data[-which.zero,6])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,1],curve.data[-which.zero,7])),cex=1.1)
    
    plot(curve.data[-which.zero,1],curve.data[-which.zero,2],xaxt='n',xlab='',ylab='',las=1,
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('S.D. of fitted\nnormal cdf',2,3.3,cex=0.7,las=1)
    hist(curve.data[-which.zero,2],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,2],curve.data[-which.zero,3])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,2],curve.data[-which.zero,4])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,2],curve.data[-which.zero,5])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,2],curve.data[-which.zero,6])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,2],curve.data[-which.zero,7])),cex=1.1)
    
    plot(curve.data[-which.zero,1],curve.data[-which.zero,3],xaxt='n',xlab='',ylab='',las=1,
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('R^2 of fitted\nnormal cdf',2,3.3,cex=0.7,las=1)
    plot(curve.data[-which.zero,2],curve.data[-which.zero,3],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    hist(curve.data[-which.zero,3],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,3],curve.data[-which.zero,4])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,3],curve.data[-which.zero,5])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,3],curve.data[-which.zero,6])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,3],curve.data[-which.zero,7])),cex=1.1)
    
    plot(curve.data[-which.zero,1],curve.data[-which.zero,4],xaxt='n',xlab='',ylab='',las=1,
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Weeks btwn.\n5% quantile\nof empirical\nand fitted cdf',2,3.3,cex=0.7,las=1)
    plot(curve.data[-which.zero,2],curve.data[-which.zero,4],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    plot(curve.data[-which.zero,3],curve.data[-which.zero,4],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    hist(curve.data[-which.zero,4],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,4],curve.data[-which.zero,5])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,4],curve.data[-which.zero,6])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,4],curve.data[-which.zero,7])),cex=1.1)
    
    plot(curve.data[-which.zero,1],curve.data[-which.zero,5],xaxt='n',xlab='',ylab='',las=1,
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Weeks btwn.\n95% quantile\nof empirical\nand fitted cdf',2,3.3,cex=0.7,las=1)
    plot(curve.data[-which.zero,2],curve.data[-which.zero,5],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    plot(curve.data[-which.zero,3],curve.data[-which.zero,5],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    plot(curve.data[-which.zero,4],curve.data[-which.zero,5],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    hist(curve.data[-which.zero,5],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,5],curve.data[-which.zero,6])),cex=1.1)
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,5],curve.data[-which.zero,7])),cex=1.1)
    
    plot(curve.data[-which.zero,1],curve.data[-which.zero,6],xaxt='n',xlab='',ylab='',las=1,
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Weeks btwn.\nfirst and last\ncase report',2,3.3,cex=0.7,las=1)
    plot(curve.data[-which.zero,2],curve.data[-which.zero,6],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    plot(curve.data[-which.zero,3],curve.data[-which.zero,6],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    plot(curve.data[-which.zero,4],curve.data[-which.zero,6],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    plot(curve.data[-which.zero,5],curve.data[-which.zero,6],xaxt='n',yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    hist(curve.data[-which.zero,6],10,xaxt='n',yaxt='n',xlab='',ylab='',main='')
    plot(-100,-100,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
    text(0.5,0.5,sprintf("%.2f",cor(curve.data[-which.zero,6],curve.data[-which.zero,7])),cex=1.1)
    
    plot(curve.data[-which.zero,1],curve.data[-which.zero,7],xlab='',ylab='',las=1,
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Propn. of weeks\nbetween first and\nlast case report\nw/ no case report',2,3.3,cex=0.7,las=1)
    mtext('Mean of\nfitted\nnormal cdf',1,4.2,cex=0.7,las=1)
    plot(curve.data[-which.zero,2],curve.data[-which.zero,7],yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('S.D. of\nfitted\nnormal cdf',1,4.2,cex=0.7,las=1)
    plot(curve.data[-which.zero,3],curve.data[-which.zero,7],yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('R^2 of\nfitted\nnormal cdf',1,4.2,cex=0.7,las=1)
    plot(curve.data[-which.zero,4],curve.data[-which.zero,7],yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Weeks btwn.\n5% quantile\nof empirical\nand fitted cdf',1,5.25,cex=0.7,las=1)
    plot(curve.data[-which.zero,5],curve.data[-which.zero,7],yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Weeks btwn.\n95% quantile\nof empirical\nand fitted cdf',1,5.25,cex=0.7,las=1)
    plot(curve.data[-which.zero,6],curve.data[-which.zero,7],yaxt='n',xlab='',ylab='',
         col=gg_color_hue(k,clusters,alpha=alpha.sil^2),pch=19,cex=0.5)
    mtext('Weeks btwn.\nfirst and last\ncase report',1,4.2,cex=0.7,las=1)
    hist(curve.data[-which.zero,7],10,yaxt='n',xlab='',ylab='',main=''); box()
    mtext('Propn. of\nweeks btwn.\nfirst and\nlast case\nreport w/ no\ncase report',1,7.4,cex=0.7,las=1)
    
  dev.off()
}



# principal components analysis on curve features
pca.data = prcomp(curve.data[-which.zero,c(2:7)],center=T,scale=T)
curve.data = cbind(curve.data,matrix(0,nrow(curve.data),6))
curve.data[-which.zero,8:13] = predict(pca.data,newdata=curve.data[-which.zero,2:7])
jpeg(filename='../../output/2_classification/muni/pca_variance_muni.jpeg',width=6.5,height=3.25,units='in',res=500)
  par(oma=rep(0,4),mar=c(2,4,0.5,0))
  barplot(summary(pca.data)$importance[3,],las=1,ylab='Cumulative proportion of variance explained')
dev.off()



# make maps with municipalities color coded by group and weighted by silhouette width
for(k in 2:4){
  colors = rep('',nrow(curve.data))
  clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
  sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
  sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
  row.names(sil.data) = sort(row.names(sil.data))
  sil.max = max(sil.data[,3])
  alpha.sil = pmax(0, sil.data[,3] / sil.max)
  colors[(1:length(colors))[-which.zero]] = gg_color_hue(k,clusters,alpha=alpha.sil)
  colors[which.zero] = '#000000'
  
  jpeg(filename=paste('../../output/2_classification/muni/cluster_map_muni_k0',k,'.jpeg',sep=''),width=6.5,height=3,units='in',res=500)
    layout(matrix(1:3,1,3))
    par(oma=rep(0,4),mar=rep(0,4))
    plot(map.muni,col=colors,border=NA,xlim=c(-79.4,-66.7),ylim=c(-4.5,12.6))
    plot(map.dept,border=1,add=T)
    plot(carto.pop.muni,col=colors,border=NA)
    plot(carto.pop.dept,border=1,add=T)
    plot(carto.cases.muni,col=colors,border=NA)
    plot(carto.cases.dept,border=1,add=T)
  dev.off()
}



# plot centered, normalized cumulative incidence curves for 2, 3, and 4 groups
jpeg(paste('../../output/2_classification/muni/pca_curves_weighted_muni.jpeg',sep=''),width=6.5,height=4.25,units='in',res=500)

layout(matrix(c(1,2,10,11,3,4,5,12,6,7,8,9),3,4,byrow=T))
par(mar=rep(0.2,4),oma=c(3.5,3.6,1.75,2.5))

k = 2
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
row.names(sil.data) = sort(row.names(sil.data))
sil.max = max(sil.data[,3])
alpha.sil = pmax(0, sil.data[,3] / sil.max)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==1)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,1,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
axis(2,cex.axis=0.6,at=seq(0,1,0.2),las=1)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==2)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,2,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
mtext('Municipal-level standardized cumulative incidence curves',3,line=0.5,cex=0.7,at=58)

k = 3
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
row.names(sil.data) = sort(row.names(sil.data))
sil.max = max(sil.data[,3])
alpha.sil = pmax(0, sil.data[,3] / sil.max)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==1)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,1,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
axis(2,cex.axis=0.6,at=seq(0,1,0.2),las=1)
mtext('Cumulative proportion of infections',2,line=2.5,cex=0.7)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==2)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,2,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==3)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,3,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}

k = 4
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
row.names(sil.data) = sort(row.names(sil.data))
sil.max = max(sil.data[,3])
alpha.sil = pmax(0, sil.data[,3] / sil.max)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==1)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,1,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
axis(2,cex.axis=0.6,at=seq(0,1,0.2),las=1)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==2)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,2,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
mtext('Week relative to mean of fitted normal',1,line=2.25,cex=0.7,at=58)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==3)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,3,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
for(jj in (1:length(munip.list))[-which.zero][which(clusters==4)]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,4,alpha=alpha.sil[jj]),lwd=0.5*alpha.sil[jj])
}
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
mtext('Number of groups = 4',4,line=1,cex=0.7)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i',bty='n')
plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i',bty='n')
mtext('Number of groups = 2',4,line=1,cex=0.7)

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i',bty='n')
mtext('Number of groups = 3',4,line=1,cex=0.7)

dev.off()



# plot the medoid centered, normalized cumulative incidence curves for each of 2 and 3 groups
jpeg(paste('../../output/2_classification/muni/pca_curves_medoid_muni.jpeg',sep=''),width=6.5,height=2.5,units='in',res=500)

layout(matrix(1:3,1,3,byrow=T))
par(mar=rep(0.2,4),oma=c(3.5,3.6,1.5,1.5))

k = 2
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
medoids = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$id.med
sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
row.names(sil.data) = sort(row.names(sil.data))
sil.max = max(sil.data[,3])

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
cluster = 1
for(jj in (1:length(munip.list))[-which.zero][medoids]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,cluster,alpha=1),lwd=2)
  cluster = cluster + 1
}
axis(2,cex.axis=0.6,at=seq(0,1,0.2),las=1)
mtext('Cumulative proportion of infections',2,line=2.5,cex=0.7)
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
mtext('Number of groups = 2',3,line=0.5,cex=0.7)

k = 3
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
medoids = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$id.med
sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
row.names(sil.data) = sort(row.names(sil.data))
sil.max = max(sil.data[,3])

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
cluster = 1
for(jj in (1:length(munip.list))[-which.zero][medoids]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,cluster,alpha=1),lwd=2)
  cluster = cluster + 1
}
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
mtext('Number of groups = 3',3,line=0.5,cex=0.7)
mtext('Week relative to mean of fitted normal',1,line=2.25,cex=0.7)

k = 4
clusters = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$clustering
medoids = pam(scale(data.frame(curve.data[-which.zero,2:7])),k)$id.med
sil.data = silhouette(pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k))
sil.data[as.numeric(row.names(sil.data[,])),] = sil.data
row.names(sil.data) = sort(row.names(sil.data))
sil.max = max(sil.data[,3])

plot(-1000,xlim=c(-52,52),ylim=c(0,1),xaxt='n',yaxt='n',yaxs='i');box()
polygon(c(-1e3,1e3,1e3,-1e3),c(-1e3,-1e3,1e3,1e3),col=rgb(0,0,0,0.6))
abline(v=seq(-52,52,13),col='gray')
abline(h=seq(0,1,0.2),col='gray')
cluster = 1
for(jj in (1:length(munip.list))[-which.zero][medoids]){
  cases.norm.cum = cumsum(munip.list[[jj]]$cases_combo) / sum(munip.list[[jj]]$cases_combo)
  lines(
    -curve.data[jj,1]+(1:58),
    cases.norm.cum,
    type='l',col=gg_color_hue(k,cluster,alpha=1),lwd=2)
  cluster = cluster + 1
}
axis(1,cex.axis=0.6,at=seq(-52,52,13),label=c('-52','','-26','','0','','26','','52'))
mtext('Number of groups = 4',3,line=0.5,cex=0.7)

dev.off()



# silhouette plots for each set of groupings
asw = pamk(scale(data.frame(curve.data[-which.zero,c(2:7)])),2:10)$crit

jpeg('../../output/2_classification/muni/silhouette_muni.jpeg',width=6.5,height=6.5,units='in',res=400)
  
  layout(matrix(1:9,3,3,byrow=T))
  par(mar=c(0.5,2.5,3,0.5),oma=c(0.25,2,0.25,0.25))
  
  for(k in 2:10){
    pam.data = pam(scale(data.frame(curve.data[-which.zero,c(2:7)])),k)
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
