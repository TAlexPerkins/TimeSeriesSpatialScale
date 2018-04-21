library(EpiEstim)

print(load('../../data/1_description/data_municipalities_new_oct7.Rdata'))

print(load("../../data/1_description/departamental_data_ZIKV_apr28_2017.RData"))
cases.dep<-cases
tot.dep<-apply(cases.dep, 2, sum)


print(load("../../data/1_description/case_matrix_ZIKV_2.RData"))

print(load("../../data/1_description/data_codes.Rdata"))
codes.munip = formatC(dat.codes$code, width=5,format='d',flag='0')
codes.dep<-sapply(codes.munip, function(x) substr(x, 1, 2))
unique.dep<-unique(codes.dep)
dep.names<-sapply(unique.dep, function(x) dat.codes$departamento[which(codes.dep==x)[1]])

order.dep<-order(tot.dep, decreasing = T)
names(tot.dep[order.dep[1:4]])


other.dep<-apply(cases.dep[, order.dep[6:33]], 1, sum)
cases.top<-cbind(cases.dep[, order.dep[1:5]], other.dep)
barplot(as.matrix(t(cases.top)), col=brewer.pal(n = 6, "Set2"), border=NA)

legend("topleft", legend=c("Valle", "N Santander", "Santander", "Tolima", "Huila", "Other"), pch=15, col=brewer.pal(6, "Set2"), bty="n")

for( i in 1:length(names(cases.dep))) {

  pdf(height=5, width=5, file=paste("../../output/1_description/figures_departments/", tolower(dep.names[which(unique.dep==names(cases.dep)[i])]), ".pdf", sep=""))
  
  print(dep.names[which(unique.dep==names(cases.dep)[i])])
  par(mar=c(4,4,4,4))
  b1<-barplot(colombia, main="", width=1, border=NA, yaxt="n")
  axis(4)
  par(new=T)
  p1<-EstimateR(colombia, T.Start=seq(which(colombia>0)[3], min(52, tail(which(colombia>0), 1)), 1), T.End=seq(which(colombia>0)[3], min(52, tail(which(colombia>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
  plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="", xaxt="n", ylab="", col="yellow")
  
  polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("yellow", alpha.f = .5), border=NA) 
  axis(2)
  mtext("R(t)", line=2.5, side = 2)
  mtext("Weekly number of cases", line=2.5, side = 4)
  
  axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))
  abline(h= c(1, 2, 5), lty=3, col="grey80")
  dev.off()
}





week.dates<-seq.Date(as.Date("2015-08-15"), by=7, length=58)
dates.plot<-seq.Date(as.Date("2015-09-01"), by="month", length.out=13)[c(1, 3, 5, 7, 9, 11, 13)]
diff.days<-dates.plot-as.Date("2015-08-15")
diff.max<-tail(week.dates, 1)-head(week.dates,1)
unit.week<-58/as.numeric(diff.max)*1.2
week.multiplier<-as.numeric(diff.days)*as.numeric(unit.week)


#### Santander
for(i in c(1,2,4:27, 29:32)) {

  pdf(height=5, width=5, file=paste("../../output/1_description/figures_departments/", tolower(dep.names[which(unique.dep==names(cases.dep)[i])]), "2.pdf", sep=""))

  print(dep.names[which(unique.dep==names(cases.dep)[i])])
  par(mar=c(4,4,4,4))
  
  b1<-barplot(cases.dep[,i], main="", width=1, border=NA, yaxt="n")
  axis(4)
  par(new=T)
  p1<-EstimateR(cases.dep[,i], T.Start=seq(which(cases.dep[,i]>0)[3], min(52, tail(which(cases.dep[,i]>0), 1)), 1), T.End=seq(which(cases.dep[,i]>0)[3], min(52, tail(which(cases.dep[,i]>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
  plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="", xaxt="n", ylab="", col="royalblue")
  
  polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("royalblue", alpha.f = .5), border=NA) 
  axis(2)
  mtext("R(t)", line=2.5, side = 2)
  mtext("Weekly number of cases", line=2.5, side = 4)
  
  axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))
  abline(h= c(1, 2, 5), lty=3, col="grey80")
  dev.off()
}


for(j in c(1,2,4:27, 29:32)) {
  
  pdf(height=5, width=5, file=paste("../../output/1_description/figures_municipalities/", tolower(dep.names[which(unique.dep==names(cases.dep)[i])]), ".pdf", sep=""))
  
  codes.cases<-sapply(dimnames(cases)[[2]], function(x) substr(x, 1, 2))
  mat.santander<-cases[,which(codes.cases==names(cases.dep)[i])]
  tot.santander<-apply(mat.santander, 2, sum)
  codes.santander<-dimnames(cases)[[2]][which(codes.cases==names(cases.dep)[i])]
  mat.include<-mat.santander[, order(tot.santander, decreasing = T)[1:4]]
  codes.include<- dimnames(mat.include)[[2]]
  
  if(sum(tot.santander)>100) {
    
  par(mfrow=c(2,2), oma=c(4,4,4,4), mar=c(.5, .5, .5, .5))
  
  for(i in 1:4) {
    b1<-barplot(mat.include[,i], main="", width=1, border=NA, xaxt="n", yaxt="n", ylim=c(0, max(mat.include[,1])))
    if(is.element(i, c(2,4))) {axis(4)}
    if(is.element(i, c(1,2))) {axis(1, tick = T, labels = F)}
    
    par(new=T)
    p1<-EstimateR(mat.include[,i], T.Start=seq(which(mat.include[,i]>0)[3], min(52, tail(which(mat.include[,i]>0), 1)), 1), T.End=seq(which(mat.include[,i]>0)[3], min(52, tail(which(mat.include[,i]>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
    
    plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="Date", xaxt="n", ylab="Weekly number of cases", col="orangered")
    
    polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("orangered", alpha.f = .5), border=NA) 
    if(is.element(i, c(1,3))) {axis(2)}
    if(is.element(i, c(3,4))) {axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))}
    if(is.element(i, c(1,2, 3, 4))) {abline(h= c(1, 2, 5), lty=3, col="grey80")}
  }
    mtext("Weekly number of cases", side = 4, outer=T, line=2)
    mtext("R(t)", side = 2, outer=T, line=2)  
  }
  
  dev.off()
}


#### Santander
for(i in c(1,2,4:27, 29:32)) {
  pdf(height=5, width=5, file=paste("../../output/1_description/figures_departments/", tolower(dep.names[which(unique.dep==names(cases.dep)[i])]), "2.pdf", sep=""))
  
  print(dep.names[which(unique.dep==names(cases.dep)[i])])
  par(mar=c(4,4,4,4))
  
  codes.cases<-sapply(dimnames(cases)[[2]], function(x) substr(x, 1, 2))
  mat.santander<-cases[,which(codes.cases==names(cases.dep)[i])]
  tot.santander<-apply(mat.santander, 2, sum)
  codes.santander<-dimnames(cases)[[2]][which(codes.cases==names(cases.dep)[i])]
  mat.include<-mat.santander[, order(tot.santander, decreasing = T)[1:5]]
  vec.rest<-apply(mat.santander[, order(tot.santander, decreasing = T)[6:ncol(mat.santander)]], 1, sum)
  mat.include<-cbind(mat.include, vec.rest)
  codes.include<- dimnames(mat.include)[[2]]
  
  
  b1<-barplot(t(mat.include), main="", width=1, border=NA, yaxt="n", col=adjustcolor(brewer.pal(6, "Dark2"), alpha.f = .6))
  axis(4)
  par(new=T)
  p1<-EstimateR(cases.dep[,i], T.Start=seq(which(cases.dep[,i]>0)[3], min(52, tail(which(cases.dep[,i]>0), 1)), 1), T.End=seq(which(cases.dep[,i]>0)[3], min(52, tail(which(cases.dep[,i]>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
  plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="", xaxt="n", ylab="", col="black")
  
  polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("grey30", alpha.f = .5), border=NA) 
  axis(2)
  mtext("R(t)", line=2.5, side = 2)
  mtext("Weekly number of cases", line=2.5, side = 4)
  axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))
  abline(h= c(1, 2, 5), lty=3, col="grey80")
  dev.off()
}


### Municipalities 2
for(j in c(2,4:27, 29:32)) {
  
  pdf(height=5, width=5, file=paste("../../output/1_description/figures_municipalities/", tolower(dep.names[which(unique.dep==names(cases.dep)[j])]), "2.pdf", sep=""))
  
  codes.cases<-sapply(dimnames(cases)[[2]], function(x) substr(x, 1, 2))
  mat.santander<-cases[,which(codes.cases==names(cases.dep)[j])]
  tot.santander<-apply(mat.santander, 2, sum)
  codes.santander<-dimnames(cases)[[2]][which(codes.cases==names(cases.dep)[j])]
  mat.include<-mat.santander[, order(tot.santander, decreasing = T)[1:4]]
  codes.include<- dimnames(mat.include)[[2]]
  print(dat.codes[match(codes.include, dat.codes$code),])
  if(sum(tot.santander)>100) {
    
    par(mfrow=c(2,2), oma=c(4,4,4,4), mar=c(.5, .5, .5, .5))
    
    for(i in 1:4) {
      b1<-barplot(mat.include[,i], main="", width=1, border=NA, xaxt="n", yaxt="n", ylim=c(0, max(mat.include[,1])), col=adjustcolor(brewer.pal(6, "Dark2")[i], .6))
      if(is.element(i, c(2,4))) {axis(4)}
      if(is.element(i, c(1,2))) {axis(1, tick = T, labels = F)}
      
      par(new=T)
      p1<-EstimateR(mat.include[,i], T.Start=seq(which(mat.include[,i]>0)[3], min(52, tail(which(mat.include[,i]>0), 1)), 1), T.End=seq(which(mat.include[,i]>0)[3], min(52, tail(which(mat.include[,i]>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
      
      plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="Date", xaxt="n", ylab="Weekly number of cases", col="black")
      
      polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("grey30", alpha.f = .5), border=NA) 
      if(is.element(i, c(1,3))) {axis(2)}
      if(is.element(i, c(3,4))) {axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))}
      if(is.element(i, c(1,2, 3, 4))) {abline(h= c(1, 2, 5), lty=3, col="grey80")}
    }
    mtext("Weekly number of cases", side = 4, outer=T, line=2)
    mtext("R(t)", side = 2, outer=T, line=2)  
  }
  
  dev.off()
}


## Huila, Quindio, Risaralda, Santander, Cundinamarca

dep.interest<-names(dep.names[which(is.element(tolower(dep.names), c("huila", "quindio", "risaralda", "santander")))])
ind.interest<-match(dep.interest, names(cases.dep))

for(j in ind.interest) {
  codes.cases<-sapply(dimnames(cases)[[2]], function(x) substr(x, 1, 2))
  mat.santander<-cases[,which(codes.cases==names(cases.dep)[j])]
  tot.santander<-apply(mat.santander, 2, sum)
  codes.santander<-dimnames(cases)[[2]][which(codes.cases==names(cases.dep)[j])]
  mat.include<-mat.santander[, order(tot.santander, decreasing = T)[1:4]]
  codes.include<- dimnames(mat.include)[[2]]
  print(dat.codes[match(codes.include, dat.codes$code),])
}


## Timing of epidemics in the top 5 departments
barplot(cases.top[,1], col=adjustcolor(brewer.pal(6, "Set2")[1], alpha.f = .5))
barplot(cases.top[,2], col=adjustcolor(brewer.pal(6, "Set2")[2], alpha.f = .5), add=T)
barplot(cases.top[,3], col=adjustcolor(brewer.pal(6, "Set2")[3], alpha.f = .5), add=T)
barplot(cases.top[,4], col=adjustcolor(brewer.pal(6, "Set2")[4], alpha.f = .5), add=T)





for(i in 1:length(names(cases.dep))) {
	par(ask=T)
	if(sum(cases.dep[,i])>10) {
	
	EstimateR(cases.dep[,i], T.Start=seq(which(cases.dep[,i]>0)[3], min(52, tail(which(cases.dep[,i]>0), 1)), 1), T.End=seq(which(cases.dep[,i]>0)[3], min(52, tail(which(cases.dep[,i]>0), 1)), 1)+6, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=T)
		print(dep.names[which(unique.dep==names(cases.dep)[i])])
	}
}

all.country<-apply(cases, 1, sum)
mun.tots<-apply(cases, 2, sum)
order.mun<-order(mun.tots, decreasing=T)


p1<-EstimateR(all.country, T.Start=seq(7, 56,1), T.End=seq(9, 58, 1), Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=T)


for(i in 1:20) {
	par(ask=T)
	EstimateR(cases[, order.mun[i]], T.Start=seq(which(cases[, order.mun[i]]>0)[1], min(54, tail(which(cases[, order.mun[i]]>0), 1)), 1), T.End=seq(which(cases[, order.mun[i]]>0)[1], min(54, tail(which(cases[, order.mun[i]]>0), 1)), 1)+4, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=T)
}


### Country plot
#save(dep.data, file="../data/Colombia_Processed/departamental_data.RData")

colombia<-cases.top
col.vec<-apply(colombia, 1, sum)

pdf(height=5, width=5, file=paste("../../output/1_description/figures_departments/","Colombia", "2.pdf", sep=""))

  print(dep.names[which(unique.dep==names(cases.dep)[i])])
  par(mar=c(4,4,4,4))
  
  b1<-barplot(t(colombia), main="", width=1, border=NA, yaxt="n", col=adjustcolor(brewer.pal(6, "Dark2"), alpha.f = .6))
  axis(4)
  par(new=T)
  p1<-EstimateR(col.vec, T.Start=seq(which(col.vec>0)[3], min(52, tail(which(col.vec>0), 1)), 1), T.End=seq(which(col.vec>0)[3], min(52, tail(which(col.vec>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
  plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="", xaxt="n", ylab="", col="black")
  
  polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("grey30", alpha.f = .5), border=NA) 
  axis(2)
  mtext("R(t)", line=2.5, side = 2)
  mtext("Weekly number of cases", line=2.5, side = 4)
  axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))
  abline(h= c(1, 2, 5), lty=3, col="grey80")
dev.off()


### Colombia_deps_2
### Municipalities 2

## Departments are Valle, Norte de Santander, Santander and Tolima
dep.names[which(is.element(names(dep.names), names(tot.dep[order.dep[1:4]])))][c(4, 1, 2, 3)]


  pdf(height=5, width=5, file=paste("../../output/1_description/figures_municipalities/", "Colombia", "2.pdf", sep=""))
  
    par(mfrow=c(2,2), oma=c(4,4,4,4), mar=c(.5, .5, .5, .5))
    
    for(i in 1:4) {
      b1<-barplot(colombia[,i], main="", width=1, border=NA, xaxt="n", yaxt="n", ylim=c(0, max(colombia[,1])), col=adjustcolor(brewer.pal(6, "Dark2")[i], .6))
      if(is.element(i, c(2,4))) {axis(4)}
      if(is.element(i, c(1,2))) {axis(1, tick = T, labels = F)}
      
      par(new=T)
      p1<-EstimateR(colombia[,i], T.Start=seq(which(colombia[,i]>0)[3], min(52, tail(which(colombia[,i]>0), 1)), 1), T.End=seq(which(colombia[,i]>0)[3], min(52, tail(which(colombia[,i]>0), 1)), 1)+5, Mean.SI=15/7, Std.SI=3/7, method=c("ParametricSI"), plot=F)
      
      plot(b1[p1$R$T.End], p1$R$`Mean(R)`, ylim=c(0.5, 10), yaxt="n", type="l", xlim=c(0, 70), log="y", xlab="Date", xaxt="n", ylab="Weekly number of cases", col="black")
      
      polygon(c(b1[p1$R$T.End], rev(b1[p1$R$T.End])), c(p1$R$`Quantile.0.025(R)`, rev(p1$R$`Quantile.0.975(R)`)), col=adjustcolor("grey30", alpha.f = .5), border=NA) 
      if(is.element(i, c(1,3))) {axis(2)}
      if(is.element(i, c(3,4))) {axis(1, at=week.multiplier, labels=c("Sept", "Nov", "Jan", "Mar", "May", "Jul", "Sept"))}
      if(is.element(i, c(1,2, 3, 4))) {abline(h= c(1, 2, 5), lty=3, col="grey80")}
    }
    mtext("Weekly number of cases", side = 4, outer=T, line=2)
    mtext("R(t)", side = 2, outer=T, line=2)  

  dev.off()
