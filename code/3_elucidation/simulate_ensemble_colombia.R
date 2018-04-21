# clear workspace
rm(list=ls())

# load required packages
if(!require(pomp)){install.packages('pomp');library(pomp)}

set.seed(126)

# load and process data
x = read.csv('../../data/3_elucidation/R0_municip_NatMicro.csv')
which.girardot = which(x$NOM_MUNICI=='GIRARDOT')
x$R0 = x$R0 * 4.61 / x$R0[which.girardot] # http://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2016.21.28.30283
x$R0[is.na(x$R0)] = 0
pop = read.csv('../../data/3_elucidation/Population_municipalities_2016.csv')
x$pop = rep(0,nrow(x))
x$pop[-which(x$ID_ESPACIA==setdiff(x$ID_ESPACIA,pop$Municipality))] =
  pop$Population[match(pop$Municipality,x$ID_ESPACIA)]
x[which(is.na(x$pop)),]$pop = 3463 # https://es.wikipedia.org/wiki/Concepción_(Antioquia)
load('../../data/2_classification/muni/data_municipalities_new_nov5.Rdata')
x$cases.data = 0
for(id in as.character(x$ID_ESPACIA)){
  if(sum(names(munip.list)==id)>0){
    x$cases.data[which(as.character(x$ID_ESPACIA)==id)] =
      sum(munip.list[[which(names(munip.list)==id)[1]]]$cases_combo)
  }
}
x$cum.infected.data = 1 - (x$pop - x$cases.data) / x$pop

# importation pattern (32.575968, 8.858822)
optim(par=c(35,10),fn=function(par){
  sum((cumsum(rowSums(sapply(munip.list,function(mm)mm$cases_combo)))/sum(sapply(munip.list,function(mm)mm$cases_combo))-pnorm(1:length(munip.list[[1]]$cases_combo),par[1],par[2]))^2)})$par

# define code that governs stochastic events in single time step
sir_step =
  Csnippet("
     float K = K_0 * N * (1 + delta * cos(2 * 3.141593 * (t + theta)));
     float epsilon = epsilon_0 / (1 + delta * cos(2 * 3.141593 * (t + theta)));
     float eta = eta_0 / (1 - delta * cos(2 * 3.141593 * (t + theta)));
     float cIsum = cI_055 + cI_060 + cI_065 + cI_070 + cI_075 + cI_080 + cI_085 + cI_090 + cI_095 + cI_100 + cI_105 + cI_110 + cI_115;
     float psi = kappa * beta_hm / N * cIsum;

     float lambda_local = kappa * beta_mh * Y / (A + H1 + H2 + Y + 1);
     float lambda_import_now = lambda_import * dnorm(t/7, 32.575968, 8.858822, 0);

     R_0 = pow(kappa,2) * beta_mh * beta_hm * 6 / (epsilon * (1 + epsilon / eta)) * (A + H1 + H2 + Y) / N;
     R_e = R_0 * (S / N);

     incidence_local = rbinom(S,1-exp(-lambda_local*dt));
     incidence_import = rbinom(S-incidence_local,1-exp(-lambda_import_now*dt));
     if(t == 1){
       incidence_import += 10;
     }

     int dcI_000 = incidence_local + incidence_import;
     int dcI_005 = rbinom(cI_000,dt/0.5);
     int dcI_010 = rbinom(cI_005,dt/0.5);
     int dcI_015 = rbinom(cI_010,dt/0.5);
     int dcI_020 = rbinom(cI_015,dt/0.5);
     int dcI_025 = rbinom(cI_020,dt/0.5);
     int dcI_030 = rbinom(cI_025,dt/0.5);
     int dcI_035 = rbinom(cI_030,dt/0.5);
     int dcI_040 = rbinom(cI_035,dt/0.5);
     int dcI_045 = rbinom(cI_040,dt/0.5);
     int dcI_050 = rbinom(cI_045,dt/0.5);
     int dcI_055 = rbinom(cI_050,dt/0.5);
     int dcI_060 = rbinom(cI_055,dt/0.5);
     int dcI_065 = rbinom(cI_060,dt/0.5);
     int dcI_070 = rbinom(cI_065,dt/0.5);
     int dcI_075 = rbinom(cI_070,dt/0.5);
     int dcI_080 = rbinom(cI_075,dt/0.5);
     int dcI_085 = rbinom(cI_080,dt/0.5);
     int dcI_090 = rbinom(cI_085,dt/0.5);
     int dcI_095 = rbinom(cI_090,dt/0.5);
     int dcI_100 = rbinom(cI_095,dt/0.5);
     int dcI_105 = rbinom(cI_100,dt/0.5);
     int dcI_110 = rbinom(cI_105,dt/0.5);
     int dcI_115 = rbinom(cI_110,dt/0.5);
     int dcI_120 = rbinom(cI_115,dt/0.5);
     
     S -= dcI_000;
     cI_000 += dcI_000 - dcI_005;
     cI_005 += dcI_005 - dcI_010;
     cI_010 += dcI_010 - dcI_015;
     cI_015 += dcI_015 - dcI_020;
     cI_020 += dcI_020 - dcI_025;
     cI_025 += dcI_025 - dcI_030;
     cI_030 += dcI_030 - dcI_035;
     cI_035 += dcI_035 - dcI_040;
     cI_040 += dcI_040 - dcI_045;
     cI_045 += dcI_045 - dcI_050;
     cI_050 += dcI_050 - dcI_055;
     cI_055 += dcI_055 - dcI_060;
     cI_060 += dcI_060 - dcI_065;
     cI_065 += dcI_065 - dcI_070;
     cI_070 += dcI_070 - dcI_075;
     cI_075 += dcI_075 - dcI_080;
     cI_080 += dcI_080 - dcI_085;
     cI_085 += dcI_085 - dcI_090;
     cI_090 += dcI_090 - dcI_095;
     cI_095 += dcI_095 - dcI_100;
     cI_100 += dcI_100 - dcI_105;
     cI_105 += dcI_105 - dcI_110;
     cI_110 += dcI_110 - dcI_115;
     cI_115 += dcI_115 - dcI_120;
     H += dcI_000;
     
     double dLin = rpois(b*(A+H1+H2+Y)*dt);
     double dLout[2];
     double ratesLout[2];
     ratesLout[0] = alpha;
     ratesLout[1] = omega * (1 + L / K);
     reulermultinom(2,L,&ratesLout[0],dt,&dLout[0]);
     
     double dAout[2];
     double ratesAout[2];
     ratesAout[0] = psi;
     ratesAout[1] = epsilon;
     reulermultinom(2,A,&ratesAout[0],dt,&dAout[0]);
     
     double dH1out[2], dH2out[2];
     double ratesHout[2];
     ratesHout[0] = 2 * eta;
     ratesHout[1] = epsilon;
     reulermultinom(2,H1,&ratesHout[0],dt,&dH1out[0]);
     reulermultinom(2,H2,&ratesHout[0],dt,&dH2out[0]);
     
     L += dLin - dLout[0] - dLout[1];
     A += dLout[0] - dAout[0] - dAout[1];
     H1 += dAout[0] - dH1out[0] - dH1out[1];
     H2 += dH1out[0] - dH2out[0] - dH2out[1];
     Y += dH2out[0] - rbinom(Y,1-exp(-epsilon*dt));
     ")

# measurement probability (not really important but required)
dmeas =
  Csnippet(
    "lik = dbinom(C,cI_000,rho_mean,give_log);
    ")

# measurement process (not really important but required)
rmeas =
  Csnippet(
    "C = rbinom(cI_000,rho_mean);
    ")

# some data (not really important but required)
num.weeks = 120
data = data.frame(
  C = rep(0,num.weeks+1),
  t = seq(1,num.weeks*7+1,7)
)

# define vectors of key parameters
N.in.vec = pmax(x$pop,sort(x$pop)[2])
m_eq.in.vec = x$R0 / (0.5 ^ 2 * 0.7 * 0.7 * 6 / (0.2 * (1 + 0.2 / (1/8.4))))

# create a data frame with all combinations of the key parameters
df = data.frame(
  N = N.in.vec,
  m_eq = m_eq.in.vec
)

# create list to store simulation outputs in
sims.list = list()

# perform simulations across different parameter combinations
for(ii in 1:nrow(df)){
  # extract values of key parameters being manipulated
  N.in = df$N[ii]
  m_eq = max(df$m_eq[ii],min(df$m_eq[!is.na(df$m_eq)][df$m_eq[!is.na(df$m_eq)]>0]))
  if(is.na(m_eq)){
    m_eq = min(df$m_eq[!is.na(df$m_eq)][df$m_eq[!is.na(df$m_eq)]>0])
  }
  
  # define values of other key parameters
  lambda_import.in = 1.55 * 10 ^ -3
  rho.mean.in = 0.004 # https://elifesciences.org/articles/29820
  rho.mean.in = 0.115 # http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726
  rho.cv.in = 1
  epsilon.in = 0.1
  epsilon_0.in = 0.2
  Rm = 2.69

  # calculate derived parameters
  b.in = Rm * epsilon_0.in * (1/19 + 0.025) / (1/19)
  K_0.in = m_eq / (1/19 / epsilon_0.in * (1/19 * (b.in - epsilon.in) / (epsilon_0.in * 0.025) - 1))
  L_eq.in = ((1/19 + 0.025) - sqrt((1/19 + 0.025) ^ 2 + 4 * (0.025 / K_0.in) * b.in * m_eq)) / (-2 * (0.025 / K_0.in))
  R_0.in = 0.5 ^ 2 * 0.7 * 0.7 * 6 / (epsilon_0.in * (1 + epsilon_0.in / (1/8.4))) * m_eq
  
  # define initial conditions for the state variables
  # (note that infections are seeded by importation)
  sir_init =
    Csnippet(paste("
       S = N;
       incidence_local = 0;
       incidence_import = 0;
       cI_000 = 0;
       cI_005 = 0;
       cI_010 = 0;
       cI_015 = 0;
       cI_020 = 0;
       cI_025 = 0;
       cI_030 = 0;
       cI_035 = 0;
       cI_040 = 0;
       cI_045 = 0;
       cI_050 = 0;
       cI_055 = 0;
       cI_060 = 0;
       cI_065 = 0;
       cI_070 = 0;
       cI_075 = 0;
       cI_080 = 0;
       cI_085 = 0;
       cI_090 = 0;
       cI_095 = 0;
       cI_100 = 0;
       cI_105 = 0;
       cI_110 = 0;
       cI_115 = 0;
       L = ", round(0.86 * L_eq.in * N.in), ";
       A = ", round(0.86 * m_eq * N.in), ";
       H1 = 0;
       H2 = 0;
       Y = 0;
       H = 0;
       R_0 = ", R_0.in, ";
       R_e = ", R_0.in, ";
       ", sep=''))
  
  # define pomp object
  sir = pomp(
    data=data,
    times='t',
    t0=0,
    rprocess=euler.sim(sir_step,delta.t=1/6),
    rmeasure=rmeas,
    dmeasure=dmeas,
    initializer=sir_init,
    paramnames=c("K_0","N","delta","theta","epsilon_0","eta_0","kappa","beta_hm","lambda_import","beta_mh","b","alpha","omega","rho_mean"),
    statenames=c("S","incidence_local","incidence_import","cI_000","cI_005","cI_010","cI_015","cI_020","cI_025","cI_030","cI_035","cI_040","cI_045","cI_050","cI_055","cI_060","cI_065","cI_070","cI_075","cI_080","cI_085","cI_090","cI_095","cI_100","cI_105","cI_110","cI_115","L","A","H1","H2","Y","H","R_0","R_e"))
  
  # perform simulations
  nsim = 100
  sims.list[[ii]] =
    simulate(
      sir,
      params=c(K_0=K_0.in,N=N.in,delta=0.05,theta=0.25,epsilon_0=epsilon_0.in,eta_0=1/8.4,kappa=0.5,beta_hm=0.7,lambda_import=lambda_import.in,beta_mh=0.7,b=b.in,alpha=1/19,omega=0.025,rho_mean=rho.mean.in),
      nsim=nsim,
      as.data.frame=TRUE,
      include.data=TRUE)
  print(ii)
}



# save sims.list for further analysis
save(data,x,sims.list,file='../../output/3_elucidation/simulations.RData')



# check how reasonable the simulation outputs are
x$cum.infected.sim = sapply(sims.list,function(ss)(ss$S[num.weeks+2]-ss$S[2*num.weeks+2])/ss$S[num.weeks+2])
x$cases.sim = sapply(sims.list,function(ss)sum(ss$C[(num.weeks+2):(2*num.weeks+2)]))


