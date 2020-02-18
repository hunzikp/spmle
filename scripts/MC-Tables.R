###################################################
### Create violin plots for Monte Carlo results
rm(list=ls())
library(xtable)
library(dplyr); library(tidyr)

if(Sys.info()["user"]=="ncb") setwd("/home/ncb/Dropbox/SpatialLogit/MonteCarlo/")
getwd()

###################################################
### Import Spatial MC results
## N=256
sp <- read.csv("output/spatial/mcresults-256n-500r-200127.csv", header=T, sep=",")
## N=1024
load("output/spatial/mcresults-sp-1024n-500r-200127.Rdata")
sp <- rbind(sp,results_tb)
rm(results_tb)
## N=4096
load("output/spatial/mcresults-sp-4096n-500r-200128.Rdata")
sp <- rbind(sp,results_tb)
rm(results_tb)
## N=16384
load("output/spatial/mcresults-sp-16384n-250r-200128.Rdata")
sp <- rbind(sp,results_tb)
rm(results_tb)

###################################################
### Compute summary stats of interest

biastab_sp <- sp %>%
                group_by(rho, N, method) %>%
                summarise(meanbias=mean(abs(rho_hat-rho), na.rm=T),
                          sdestim=sd(rho_hat,na.rm=T),
                          rmse=sqrt(mean((rho_hat-rho)^2, na.rm=T)),
                          nfail = sum(is.na(time))
                          )
print(biastab_sp, n=100)
biastab_sp <- biastab_sp[biastab_sp$rho!=.75,]

# bias statistics to long format
biastab_sp_long <- gather(biastab_sp, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_sp_long

# MC experiments to wide-format
biastab_sp_long$rhoN <- paste(biastab_sp_long$rho, biastab_sp_long$N, sep="-")

tab <- matrix(NA_real_, nrow=20, 12)

for(i in 1:length(unique(biastab_sp_long$rhoN))){
  temp <- biastab_sp_long[biastab_sp_long$rhoN==unique(biastab_sp_long$rhoN)[i], c("method", "result")]
  if("bayes"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="bayes", "result"])
  if("gmm"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="gmm", "result"])
  if("naiveprobit"%in%temp$method) tab[9:12,i] <- as.matrix(temp[temp$method=="naiveprobit", "result"])
  if("ris"%in%temp$method) tab[13:16,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[17:20,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "SD", "RMSE", "No convergence"),as.data.frame(tab)), digits = 3)


###################################################
### Import Spatio-temporal MC results
## N=256
sp <- read.csv("output/spatial/mcresults-256n-500r-200127.csv", header=T, sep=",")
## N=1024
load("output/spatial/mcresults-sp-1024n-500r-200127.Rdata")
sp <- rbind(sp,results_tb)
rm(results_tb)
## N=4096
load("output/spatial/mcresults-sp-4096n-500r-200128.Rdata")
sp <- rbind(sp,results_tb)
rm(results_tb)
## N=16384
load("output/spatial/mcresults-sp-16384n-250r-200128.Rdata")
sp <- rbind(sp,results_tb)
rm(results_tb)
