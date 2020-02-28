###################################################
### Export Monte Carlo Summary Statistics
rm(list=ls())
library(xtable)
library(dplyr); library(tidyr)

if(Sys.info()["user"]=="ncb") setwd("/home/ncb/Dropbox/SpatialLogit/MonteCarlo/")
getwd()

###################################################
### Import Spatial MC results
files <- paste("output/spatial/", dir("output/spatial"), sep="")
files
files <- files[c(1:3,5:6)]
sp <- data.frame()
for(f in files) {
  print(f)
  load(f)
  sp <- rbind(sp,
              results_tb[,
                         c("N", "method", "rho", "beta0", "beta1", "time",
                           "beta0_hat", "beta1_hat", "rho_hat",
                           "beta0_se", "beta1_se", "rho_se")]
              )
}
rm(files, f, results_tb)
unique(paste(sp$method, sp$N, sep="-"))

# set time and estimates to NA if error is missing
sp[is.na(sp$rho_se), c("time", "beta0_hat", "beta1_hat", "beta1_se", "rho_hat")] <- NA

### Import Temporal MC results
files <- paste("output/temporal/", dir("output/temporal"), sep="")
te <- data.frame()
for(f in files) {
  print(f)
  load(f)
  te <- rbind(te,
              results_tb[,
                         c("N", "TT", "method", "gamma", "beta0", "beta1", "time",
                           "beta0_hat", "beta1_hat", "gamma_hat",
                           "beta0_se", "beta1_se", "gamma_se")]
            )
}
rm(files, f, results_tb)
unique(paste0(te$method, te$N, te$TT))

# set time and estimates to NA if error is missing
te[is.na(te$gamma_se), c("time", "beta0_hat", "beta1_hat", "beta1_se", "gamma_hat")] <- NA


### Import Spatio-Temporal MC results
files <- paste("output/spatio-temporal/", dir("output/spatio-temporal"), sep="")
files <- files[1:4]
st <- data.frame()
for(f in files) {
    print(f)
    load(f)
    st <- rbind(st,
                results_tb[, # drops non-converged spmles
                           c("N", "TT", "method", "gamma", "rho", "beta0", "beta1", "time",
                             "beta0_hat", "beta1_hat", "gamma_hat", "rho_hat",
                             "beta0_se", "beta1_se", "gamma_se", "rho_se")]

                )
}

rm(files, f, results_tb)
unique(paste0(st$method, st$N, st$TT))

# set time and estimates to NA if error is missing
st[is.na(st$gamma_se) | is.na(st$rho_se),
   c("time", "beta0_hat", "beta1_hat", "beta1_se", "gamma_hat", "gamma_se", "rho_hat", "rho_se")] <- NA

###################################################
### Compute and export SPATIAL summary stats
## rho
biastab_sp <- sp %>%
                group_by(rho, N, method) %>%
                summarise(meanbias=mean(abs(rho_hat-rho), na.rm=T),
                          #sdestim=sd(rho_hat,na.rm=T),
                          rmse=sqrt(mean((rho_hat-rho)^2, na.rm=T)),
                          oc = sd(rho_hat,na.rm=T)/mean(rho_se, na.rm=T),  # overconfidence
                          nfail = sum(is.na(time))
                          )
print(biastab_sp, n=100)

# bias statistics to long format
biastab_sp_long <- gather(biastab_sp, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_sp_long

# MC experiments to wide-format
biastab_sp_long$rhoN <- paste(biastab_sp_long$rho, biastab_sp_long$N, sep="-")

tab <- matrix(NA_real_, nrow=20, 9)

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

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)


############
## beta_0
biastab_sp <- sp %>%
  group_by(rho, N, method) %>%
  summarise(meanbias=mean(abs(beta0_hat-beta0), na.rm=T),
            #sdestim=sd(rho_hat,na.rm=T),
            rmse=sqrt(mean((beta0_hat-beta0)^2, na.rm=T)),
            oc = sd(beta0_hat,na.rm=T)/mean(beta0_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
print(biastab_sp, n=100)

# bias statistics to long format
biastab_sp_long <- gather(biastab_sp, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_sp_long

# MC experiments to wide-format
biastab_sp_long$rhoN <- paste(biastab_sp_long$rho, biastab_sp_long$N, sep="-")

tab <- matrix(NA_real_, nrow=20, 9)

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

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)



############
## beta_1
biastab_sp <- sp %>%
  group_by(rho, N, method) %>%
  summarise(meanbias=mean(abs(beta1_hat-beta1), na.rm=T),
            #sdestim=sd(beta1_hat,na.rm=T),
            rmse=sqrt(mean((beta1_hat-beta1)^2, na.rm=T)),
            oc = sd(beta1_hat,na.rm=T)/mean(beta1_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
print(biastab_sp, n=100)

# bias statistics to long format
biastab_sp_long <- gather(biastab_sp, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_sp_long

# MC experiments to wide-format
biastab_sp_long$rhoN <- paste(biastab_sp_long$rho, biastab_sp_long$N, sep="-")

tab <- matrix(NA_real_, nrow=20, 9)

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

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)


###################################################
### Compute and export TEMPORAL summary stats
############
## gamma
biastab_te <- te %>%
  group_by(gamma, N, TT, method) %>%
  summarise(meanbias=mean(abs(gamma_hat-gamma), na.rm=T),
            #sdestim=sd(rho_hat,na.rm=T),
            rmse=sqrt(mean((gamma_hat-gamma)^2, na.rm=T)),
            oc = sd(gamma_hat,na.rm=T)/mean(gamma_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
biastab_te

# bias statistics to long format
biastab_te_long <- gather(biastab_te, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_te_long

# MC experiments to wide-format
biastab_te_long$gammaN <- paste(biastab_te_long$gamma, biastab_te_long$N, biastab_te_long$TT, sep="-")

tab <- matrix(NA_real_, nrow=8, 9)

for(i in 1:length(unique(biastab_te_long$gammaN))){
  temp <- biastab_te_long[biastab_te_long$gammaN==unique(biastab_te_long$gammaN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)


############
## beta_0
biastab_te <- te %>%
  group_by(gamma, N, TT, method) %>%
  summarise(meanbias=mean(abs(beta0_hat-beta0), na.rm=T),
            #sdestim=sd(rho_hat,na.rm=T),
            rmse=sqrt(mean((beta0_hat-beta0)^2, na.rm=T)),
            oc = sd(beta0_hat,na.rm=T)/mean(beta0_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
biastab_te

# bias statistics to long format
biastab_te_long <- gather(biastab_te, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_te_long

# MC experiments to wide-format
biastab_te_long$gammaN <- paste(biastab_te_long$gamma, biastab_te_long$N, biastab_te_long$TT, sep="-")

tab <- matrix(NA_real_, nrow=8, 9)

for(i in 1:length(unique(biastab_te_long$gammaN))){
  temp <- biastab_te_long[biastab_te_long$gammaN==unique(biastab_te_long$gammaN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)

############
## beta_1
biastab_te <- te %>%
  group_by(gamma, N, TT, method) %>%
  summarise(meanbias=mean(abs(beta1_hat-beta1), na.rm=T),
            #sdestim=sd(rho_hat,na.rm=T),
            rmse=sqrt(mean((beta1_hat-beta1)^2, na.rm=T)),
            oc = sd(beta1_hat,na.rm=T)/mean(beta1_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
biastab_te

# bias statistics to long format
biastab_te_long <- gather(biastab_te, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_te_long

# MC experiments to wide-format
biastab_te_long$gammaN <- paste(biastab_te_long$gamma, biastab_te_long$N, biastab_te_long$TT, sep="-")

tab <- matrix(NA_real_, nrow=8, 9)

for(i in 1:length(unique(biastab_te_long$gammaN))){
  temp <- biastab_te_long[biastab_te_long$gammaN==unique(biastab_te_long$gammaN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)


###################################################
### Compute and export SPATIO-TEMPORAL summary stats

############
## RHO
biastab_st <- st %>%
  group_by(gamma, rho, N, TT, method) %>%
  summarise(meanbias=mean(abs(rho_hat-rho), na.rm=T),
            #sdestim=sd(rho_hat,na.rm=T),
            rmse=sqrt(mean((rho_hat-rho)^2, na.rm=T)),
            oc = sd(rho_hat,na.rm=T)/mean(rho_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
print(biastab_st, n=100)

# bias statistics to long format
biastab_st_long <- gather(biastab_st, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_st_long

# MC experiments to wide-format
biastab_st_long$rhoN <- paste(biastab_st_long$gamma, biastab_st_long$rho, biastab_st_long$N, biastab_st_long$TT, sep="-")

tab <- matrix(NA_real_, nrow=8, 9)

for(i in 1:length(unique(biastab_st_long$rhoN))){
  temp <- biastab_st_long[biastab_st_long$rhoN==unique(biastab_st_long$rhoN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)

############
## GAMMA
biastab_st <- st %>%
  group_by(gamma, rho, N, TT, method) %>%
  summarise(meanbias=mean(abs(gamma_hat-gamma), na.rm=T),
            #sdestim=sd(rho_hat,na.rm=T),
            rmse=sqrt(mean((gamma_hat-gamma)^2, na.rm=T)),
            oc = sd(gamma_hat,na.rm=T)/mean(gamma_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
biastab_st

# bias statistics to long format
biastab_st_long <- gather(biastab_st, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_st_long

# MC experiments to wide-format
biastab_st_long$gammaN <- paste(biastab_st_long$gamma, biastab_st_long$rho, biastab_st_long$N, biastab_st_long$TT, sep="-")

tab <- matrix(NA_real_, nrow=8, 9)

for(i in 1:length(unique(biastab_st_long$gammaN))){
  temp <- biastab_st_long[biastab_st_long$gammaN==unique(biastab_st_long$gammaN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)



############
## BETA_0
biastab_st <- st %>%
  group_by(gamma, rho, N, TT, method) %>%
  summarise(meanbias=mean(abs(beta0_hat-beta0), na.rm=T),
            #sdestim=sd(beta0_hat,na.rm=T),
            rmse=sqrt(mean((beta0_hat-beta0)^2, na.rm=T)),
            oc = sd(beta0_hat,na.rm=T)/mean(beta0_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
print(biastab_st, n=100)

# bias statistics to long format
biastab_st_long <- gather(biastab_st, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_st_long

# MC experiments to wide-format
biastab_st_long$rhoN <- paste(biastab_st_long$gamma, biastab_st_long$rho, biastab_st_long$N, biastab_st_long$TT, sep="-")

tab <- matrix(NA_real_, nrow=8, 9)

for(i in 1:length(unique(biastab_st_long$rhoN))){
  temp <- biastab_st_long[biastab_st_long$rhoN==unique(biastab_st_long$rhoN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)


############
## BETA_1
biastab_st <- st %>%
  group_by(gamma, rho, N, TT, method) %>%
  summarise(meanbias=mean(abs(beta1_hat-beta1), na.rm=T),
            #sdestim=sd(beta1_hat,na.rm=T),
            rmse=sqrt(mean((beta1_hat-beta1)^2, na.rm=T)),
            oc = sd(beta1_hat,na.rm=T)/mean(beta1_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
print(biastab_st, n=100)

# bias statistics to long format
biastab_st_long <- gather(biastab_st, statistic, result, meanbias:nfail, factor_key=FALSE)
biastab_st_long

# Index
biastab_st_long$rhoN <- paste(biastab_st_long$gamma, biastab_st_long$rho, biastab_st_long$N, biastab_st_long$TT, sep="-")

# input matrix
tab <- matrix(NA_real_, nrow=8, 9)

# fill matrix
for(i in 1:length(unique(biastab_st_long$rhoN))){
  temp <- biastab_st_long[biastab_st_long$rhoN==unique(biastab_st_long$rhoN)[i], c("method", "result")]
  if("ris"%in%temp$method) tab[1:4,i] <- as.matrix(temp[temp$method=="ris", "result"])
  if("spmle"%in%temp$method) tab[5:8,i] <- as.matrix(temp[temp$method=="spmle", "result"])
  rm(temp)
  print(i)
}

# print for export
xtable(cbind(c("Mean Bias", "RMSE", "Overconfidence", "No convergence"),as.data.frame(tab)), digits = 3)


###################################################
### Compute and export SPATIO-TEMPORAL summary stats for PMLE
biastab_st <- st %>%
  filter(method=="spmle" & N==64 & TT==16) %>% # only 64 x 16
  group_by(gamma, rho, N, TT) %>%
  summarise(# beta_0
            mean_b0=mean(beta0_hat, na.rm=T),
            meanbias_b0=mean(abs(beta0_hat-beta0), na.rm=T),
            rmse_b0=sqrt(mean((beta0_hat-beta0)^2, na.rm=T)),
            sdestim_b0=sd(beta0_hat,na.rm=T),
            oc_b0 = sd(beta1_hat,na.rm=T)/mean(beta0_se, na.rm=T),  # overconfidence
            # beta_1
            mean_b1=mean(beta1_hat, na.rm=T),
            meanbias_b1=mean(abs(beta0_hat-beta0), na.rm=T),
            rmse_b1=sqrt(mean((beta1_hat-beta1)^2, na.rm=T)),
            sdestim_b1=sd(beta1_hat,na.rm=T),
            oc_b1 = sd(beta1_hat,na.rm=T)/mean(beta1_se, na.rm=T),  # overconfidence
            # gamma
            mean_gamma=mean(gamma_hat, na.rm=T),
            meanbias_gamma=mean(abs(gamma_hat-gamma), na.rm=T),
            rmse_gamma=sqrt(mean((gamma_hat-gamma)^2, na.rm=T)),
            sdestim_gamma=sd(gamma_hat,na.rm=T),
            oc_gamma = sd(gamma_hat,na.rm=T)/mean(gamma_se, na.rm=T),  # overconfidence
            # rho
            mean_rho=mean(rho_hat, na.rm=T),
            meanbias_rho=mean(abs(rho_hat-rho), na.rm=T),
            rmse_rho=sqrt(mean((rho_hat-rho)^2, na.rm=T)),
            sdestim_rho=sd(gamma_hat,na.rm=T),
            oc_rho = sd(rho_hat,na.rm=T)/mean(rho_se, na.rm=T),  # overconfidence
            nfail = sum(is.na(time))
  )
print(biastab_st, n=100)

# bias statistics to long format
biastab_st_long <- gather(biastab_st, statistic, result, mean_b0:nfail, factor_key=FALSE)
biastab_st_long

# Index
biastab_st_long$experiment <- paste(biastab_st_long$gamma, biastab_st_long$rho, biastab_st_long$N, biastab_st_long$TT, sep="-")

# input matrix
tab <- matrix(NA_real_, nrow=15, 4)

# fill matrix
for(i in 1:length(unique(biastab_st_long$experiment))){
  temp <- biastab_st_long[biastab_st_long$experiment==unique(biastab_st_long$experiment)[i],
                          c("statistic", "result")]
  tab[1:5+(5*(i-1)),1] <- as.matrix(temp$result[grepl("b0", temp$statistic)])
  tab[1:5+5*(i-1),2] <- as.matrix(temp$result[grepl("b1", temp$statistic)])
  tab[1:5+5*(i-1),3] <- as.matrix(temp$result[grepl("gamma", temp$statistic)])
  tab[1:5+5*(i-1),4] <- as.matrix(temp$result[grepl("rho", temp$statistic)])
  rm(temp)
  print(i)
}

# print for export
print(
  xtable(cbind(c("Mean Coefficient Estimate", "Mean Bias", "RMSE", "Actual SD of estimates", "Overconfidence"),
               as.data.frame(tab))
         , digits = 3),
  include.rownames = F
)

#################################################################################
### Compute and export SPATIO-TEMPORAL summary stats for PMLE (all sample sizes)
biastab_st <- st %>%
  filter(method=="spmle") %>%
  group_by(gamma, rho, N, TT) %>%
  summarise(# beta_0
    mean_b0=mean(beta0_hat, na.rm=T),
    meanbias_b0=mean(abs(beta0_hat-beta0), na.rm=T),
    rmse_b0=sqrt(mean((beta0_hat-beta0)^2, na.rm=T)),
    sdestim_b0=sd(beta0_hat,na.rm=T),
    oc_b0 = sd(beta1_hat,na.rm=T)/mean(beta0_se, na.rm=T),  # overconfidence
    # beta_1
    mean_b1=mean(beta1_hat, na.rm=T),
    meanbias_b1=mean(abs(beta0_hat-beta0), na.rm=T),
    rmse_b1=sqrt(mean((beta1_hat-beta1)^2, na.rm=T)),
    sdestim_b1=sd(beta1_hat,na.rm=T),
    oc_b1 = sd(beta1_hat,na.rm=T)/mean(beta1_se, na.rm=T),  # overconfidence
    # gamma
    mean_gamma=mean(gamma_hat, na.rm=T),
    meanbias_gamma=mean(abs(gamma_hat-gamma), na.rm=T),
    rmse_gamma=sqrt(mean((gamma_hat-gamma)^2, na.rm=T)),
    sdestim_gamma=sd(gamma_hat,na.rm=T),
    oc_gamma = sd(gamma_hat,na.rm=T)/mean(gamma_se, na.rm=T),  # overconfidence
    # rho
    mean_rho=mean(rho_hat, na.rm=T),
    meanbias_rho=mean(abs(rho_hat-rho), na.rm=T),
    rmse_rho=sqrt(mean((rho_hat-rho)^2, na.rm=T)),
    sdestim_rho=sd(gamma_hat,na.rm=T),
    oc_rho = sd(rho_hat,na.rm=T)/mean(rho_se, na.rm=T),  # overconfidence
    nfail = sum(is.na(time))
  )
print(biastab_st, n=100)

# bias statistics to long format
biastab_st_long <- gather(biastab_st, statistic, result, mean_b0:nfail, factor_key=FALSE)
biastab_st_long

# Index
biastab_st_long <- biastab_st_long[order(biastab_st_long$N,biastab_st_long$N, biastab_st_long$TT, biastab_st_long$gamma, biastab_st_long$rho),]
biastab_st_long$experiment <- paste(biastab_st_long$gamma, biastab_st_long$rho, sep="-")
biastab_st_long$samplesize <- paste(biastab_st_long$N, biastab_st_long$TT, sep="-")
unique(biastab_st_long$experiment)
unique(biastab_st_long$samplesize)

# input matrix
tab <- matrix(NA_real_, nrow=15, 12)

# fill matrix
for(i in 1:length(unique(biastab_st_long$samplesize))){
  for(j in 1:length(unique(biastab_st_long$experiment))){
    temp <- biastab_st_long[biastab_st_long$samplesize==unique(biastab_st_long$samplesize)[i] &
                            biastab_st_long$experiment==unique(biastab_st_long$experiment)[j],
                            c("statistic", "result")]
    tab[1:5+(5*(j-1)), 1+4*(i-1)] <- as.matrix(temp$result[grepl("b0", temp$statistic)])
    tab[1:5+(5*(j-1)), 2+4*(i-1)] <- as.matrix(temp$result[grepl("b1", temp$statistic)])
    tab[1:5+(5*(j-1)), 3+4*(i-1)] <- as.matrix(temp$result[grepl("gamma", temp$statistic)])
    tab[1:5+(5*(j-1)), 4+4*(i-1)] <- as.matrix(temp$result[grepl("rho", temp$statistic)])
    rm(temp)
    }
}

# print for export
print(
  xtable(cbind(c("Mean Coefficient Estimate", "Mean Bias", "RMSE", "Actual SD of estimates", "Overconfidence"),
               as.data.frame(tab))
         , digits = 3),
  include.rownames = F
)
