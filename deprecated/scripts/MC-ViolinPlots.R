###################################################
### Create violin plots for Monte Carlo results
rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)

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
              results_tb[!is.na(results_tb$rho_se),# drops non-converged spmles
                         c("N", "method", "rho", "beta0_hat", "beta1_hat", "rho_hat")])
}
rm(files, f, results_tb)
unique(paste(sp$method, sp$N, sep="-"))

### Import Temporal MC results
files <- paste("output/temporal/", dir("output/temporal"), sep="")
te <- data.frame()
for(f in files) {
  print(f)
  load(f)
  te <- rbind(te, results_tb[!is.na(results_tb$gamma_se) & !is.na(results_tb$rho_se),# drops non-converged spmles
                             c("N", "TT", "method", "gamma", "beta0_hat", "beta1_hat", "gamma_hat")]
              )
}
rm(files, f, results_tb)
unique(paste0(te$method, te$N, te$TT))

### Import Spatio-Temporal MC results
files <- paste("output/spatio-temporal/", dir("output/spatio-temporal"), sep="")
st <- data.frame()
for(f in files) {
  if(f!="mcresults-ST-49n20t-RIS-100r-200212.Rdata"){
    print(f)
    load(f)
    st <- rbind(st, results_tb[!is.na(results_tb$rho_se) & !is.na(results_tb$gamma_se), # drops non-converged spmles
                               c("N", "TT", "method", "rho", "gamma", "beta0_hat", "beta1_hat", "rho_hat", "gamma_hat")])
  }

}
rm(files, f, results_tb)
unique(paste0(st$method, st$N, st$TT))


###################################################
### Plot spatial results
## Prep
sp$Estimator <- as.factor(sp$method)
levels(sp$Estimator) <- list(Bayes="bayes", GMM="gmm", MLE="naiveprobit", RIS="ris", SPMLE="spmle")

spplots <- expand.grid(rho=unique(sp$rho),
                       N=unique(sp$N))
spplots <- spplots[order(spplots$N, spplots$rho),]
spplots$index <- 1:dim(spplots)[1]

plotsout <- list()

## Plot
for(i in spplots$index){
  plotsout[[i]] <- ggplot(subset(sp, rho==spplots$rho[i] & N==spplots$N[i])) +
    geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
    scale_y_continuous(limits = c(-1, 1), name=expression("Distribution of " *rho* " estimates")) +
    geom_hline(yintercept=spplots$rho[i]) +
    ggtitle(bquote(paste("N=", .(spplots[i,2]), ", ", rho, "=", .(spplots[i,1])))) +
    theme_minimal()
}

pdf("plots/MC_sp-ViolinPlots-200228.pdf")
  grid.arrange(plotsout[[1]], plotsout[[2]], plotsout[[3]],
               plotsout[[4]], plotsout[[5]], plotsout[[6]],
               plotsout[[7]], plotsout[[8]], plotsout[[9]],
               ncol=3, as.table = FALSE)
dev.off()

rm(plotsout, spplots)

###################################################
### Plot temporal results
## Prep
te$Estimator <- as.factor(te$method)
levels(te$Estimator) <- list(RIS="ris", SPMLE="spmle")

teplots <- expand.grid(gamma=unique(te$gamma),
                       TT=unique(te$TT))
teplots <- teplots[order(teplots$T, teplots$gamma),]
teplots$index <- 1:dim(teplots)[1]

plotsout <- list()

## Plot
for(i in teplots$index){
  plotsout[[i]] <- ggplot(subset(te, gamma==teplots$gamma[i] & TT==teplots$TT[i])) +
                          geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                          scale_y_continuous(limits = c(-1, 1), name=expression("Distribution of " *gamma* " estimates")) +
                          geom_hline(yintercept=teplots$gamma[i]) +
                          ggtitle(bquote(paste("N=64, T=", .(teplots[i,2]), ", ", gamma, "=", .(teplots[i,1])))) +
                          theme_minimal()
}

pdf("plots/MC_te-ViolinPlots-200228.pdf")
  grid.arrange(plotsout[[1]], plotsout[[2]], plotsout[[3]],
               plotsout[[4]], plotsout[[5]], plotsout[[6]],
               ncol=2, as.table = FALSE)
dev.off()

rm(plotsout, teplots)

###################################################
### Plot spatio-temporal results (gamma)
## Prep
st$Estimator <- as.factor(st$method)
levels(st$Estimator) <- list(RIS="ris", SPMLE="spmle")

stplots <- expand.grid(gamma=unique(st$gamma),
                       rho = unique(st$rho),
                       TT=unique(st$TT))
stplots <- stplots[!(stplots$rho==.5 & stplots$gamma==.5),]
stplots <- stplots[order(stplots$TT, stplots$gamma, stplots$rho),]
stplots$index <- 1:dim(stplots)[1]

plotsout <- list()

## Plot
for(i in stplots$index){
  plotsout[[i]] <- ggplot(subset(st, gamma==stplots$gamma[i] & rho==stplots$rho[i] & TT==stplots$TT[i])) +
    geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
    scale_y_continuous(limits = c(-1, 1), name=expression("Distribution of " *gamma* " estimates")) +
    geom_hline(yintercept=stplots$gamma[i]) +
    ggtitle(bquote(paste("N=64, T=", .(stplots[i,3]), ", ", gamma, "=", .(stplots[i,1]), ", ", rho, "=", .(stplots[i,2])))) +
    theme_minimal()
}

pdf("plots/MC_st-ViolinPlots-gamma-200228.pdf")
  grid.arrange(plotsout[[1]], plotsout[[2]], plotsout[[3]],
               plotsout[[4]], plotsout[[5]], plotsout[[6]],
               ncol=2, as.table = FALSE)
dev.off()

rm(plotsout)


###################################################
### Plot spatio-temporal results (rho)
## Prep
stplots <- stplots[order(stplots$TT, stplots$rho, stplots$gamma),]
stplots$index <- 1:dim(stplots)[1]

plotsout <- list()

## Plot
for(i in stplots$index){
  plotsout[[i]] <- ggplot(subset(st, gamma==stplots$gamma[i] & rho==stplots$rho[i] & TT==stplots$TT[i])) +
    geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
    scale_y_continuous(limits = c(-1, 1), name=expression("Distribution of " *rho* " estimates")) +
    geom_hline(yintercept=stplots$rho[i]) +
    ggtitle(bquote(paste("N=64, T=", .(stplots[i,3]), ", ", gamma, "=", .(stplots[i,1]), ", ", rho, "=", .(stplots[i,2])))) +
    theme_minimal()
}

pdf("plots/MC_st-ViolinPlots-rho-200228.pdf")
  grid.arrange(plotsout[[1]], plotsout[[2]], plotsout[[3]],
               plotsout[[4]], plotsout[[5]], plotsout[[6]],
               ncol=2, as.table = FALSE)
dev.off()

rm(plotsout, stplots)



