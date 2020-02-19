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
## N=256
sp1 <- read.csv("output/spatial/mcresults-256n-500r-200127.csv", header=T, sep=",")
## N=1024
load("output/spatial/mcresults-sp-1024n-500r-200127.Rdata")
sp2 <- results_tb
rm(results_tb)
## N=4096
load("output/spatial/mcresults-sp-4096n-500r-200128.Rdata")
sp3 <- results_tb
rm(results_tb)
## N=16384
load("output/spatial/mcresults-sp-16384n-250r-200128.Rdata")
sp4 <- results_tb
rm(results_tb)

### Import Temporal MC results
files <- paste("output/temporal/", dir("output/temporal"), sep="")
te <- data.frame()
for(f in files) {
  print(f)
  load(f)
  te <- rbind(te, results_tb[,c("N", "TT", "method", "gamma", "beta0_hat", "beta1_hat", "gamma_hat")])
}
rm(files, f, results_tb)

### Import Spatio-Temporal MC results
files <- paste("output/spatio-temporal/", dir("output/spatio-temporal"), sep="")
st <- data.frame()
for(f in files) {
  print(f)
  load(f)
  st <- rbind(st, results_tb[,c("N", "TT", "method", "gamma", "rho",
                                "beta0_hat", "beta1_hat", "gamma_hat", "rho_hat")])
}
rm(files, f, results_tb)
st <- st[st$N!=49,]


###################################################
### Plot spatial results
### N = 256
sp1_rho0 <- ggplot(subset(sp1, rho==0)) +
                     geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                     geom_hline(yintercept=0) +
                     theme_bw()

sp1_rho025 <- ggplot(subset(sp1, rho==0.25)) +
                      geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                      geom_hline(yintercept=.25) +
                      theme_bw()

sp1_rho05 <- ggplot(subset(sp1, rho==0.5)) +
                      geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                      geom_hline(yintercept=.5) +
                      theme_bw()

### N = 1024
sp2_rho0 <- ggplot(subset(sp2, rho==0)) +
            geom_violin(aes(x=as.factor(method), y=rho_hat)) +
            geom_hline(yintercept=0) +
            theme_bw()

sp2_rho025 <- ggplot(subset(sp2, rho==0.25)) +
              geom_violin(aes(x=as.factor(method), y=rho_hat)) +
              geom_hline(yintercept=.25) +
              theme_bw()

sp2_rho05 <- ggplot(subset(sp2, rho==0.5)) +
              geom_violin(aes(x=as.factor(method), y=rho_hat)) +
              geom_hline(yintercept=.5) +
              theme_bw()

### N = 4096
sp3_rho0 <- ggplot(subset(sp3, rho==0)) +
                  geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                  geom_hline(yintercept=0) +
                  theme_bw()

sp3_rho025 <- ggplot(subset(sp3, rho==0.25)) +
                    geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                    geom_hline(yintercept=.25) +
                    theme_bw()

sp3_rho05 <- ggplot(subset(sp3, rho==0.5)) +
                    geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                    geom_hline(yintercept=.5) +
                    theme_bw()

### N = 16384
sp4_rho0 <- ggplot(subset(sp4, rho==0)) +
                  geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                  geom_hline(yintercept=0) +
                  theme_bw()

sp4_rho025 <- ggplot(subset(sp4, rho==0.25)) +
                    geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                    geom_hline(yintercept=.25) +
                    theme_bw()

sp4_rho05 <- ggplot(subset(sp4, rho==0.5)) +
                    geom_violin(aes(x=as.factor(method), y=rho_hat)) +
                    geom_hline(yintercept=.5) +
                    theme_bw()


grid.arrange(sp1_rho0, sp2_rho0, sp3_rho0, sp4_rho0,
             sp1_rho025, sp2_rho025, sp3_rho025, sp4_rho025,
             sp1_rho05, sp2_rho05, sp3_rho05, sp4_rho05,
            ncol=4)



###################################################
### Plot temporal results
## Prep
te$Estimator <- as.factor(te$method)
levels(te$Estimator) <- list(RIS="ris", SPMLE="spmle")

teplots <- expand.grid(gamma=unique(te$gamma),
                       TT=unique(te$TT))
teplots <- teplots[order(teplots$T, teplots$gamma),]
teplots$index <- 1:dim(teplots_index)[1]

plotsout <- list()

## Plot
for(i in teplots$index){
  plotsout[[i]] <- ggplot(subset(te, gamma==teplots$gamma[i] & TT==teplots$TT[i])) +
                          geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                          scale_y_continuous(limits = c(-1, 1), name=expression("Distribution of " *gamma* " estimates")) +
                          geom_hline(yintercept=teplots$gamma[i]) +
                          ggtitle(paste0("N=64, T=", teplots$TT[i], ", ", expression(gamma), "=", teplots$gamma[i])) +
                          theme_minimal()
}

pdf("plots/MC_te-ViolinPlots-200219.pdf")
  grid.arrange(plotsout[[1]], plotsout[[2]], plotsout[[3]],
               plotsout[[4]], plotsout[[5]], plotsout[[6]],
               ncol=2, as.table = FALSE)
dev.off()

rm(plotsout, teplots)

###################################################
### Plot spatio-temporal results
## Prep
st$Estimator <- as.factor(st$method)
levels(st$Estimator) <- list(RIS="ris", SPMLE="spmle")

plotsout <- list()

## TT=4, gamma=0.25, rho=.25 - gamma estimates
plotsout[[1]] <- ggplot(subset(st, TT==4 & gamma==.25 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.25) +
#                        ggtitle(expression("" *gamma* "=0.25, " *rho* "=.25")) +
                        annotate(geom="text", x=1.5, y=-.5, label=expression("" *gamma* "=0.25, " *rho* "=.25")) +
                        theme_minimal()

## TT=4, gamma=0.25, rho=.5 - gamma estimates
plotsout[[2]] <- ggplot(subset(st, TT==4 & gamma==.25 & rho==.5)) +
                        geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name=expression("Distribution of " *gamma* " estimates")) +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.25) +
#                        ggtitle(expression("" *gamma* "=0.25, " *rho* "=.5")) +
                        annotate(geom="text", x=1.5, y=-.5, label=expression("" *gamma* "=0.25, " *rho* "=.5")) +
                        theme_minimal()

## TT=4, gamma=0.5, rho=.25 - gamma estimates
plotsout[[3]] <- ggplot(subset(st, TT==4 & gamma==.5 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        geom_hline(yintercept=.5) +
#                        ggtitle(expression("" *gamma* "=0.25, " *rho* "=.25")) +
                        annotate(geom="text", x=1.5, y=-.5, label=expression("" *gamma* "=0.5, " *rho* "=.25")) +
                        theme_minimal()

## TT=16, gamma=0.25, rho=.25 - gamma estimates
plotsout[[4]] <- ggplot(subset(st, TT==16 & gamma==.25 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.25) +
#                        ggtitle(expression("" *gamma* "=0.25, " *rho* "=.25")) +
                        annotate(geom="text", x=1.5, y=-.5, label=expression("" *gamma* "=0.25, " *rho* "=.25")) +
                        theme_minimal()

## TT=16, gamma=0.25, rho=.5 - gamma estimates
plotsout[[5]] <- ggplot(subset(st, TT==16 & gamma==.25 & rho==.5)) +
                        geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.25) +
#                        ggtitle(expression("" *gamma* "=0.25, " *rho* "=.5")) +
                        annotate(geom="text", x=1.5, y=-.5, label=expression("" *gamma* "=0.25, " *rho* "=.5")) +
                        theme_minimal()

## TT=16, gamma=0.25, rho=.5 - gamma estimates
plotsout[[6]] <- ggplot(subset(st, TT==16 & gamma==.5 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=gamma_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        geom_hline(yintercept=.5) +
#                        ggtitle(expression("" *gamma* "=0.5, " *rho* "=.25")) +
                        annotate(geom="text", x=1.5, y=-.5, label=expression("" *gamma* "=0.5, " *rho* "=.25")) +
                        theme_minimal()

pdf("plots/MC_st-ViolinPlots-gamma-200219.pdf")
  grid.arrange(arrangeGrob(plotsout[[1]], plotsout[[2]], plotsout[[3]],
                           top=textGrob("N=64, T=4", gp=gpar(fontface="bold")),
                           ncol=1),
               arrangeGrob(plotsout[[4]], plotsout[[5]], plotsout[[6]],
                           top=textGrob("N=64, T=16", gp=gpar(fontface="bold")),
                           ncol=1),
               ncol=2
               )
dev.off()





## TT=4, gamma=0.25, rho=.25 - rho estimates
plotsout[[1]] <- ggplot(subset(st, TT==4 & gamma==.25 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name=expression("" *gamma* "=0.25, " *rho* "=.25")) +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.25) +
                        theme_minimal()

## TT=4, gamma=0.25, rho=.5 - rho estimates
plotsout[[2]] <- ggplot(subset(st, TT==4 & gamma==.25 & rho==.5)) +
                        geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name=expression("" *gamma* "=0.25, " *rho* "=.5")) +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.5) +
                        theme_minimal()

## TT=4, gamma=0.5, rho=.25 - rho estimates
plotsout[[3]] <- ggplot(subset(st, TT==4 & gamma==.5 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name=expression("" *gamma* "=0.5, " *rho* "=.25")) +
                        geom_hline(yintercept=.25) +
                        theme_minimal()

## TT=16, gamma=0.25, rho=.25 - rho estimates
plotsout[[4]] <- ggplot(subset(st, TT==16 & gamma==.25 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.25) +
                        theme_minimal()

## TT=16, gamma=0.25, rho=.5 - rho estimates
plotsout[[5]] <- ggplot(subset(st, TT==16 & gamma==.25 & rho==.5)) +
                        geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        scale_x_discrete(name="") +
                        geom_hline(yintercept=.5) +
                        theme_minimal()

## TT=16, gamma=0.25, rho=.5 - rho estimates
plotsout[[6]] <- ggplot(subset(st, TT==16 & gamma==.5 & rho==.25)) +
                        geom_violin(aes(x=Estimator, y=rho_hat), fill='#A4A4A4', color="#A4A4A4", alpha = 0.5) +
                        scale_y_continuous(limits = c(-1, 1), name="") +
                        geom_hline(yintercept=.25) +
                        theme_minimal()

pdf("plots/MC_st-ViolinPlots-rho-200219.pdf")
  grid.arrange(arrangeGrob(plotsout[[1]], plotsout[[2]], plotsout[[3]],
                           top=textGrob("N=64, T=4", gp=gpar(fontface="bold")),
                           ncol=1),
               arrangeGrob(plotsout[[4]], plotsout[[5]], plotsout[[6]],
                           top=textGrob("N=64, T=16", gp=gpar(fontface="bold")),
                           ncol=1),
               ncol=2
  )
dev.off()
