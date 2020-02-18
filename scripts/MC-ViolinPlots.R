###################################################
### Create violin plots for Monte Carlo results
rm(list=ls())
library(ggplot2)
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

###################################################
### Plot rho estimates
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

### Plot rho estimates
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

### Plot
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

### Plot
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
### Plot times

sp1_t <- ggplot(sp1) +
          geom_boxplot(aes(x=as.factor(method), y=time)) +
          theme_bw()

sp2_t <- ggplot(sp2) +
          geom_boxplot(aes(x=as.factor(method), y=time)) +
          theme_bw()

sp3_t <- ggplot(sp3) +
          geom_boxplot(aes(x=as.factor(method), y=time)) +
          theme_bw()

sp4_t <- ggplot(sp4) +
          geom_boxplot(aes(x=as.factor(method), y=time)) +
          theme_bw()






