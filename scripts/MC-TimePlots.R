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
### Plot times

sp_t <- ggplot(sp) +
          geom_boxplot(aes(x=as.factor(method), y=time)) +
          theme_bw()

ggplot(aes(y = time, x = as.factor(method), fill = as.factor(N), dodge=as.factor(N)), data=sp) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()

