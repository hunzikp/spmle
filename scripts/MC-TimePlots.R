###################################################
### Create violin plots for Monte Carlo results
rm(list=ls())
library(ggplot2)
library(gridExtra)

if(Sys.info()["user"]=="ncb") setwd("/home/ncb/Dropbox/SpatialLogit/MonteCarlo/")
getwd()

###################################################
### Import Spatial MC results
files <- paste("output/spatial/", dir("output/spatial"), sep="")
sp <- data.frame()
for(f in files) {
  print(f)
  load(f)
  sp <- rbind(te, results_tb[,c("method", "N", "TT", "time")])
}
rm(files, f, results_tb)

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

### Import Temporal MC results
files <- paste("output/temporal/", dir("output/temporal"), sep="")
te <- data.frame()
for(f in files) {
  print(f)
  load(f)
  te <- rbind(te, results_tb[,c("method", "N", "TT", "time")])
}
rm(files, f, results_tb)

### Import Spatio-Temporal MC results
files <- paste("output/spatio-temporal/", dir("output/spatio-temporal"), sep="")
st <- data.frame()
for(f in files) {
  print(f)
  load(f)
  st <- rbind(st, results_tb[,c("method", "N", "TT", "time")])
}
rm(files, f, results_tb)
st <- st[st$N!=49,]


###################################################
### Plot times
## Spatial

# Prep
sp$N <- as.factor(sp$N)
sp$Estimator <- as.factor(sp$method)
levels(sp$Estimator)
levels(sp$Estimator) <- list(Bayes="bayes", GMM="gmm", NaiveProbit="naiveprobit",RIS="ris", SPMLE="spmle")

# Plot
ggplot(aes(y = time, x = Estimator, fill = N, dodge=N), data=sp) +
  geom_boxplot() +
  scale_y_log10(name="Log(Seconds)") +
#  scale_fill_manual(values=c("#F8766D", "#00BA38")) + # change colour of fill
  ggtitle("Mean Estimation Time for Spatial Monte Carlo Experiments") +
  theme_minimal()
ggsave(filename="plots/MC_sp-time-200219.pdf")


## Temporal
# Prep
#te$NT <- as.factor(paste(te$N*te$TT, " (", te$N, " x ", te$TT, ")", sep=""))
te$NT <- as.factor(paste(te$N, te$TT, sep=" x "))
te$NT <- ordered(te$NT, levels = c("64 x 4", "64 x 16"))
te$Estimator <- as.factor(te$method)
levels(te$Estimator)
levels(te$Estimator) <- list(RIS="ris", SPMLE="spmle")

# Plot
ggplot(aes(y = time, x = Estimator, fill = NT, dodge=NT), data=te) +
  geom_boxplot() +
  scale_y_log10(name="Log(Seconds)") +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) + # change colour of fill
  ggtitle("Mean Estimation Time for Temporal Monte Carlo Experiments") +
  theme_minimal()
ggsave(filename="plots/MC_te-time-200219.pdf")


## Spatio-Temporal
# Prep
#st$NT <- as.factor(passt(st$N*st$TT, " (", st$N, " x ", st$TT, ")", sep=""))
st$NT <- as.factor(paste(st$N, st$TT, sep=" x "))
st$NT <- ordered(st$NT, levels = c("64 x 4", "64 x 16"))
st$Estimator <- as.factor(st$method)
levels(st$Estimator)
levels(st$Estimator) <- list(RIS="ris", SPMLE="spmle")

# Plot
ggplot(aes(y = time, x = Estimator, fill = NT, dodge=NT), data=st) +
  geom_boxplot() +
  scale_y_log10(name="Log(Seconds)") +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) + # change colour of fill
  ggtitle("Mean Estimation Time for Spatio-Temporal Monte Carlo Experiments") +
  theme_minimal()
ggsave(filename="plots/MC_st-time-200219.pdf")








