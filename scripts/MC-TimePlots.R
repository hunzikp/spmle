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
files
files <- files[c(1:3,5:6)]
sp <- data.frame()
for(f in files) {
  print(f)
  load(f)
  sp <- rbind(sp, results_tb[!is.na(results_tb$rho_se),c("method", "N", "time", "rho_se")])
}
rm(files, f, results_tb)
unique(paste(sp$method, sp$N, sep="-"))

### Import Temporal MC results
files <- paste("output/temporal/", dir("output/temporal"), sep="")
te <- data.frame()
for(f in files) {
  print(f)
  load(f)
  te <- rbind(te, results_tb[!is.na(results_tb$rho_se),c("method", "N", "TT", "time")])
}
rm(files, f, results_tb)
unique(paste0(te$method, st$N, st$TT))

### Import Spatio-Temporal MC results
files <- paste("output/spatio-temporal/", dir("output/spatio-temporal"), sep="")
st <- data.frame()
for(f in files) {
  if(f!="mcresults-ST-49n20t-RIS-100r-200212.Rdata"){
    print(f)
    load(f)
    st <- rbind(st, results_tb[!is.na(results_tb$rho_se), c("method", "N", "TT", "time")])
  }
}
rm(files, f, results_tb)
unique(paste0(st$method, st$N, st$TT))

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
te$NT <- ordered(te$NT, levels = c("64 x 4", "64 x 16", "256 x 16"))
te$Estimator <- as.factor(te$method)
levels(te$Estimator)
levels(te$Estimator) <- list(RIS="ris", SPMLE="spmle")

# Plot
ggplot(aes(y = time, x = Estimator, fill = NT, dodge=NT), data=te) +
  geom_boxplot() +
  scale_y_log10(name="Log(Seconds)") +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#4E84C4")) + # change colour of fill
  ggtitle("Mean Estimation Time for Temporal Monte Carlo Experiments") +
  theme_minimal()
ggsave(filename="plots/MC_te-time-200228.pdf")


## Spatio-Temporal
# Prep
#st$NT <- as.factor(passt(st$N*st$TT, " (", st$N, " x ", st$TT, ")", sep=""))
st$NT <- as.factor(paste(st$N, st$TT, sep=" x "))
st$NT <- ordered(st$NT, levels = c("64 x 4", "64 x 16", "256 x 16"))
st$Estimator <- as.factor(st$method)
levels(st$Estimator)
levels(st$Estimator) <- list(RIS="ris", SPMLE="spmle")

# Plot
ggplot(aes(y = time, x = Estimator, fill = NT, dodge=NT), data=st) +
  geom_boxplot() +
  scale_y_log10(name="Log(Seconds)") +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#4E84C4")) + # change colour of fill
  ggtitle("Mean Estimation Time for Spatio-Temporal Monte Carlo Experiments") +
  theme_minimal()
ggsave(filename="plots/MC_st-time-200228.pdf")








