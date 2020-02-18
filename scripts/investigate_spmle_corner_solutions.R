###################################################
### Investigate 'corner-solution' results for SPMLE

###################################################
### Header
rm(list=ls())
library(spmle)
library(spatialprobit)
library(risprobit)
library(spatialreg)
library(McSpatial)
library(car)  # Silent dependency of McSpatial

library(dplyr)
library(tidyr)
library(doParallel)

library(reticulate)
if(Sys.info()["sysname"]=="Linux") reticulate::use_virtualenv("r-reticulate")  # In Unix
if(Sys.info()["sysname"]=="Windows") reticulate::use_condaenv("r-reticulate")  # In Windows

if(Sys.info()["user"]=="ncb") setwd("/home/ncb/Dropbox/SpatialLogit/MonteCarlo/")
getwd()

###################################################
### Functions

# Constants ----------------------------------------------------------------------------------

n_cores <- 10
save_results <- FALSE

# Functions ----------------------------------------------------------------------------------

fit_spprobit <- function(data, method = c('bayes', 'spmle', 'ris', 'gmm', 'naiveprobit')) {
  ## Fits Spatial Probit & times result

  method <- match.arg(method)

  # tryCatch block returns error object if expr fails
  err <- tryCatch(
    expr = {
      if (method == 'bayes') {
        perf <- system.time({
          suppressWarnings({
            bayes.fit <- sar_probit_mcmc(y = data$y, X = data$X, W = data$W_t)
          })
        })
        elapsed <- perf[3]
        theta <- coef(bayes.fit)
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
      } else if (method == 'spmle') {
        mbdm <- MBDM$new(y_list = data$y_list, X_list = data$X_list,
                         N = data$N, TT = 1,
                         W_t = data$W_t,
                         has_temporal_lag = FALSE)
        theta_init <- runif(mbdm$n_parameters, min = -.75, max = 0.75)
        perf <- system.time({
          ml.fit <- mbdm$train(theta_init)
        })
        elapsed <- perf[3]
        theta <- coef(ml.fit)
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
      } else if (method == 'ris') {
        ris <- RisSpatialProbit$new(X = data$X,
                                    y = data$y,
                                    W_t = data$W_t,
                                    N = data$N, TT = 1,
                                    temporal_dep = FALSE, spatial_dep = TRUE)
        perf <- system.time({
          ris$train()
        })
        elapsed <- perf[3]
        rho <- ris$rho
        beta <- ris$beta
      } else if (method == 'gmm') {
        perf <- system.time({
          gmm.fit <- gmmprobit(data$y ~ data$X[,2], wmat=as.matrix(data$W_t), silent=T)
        })
        elapsed <- perf[3]
        theta <- gmm.fit$coef
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
      } else if (method == 'naiveprobit') {
        Wy <- as.matrix(data$W_t) %*% data$y
        perf <- system.time({
          glm.fit <- glm(data$y~ data$X[,2] + Wy, family=binomial(link="probit"))
        })
        elapsed <- perf[3]
        theta <- coef(glm.fit)
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
      }
    },
    error = function(e) {
      return(e)
    }
  )

  K <- ncol(data$X)

  # Deal with error
  if (inherits(err, 'error')) {
    elapsed <- NA
    rho <- NA
    beta <- rep(NA, K)
  }

  # Prep return value
  beta_names <- paste('beta', 0:(K-1), "_hat", sep = "")
  res_names <- c('time', beta_names, 'rho_hat')
  res_ls <- as.list(c(elapsed, beta, rho))
  names(res_ls) <- res_names

  return(res_ls)
}

simulate_spprobit <- function(params) {
  ## Simulates some spatial probit data
  set.seed(params$seed)
  beta <- get_beta(params)
  data <- spmle::simulate_data(G = 1, N = params$N, TT = 1, K_outcome = 1,
                               has_temporal_lag = FALSE, has_spatial_lag = TRUE,
                               rho = params$rho, beta = beta)
  return(data)
}

get_beta <- function(params) {
  as.numeric(params[grepl('beta', x = names(params))])
}


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
### Identify corner solutions
## n=256
sp1[sp1$method=="spmle" & sp1$rho_hat<=-0.75, c("seed", "rho", "rho_hat")]
# export seeds
seeds1 <- sp1[sp1$method=="spmle" & sp1$rho_hat<=-0.75, c("seed")]

## n=4096
sp3[sp3$method=="spmle" & sp3$rho_hat<=-0.75, c("seed", "rho", "rho_hat")]
sp3[sp3$method=="spmle" & sp3$rho_hat>=0.9, c("seed", "rho", "rho_hat")]
# export seeds
seeds3 <- unlist(sp3[sp3$method=="spmle" & sp3$rho_hat<=-0.75, c("seed")])
length(seeds3)
## n=16384
sp4[sp4$method=="spmle" & sp4$rho_hat<=-0.75, c("seed", "rho", "rho_hat")]
sp4[sp4$method=="spmle" & sp4$rho_hat>=0.9, c("seed", "rho", "rho_hat")]
# export seeds
seeds4 <- unlist(sp4[sp4$method=="spmle" & sp4$rho_hat<=-0.75, c("seed")])
length(seeds4)

###################################################
### Rerun MCs
## Parameter tibble
param_tb <- expand.grid(N = c(2^8),
                        rho = c(0, 0.5),
                        beta0 = -0.5,
                        beta1 = 1,
                        seed = c(seeds1), # Replications per config (number of MCs)
                        method = c('spmle'),
                        stringsAsFactors = FALSE) %>%
  as_tibble()


## Prep parallelization
cl <- makeCluster(n_cores)
registerDoParallel(cl)

## Monte Carlos
M <- nrow(param_tb)
results_tb <- foreach (i = 1:M, .packages = c("spmle",
                                              "spatialprobit",
                                              "risprobit",
                                              "spatialreg",
                                              "McSpatial"),
                       .combine = bind_rows) %dopar% {

                         params <- param_tb[i,]
                         data <- simulate_spprobit(params)
                         results <- fit_spprobit(data, method = params$method)

                         c(params, results)
                       }
stopCluster(cl)

print(results_tb)
dim(results_tb[results_tb$rho_hat<=-0.75, c("seed", "rho", "rho_hat")])
dim(results_tb[results_tb$rho_hat>=.9, c("seed", "rho", "rho_hat")])


## Results for setting starting values from [-.75,.75] to [0,0.5]
# n=256 # 7 instead of 11 trials with rho <= -0.75
# n=1024 # ?
# n=4096 # 7 instead of 36 trials with rho <= -0.75 (includes estimates of rho==.75)
# n=19384 #  instead of 17 trials with rho <= -0.75 (includes estimates of rho==.75)


###################################################
### Investigate individual runs
####################################
## n=256
set.seed(104)

## Set data parameters
G <- 1
N <- 2^8
TT <- 1
K_outcome <- 1
has_temporal_lag = FALSE
has_spatial_lag = TRUE

## Make the data
data.ls <- spmle::simulate_data(G, N, TT, K_outcome,
                                has_temporal_lag = has_temporal_lag,
                                has_spatial_lag = has_spatial_lag,
                                rho = 0, beta = rep(c(-.5,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
summary(ml.fit)
ml.fit$iterations
ml.fit$message


## seed=70, rho=0.5
# 129 iterations, successful convergence
# Warning messages:
# 1: In sqrt(diag(vc)) : NaNs produced
# 2: In sqrt(diag(vc)) : NaNs produced
# no SEs for beta0 and rho

## seed=100, rho=0.5
# 114 iterations, successful convergence
# same warning messages and SEs as above

## seed=104, rho=0
# 80 iterations, successful convergence
# same warning messsages and SEs as above

####################################
## n=1024
set.seed(83)

## Set data parameters
G <- 1
N <- 2^12
TT <- 1
K_outcome <- 1
has_temporal_lag = FALSE
has_spatial_lag = TRUE

## Make the data
data.ls <- spmle::simulate_data(G, N, TT, K_outcome,
                                has_temporal_lag = has_temporal_lag,
                                has_spatial_lag = has_spatial_lag,
                                rho = 0.25, beta = rep(c(-.5,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
summary(ml.fit)
ml.fit$iterations
ml.fit$message


## seed=63, rho=0.5
# 110 iterations, successful convergence
# Warning messages:
  # 1: In sqrt(diag(vc)) : NaNs produced
  # 2: In sqrt(diag(vc)) : NaNs produced
  # no SEs for beta0 and rho

## seed=67, rho=0.75
# 179 iterations, successful convergence
# same warning messages and SEs as above

## seed=83, rho=0.25
# 130 iterations, successful convergence
# same warning messages and SEs as above


