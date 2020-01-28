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

# Constants ----------------------------------------------------------------------------------

n_cores <- 11
save_results <- TRUE
save_path <- "mcresults-T-256n-500r-200128.Rdata"
n_repetitions <- 500
has_temporal_lag = TRUE
has_spatial_lag = TRUE

# Functions ----------------------------------------------------------------------------------

fit_spprobit <- function(data, method = c('bayes', 'spmle', 'ris', 'gmm', 'naiveprobit')) {
  ## Fits Spatial Probit & times result

  method <- match.arg(method)

  # tryCatch block returns error object if expr fails
  err <- tryCatch(
    expr = {
      if (method == 'spmle') {
        mbdm <- MBDM$new(y_list = data$y_list, X_list = data$X_list,
                         N = data$N, TT = data$TT,
                         W_t = data$W_t,
                         has_temporal_lag = has_temporal_lag)
        theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
        perf <- system.time({
          ml.fit <- mbdm$train(theta_init)
        })
        elapsed <- perf[3]
        theta <- coef(ml.fit)
        rho <- theta[length(theta)]
        gamma <- theta[3]
        beta <- theta[1:2]
      } else if (method == 'ris') {
        ris <- RisSpatialProbit$new(X = data$X,
                                    y = data$y,
                                    W_t = data$W_t,
                                    N = data$N, TT = data$TT,
                                    temporal_dep = has_temporal_lag,
                                    spatial_dep = has_spatial_lag)
        perf <- system.time({
          ris$train()
        })
        elapsed <- perf[3]
        gamma <- ris$gamma
        rho <- ris$rho
        beta <- ris$beta
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
    gamma <- NA
    beta <- rep(NA, K)
  }

  # Prep return value
  beta_names <- paste('beta', 0:(K-1), "_hat", sep = "")
  res_names <- c('time', beta_names, 'gamma_hat', 'rho_hat')
  res_ls <- as.list(c(elapsed, beta, gamma, rho))
  names(res_ls) <- res_names

  return(res_ls)
}

simulate_spprobit <- function(params) {
  ## Simulates some spatial probit data
  set.seed(params$seed)
  beta <- get_beta(params)
  data <- spmle::simulate_data(G = 1, N = params$N, TT = params$TT, K_outcome = 1,
                               has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                               rho = params$rho, gamma = params$gamma, beta = beta)
  return(data)
}


get_beta <- function(params) {
  as.numeric(params[grepl('beta', x = names(params))])
}


# Simulations ----------------------------------------------------------------------------------

## Parameter tibble
param_tb <- expand.grid(N = c(2^6),
                        TT = c(2^2),
                        rho = 0,
                        gamma = c(0, 0.25, 0.5),
                        beta0 = -0.5,
                        beta1 = 1,
                        seed = c(1:n_repetitions), # Replications per config (number of MCs)
                        method = c('spmle', 'ris'),
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

  ## Process / save results
if (save_results) {
  save(results_tb, file=save_path, version=3)
}

## Summarize results
results_tb %>%
  group_by(method, rho) %>%
  summarize(
    n_sim = sum(!is.na(time)),
    n_fail = sum(is.na(time)),
    time_mean=mean(time, na.rm=T),
    time_sd=sd(time, na.rm=T),
    beta0_hat_mean=mean(beta0_hat, na.rm=T),
    beta0_hat_sd=sd(beta0_hat, na.rm=T),
    beta1_hat_mean=mean(beta1_hat, na.rm=T),
    beta1_hat_sd=sd(beta1_hat, na.rm=T),
    gamma_hat_mean = mean(gamma_hat, na.rm=T),
    gamma_hat_sd = sd(gamma_hat, na.rm=T)
  )
