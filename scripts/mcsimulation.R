library(spmle)
library(spatialprobit)
library(risprobit)
library(spatialreg)
library(McSpatial)
library(dplyr)
library(tidyr)

library(reticulate)
if(Sys.info()["sysname"]=="Linux") reticulate::use_virtualenv("r-reticulate")  # In Unix
if(Sys.info()["sysname"]=="Windows") reticulate::use_condaenv("r-reticulate")  # In Windows



# Functions ----------------------------------------------------------------------------------

fit_spprobit <- function(data, method = c('bayes', 'spmle', 'ris', 'gmm')) {
  ## Fits Spatial Probit & times result

  method <- match.arg(method)

  if (method == 'bayes') {
    perf <- system.time({
      suppressWarnings({
        bayes.fit <- sar_probit_mcmc(y = data$y, X = data$X, W = data$W_t)
      })
    })
    theta <- coef(bayes.fit)
    rho <- theta[length(theta)]
    beta <- theta[-length(theta)]
  } else if (method == 'spmle') {
    mbdm <- MBDM$new(y_list = data$y_list, X_list = data$X_list,
                     N = data$N, TT = 1,
                     W_t = data$W_t,
                     has_temporal_lag = FALSE)
    theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
    perf <- system.time(ml.fit <- mbdm$train(theta_init))
    theta <- coef(ml.fit)
    rho <- theta[length(theta)]
    beta <- theta[-length(theta)]
  } else if (method == 'ris') {
    ris <- RisSpatialProbit$new(X = data$X,
                                y = data$y,
                                W_t = data$W_t,
                                N = data$N, TT = 1,
                                temporal_dep = FALSE, spatial_dep = TRUE)
    perf <- system.time(ris$train())
    rho <- ris$rho
    beta <- ris$beta
  } else if (method == 'gmm') {
    perf <- system.time({
      tryCatch(
        expr = {
          gmm.fit <- gmmprobit(data$y ~ data$X[,2], wmat=as.matrix(data$W_t), silent=T)
          theta <- gmm.fit$coef
        },
        error = function(e){
          message("Caught an error on iteration ", i)
          print(e)
          theta <- rep(NA_real_, 3)
        }
      )

      rho <- theta[length(theta)]
      beta <- theta[-length(theta)]

    })
  }

  K <- length(beta)
  beta_names <- paste('beta', 0:(K-1), "_hat", sep = "")
  res_names <- c('time', beta_names, 'rho_hat')
  res_ls <- as.list(c(perf[3], beta, rho))
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


# Simulations ----------------------------------------------------------------------------------

## Parameter tibble
param_tb <- expand.grid(N = c(2^8),
                        rho = c(0, 0.25, 0.5),
                        beta0 = 0,
                        beta1 = 1,
                        seed = c(1:500), # Replications per config # number of MCs
                        method = c('gmm'),
                        stringsAsFactors = FALSE) %>%
  as_tibble()

## Monte Carlos
M <- nrow(param_tb)
results_ls <- vector('list', M)
for (i in 1:M) {
  params <- param_tb[i,]
  data <- simulate_spprobit(params)
  results <- fit_spprobit(data, method = params$method)
  results_ls[[i]] <- c(params, results)
}
results_tb <- bind_rows(results_ls)
print(results_tb)

## Summarise results
results_tb %>%
  group_by(method, rho) %>%
  summarise(time_mean=mean(time, na.rm=T),
            time_sd=sd(time, na.rm=T),
            beta0_hat_mean=mean(beta0_hat, na.rm=T),
            beta0_hat_sd=sd(beta0_hat, na.rm=T),
            beta1_hat_mean=mean(beta1_hat, na.rm=T),
            beta1_hat_sd=sd(beta1_hat, na.rm=T),
            rho_hat_mean=mean(rho_hat, na.rm=T),
            rho_hat_sd=sd(rho_hat, na.rm=T)
  )
