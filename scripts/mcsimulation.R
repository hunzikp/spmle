library(spmle)
library(spatialprobit)
library(risprobit)
library(spatialreg)
library(dplyr)
library(tidyr)


# Functions ----------------------------------------------------------------------------------

fit_spprobit <- function(data, method = c('bayes', 'spmle', 'ris')) {
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
param_tb <- expand.grid(N = c(10^2),
                        rho = c(-0.5, 0.5),
                        beta0 = 0,
                        beta1 = 1,
                        seed = c(1:2), # Replications per config
                        method = c('bayes', 'spmle', 'ris'),
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
