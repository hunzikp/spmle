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

n_cores <- 10
save_results <- FALSE
save_path <- "mcresults-256n-500r-200123.Rdata"
n_repetitions <- 2


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
        theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
        perf <- system.time({
          ml.fit <- mbdm$train(theta_init)
        })
        elapsed <- perf[3]
        theta <- coef(ml.fit)
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
        ses <- sqrt(diag(solve(-ml.fit$hessian)))
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
        ses <- sqrt(diag(ris$VC_mat))
      } else if (method == 'gmm') {
        perf <- system.time({
          gmm.fit <- gmmprobit(data$y ~ data$X[,2], wmat=as.matrix(data$W_t), silent=T)
        })
        elapsed <- perf[3]
        theta <- gmm.fit$coef
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
        ses <- gmm.fit$se
      } else if (method == 'naiveprobit') {
        Wy <- as.matrix(data$W_t) %*% data$y
        perf <- system.time({
          glm.fit <- glm(data$y~ data$X[,2] + Wy, family=binomial(link="probit"))
        })
        elapsed <- perf[3]
        theta <- coef(glm.fit)
        rho <- theta[length(theta)]
        beta <- theta[-length(theta)]
        ses <- sqrt(diag(vcov(glm.fit)))
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
rm(beta_names)
  # Prep return value
  beta_names <- paste('beta', 0:(K-1), rep(c("_hat", "_sd"), each=2), sep = "")
  res_names <- c('time', beta_names, 'rho_hat', 'rho_sd')
  res_ls <- as.list(c(elapsed, beta, ses[1:2], rho, ses[3]))
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
                        beta0 = -0.5,
                        beta1 = 1,
                        seed = c(1:n_repetitions), # Replications per config (number of MCs)
                        method = c('naiveprobit', 'gmm'),
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
  save(results_tb, file=save_path)
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
    rho_hat_mean=mean(rho_hat, na.rm=T),
    rho_hat_sd=sd(rho_hat, na.rm=T)
  )
