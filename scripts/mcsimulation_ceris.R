n_cores <- 11
save_results <- TRUE
save_path <- "mcresults-S-256n-500r-CERIS-200207.Rdata"
n_repetitions <- 500
has_temporal_lag = FALSE
has_spatial_lag = TRUE


ris_estimator <- function(y, X, W, control = list(...), debug = FALSE)
{
  con <- list(
    R = 1000,
    store.density = FALSE,
    optim.max.iterations = 1000,
    optim.reltol = .0025
  )
  con[names(control)] <- control
  with(con, {
    n <- length(y)
    k <- dim(X)[2]
    Z <- diag(as.vector(1 - 2 * y))
    In <- diag(n)
    ## Random draw from importance density function, given
    ## an upper bound, using antithetical sampling
    random.draw <- function(upper.bound)
    {
      q <- runif(R/2)
      q <- c(q, 1-q)
      qnorm(q * pnorm(upper.bound))
    }
    iter <- 0
    densities <- NULL
    ll <- function(par)
    {
      iter <<- iter + 1
      rho <-
        -1 + 2 * pnorm(par[k+1])
      beta <- par[1:k]
      Ai <- -rho * W
      diag(Ai) <- 1
      A<-solve(Ai)
      omega <- Z %*% A %*% t(A) %*% t(Z)
      V <- -Z %*% A %*% X %*% beta
      B <- solve(chol(solve(omega)))
      Eta0 <- Eta <- matrix(NA, nrow = n, ncol = R)
      Eta0[n,] <- rep(1/B[n,n] * V[n], R)
      Eta[n,] <- random.draw(Eta0[n,])
      for (i in (n-1):1)
      {
        Eta0[i,] <- 1/B[i,i] * (V[i] - t(B[i,(i+1):n]) %*% Eta[(i+1):n,])
        if (i > 1)
          Eta[i,] <- random.draw(Eta0[i,])
      }
      pnEta0 <- pnorm(Eta0)
      a <- 1 / mean(pnEta0)
      ll <- log(mean(apply(pnEta0 * a, 2, prod))) - log(1/R) - n * log(a)
      ll
    }
    par <- optim(c(rep(0,k), 0), ll, control=list(fnscale = -1,
                                                  reltol = optim.reltol, maxit = optim.max.iterations))
    list(beta = par$par[1:k],
         rho = -1 + 2 * pnorm(par$par[k+1]),
         beta.se = rep(NA, k),
         rho.se = NA,
         densities = densities)
  })
}

fit_spprobit <- function(data, method = c('bayes', 'spmle', 'ris', 'gmm', 'naiveprobit', 'ceris')) {
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
      } else if (method == 'ceris') {
        perf <- system.time({
          ceris.fit <- ris_estimator(y=data$y,
                                     X=data$X,
                                     W=as.matrix(data$W_t),
                                     control=list(fnscale = -1,
                                                  reltol = 0.025,
                                                  maxit = 1000),
                                     debug=FALSE
          )
        })
        elapsed <- perf[3]
        rho <- ceris.fit$rho
        beta <- ceris.fit$beta
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


# Simulations ----------------------------------------------------------------------------------

## Parameter tibble
param_tb <- expand.grid(N = c(2^8),
                        rho = c(0.25, 0.5),
                        beta0 = -0.5,
                        beta1 = 1,
                        seed = c(1:n_repetitions), # Replications per config (number of MCs)
                        method = c('ceris'),
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
#    n_fail = sum(is.na(time)),
    time_mean=mean(time, na.rm=T),
    time_sd=sd(time, na.rm=T),
    beta0_hat_mean=mean(beta0_hat, na.rm=T),
    beta0_hat_sd=sd(beta0_hat, na.rm=T),
    beta1_hat_mean=mean(beta1_hat, na.rm=T),
    beta1_hat_sd=sd(beta1_hat, na.rm=T),
    rho_hat_mean=mean(rho_hat, na.rm=T),
    rho_hat_sd=sd(rho_hat, na.rm=T)
  )

