##############################################################
# Main class definition
##############################################################

#### Fields & Constructor
STModel <- R6Class("STModel",
                 public = list(

                   #### CLASS FIELDS
                   ## Data param fields
                   N = NULL,
                   TT = NULL,
                   NT = NULL,
                   K = NULL,

                   ## Data fields
                   y = NULL,
                   X = NULL,

                   ## Weights matrices
                   TM = NULL,  # temp lag matrix
                   W_t = NULL,
                   W = NULL,

                   ## Parameter indices
                   beta_indices = NULL,
                   gamma_index = NULL,
                   rho_index = NULL,
                   n_parameters = NULL,
                   n_auto_params = NULL,

                   ## Parameter names
                   beta_names = NULL,
                   auto_param_names = NULL,

                   ## Switches
                   has_temporal_lag = NULL,
                   has_spatial_lag = NULL,

                   ## Neumann-related objects
                   L = NULL,
                   neumann_diag_list = NULL,

                   ## Zero period objects
                   zero_period = NULL,
                   X_expectation = NULL,

                   ## Whether to use python for linalg speed-up
                   use_py = NULL,

                   #### CONSTRUCTOR
                   initialize = function(y, X, N, TT,
                                         W_t = NULL,
                                         has_temporal_lag = TRUE,
                                         has_spatial_lag = TRUE,
                                         zero_period = TRUE,
                                         use_py = FALSE) {

                     ## Some argument checking
                     if (length(y) != nrow(X)) {
                       stop('y and X dimensions do not match.')
                     }
                     if (length(y) != N*TT) {
                       stop('y has wrong length (check N and TT).')
                     }
                     if (TT == 1 & has_temporal_lag) {
                       stop('Estimating a temporal lag requires more than one period.')
                     }
                     if (has_spatial_lag) {
                       if (is.null(W_t)) {
                         stop('W_t missing.')
                       }
                       if (nrow(W_t) != N | ncol(W_t) != N) {
                         stop('W_t has wrong dimensions.')
                       }
                     }

                     ## Set py preference
                     self$use_py = use_py

                     ## Set data parameters
                     self$N <- N
                     self$TT <- TT
                     self$NT <- N*TT

                     ## Set y, X
                     self$y <- y
                     self$X <- X
                     self$K <- ncol(self$X)

                     ## Get (or make) beta names
                     self$K <- ncol(self$X)
                     beta_names <- colnames(X)
                     if (is.null(beta_names)) {
                       beta_names <- paste0('beta_', 1:self$K)
                     }
                     self$beta_names <- beta_names

                     ## Set spatial weigths
                     self$has_spatial_lag <- has_spatial_lag
                     self$W_t <- W_t
                     if (has_spatial_lag & is.null(W_t)) {
                       stop('Missing spatial weights matrix.')
                     }

                     # Create the temporal lag weights matrix
                     self$has_temporal_lag <- has_temporal_lag
                     if (has_temporal_lag) {
                       self$TM <- self$make_temporal_lag_matrix()
                     }

                     # Create the block-diagonal spatial weights matrix
                     if (self$has_spatial_lag) {
                       self$W <- bdiag(rep(list(self$W_t), TT))
                     }

                     ## Precompute diagonals for Neumann approximation to d
                     if (self$has_spatial_lag) {
                       self$L <-  8L  # We use the 8th order Neumann series to approximate
                       neumann_list <- vector('list', self$L)
                       neumann_diag_list <- vector('list', self$L)
                       R <- self$W_t
                       neumann_list[[1]] <- R
                       neumann_diag_list[[1]] <- diag(R)
                       for (i in 2:self$L) {
                         neumann_list[[i]] <- neumann_list[[i-1]]%*%R
                         neumann_diag_list[[i]] <- diag(neumann_list[[i]])
                       }
                       self$neumann_diag_list <- neumann_diag_list
                       rm(neumann_list)
                     }

                     ## Determine location of beta/rho/gamma in theta parameter vector
                     n_parameters <- self$K
                     auto_param_names <- c()
                     self$beta_indices <- 1:self$K
                     idx <- self$K + 1
                     if (self$has_spatial_lag) {
                       self$rho_index <- idx
                       idx <- idx + 1
                       n_parameters <- n_parameters + 1
                       auto_param_names <- c(auto_param_names, 'rho')
                     }
                     if (self$has_temporal_lag) {
                       self$gamma_index <- idx
                       n_parameters <- n_parameters + 1
                       auto_param_names <- c(auto_param_names, 'gamma')
                     }
                     self$n_parameters <- n_parameters
                     self$n_auto_params <- length(auto_param_names)
                     self$auto_param_names <- auto_param_names

                     ## Compute X_expectation
                     self$zero_period <- zero_period
                     if (zero_period) {
                       X_t_list <- vector('list', TT)
                       for (t in 1:TT) {
                         start <- (t-1)*N + 1
                         end <- start + N - 1
                         X_t_list[[t]] <- X[start:end,]
                       }
                       self$X_expectation <- Reduce('+', X_t_list)/TT
                     }
                   }

                 ))


#### Returns temporal lag matrix
STModel$set("public", "make_temporal_lag_matrix",
         function() {

           index.mat <- cbind(1:self$NT, 1:self$NT)
           index.mat[,1] <- index.mat[,1]+self$N
           index.mat <- index.mat[index.mat[,1] <= self$NT,]
           TM <- sparseMatrix(i = index.mat[,1], j = index.mat[,2], dims = c(self$NT, self$NT))

           return(TM)
         })


#### Methods extracting specific parameter vectors from theta
STModel$set("public", "extract_beta",
         function(theta) {
           beta <- theta[self$beta_indices]
           return(beta)
           })

STModel$set("public", "extract_rho",
         function(theta) {
           rho <- theta[self$rho_index]
           names(rho) <- NULL
           return(rho)
         })

STModel$set("public", "extract_gamma",
         function(theta) {
           gamma <- theta[self$gamma_index]
           names(gamma) <- NULL
           return(gamma)
         })


#### Method for returning E(ystar | W) - for computation of initial ystar
STModel$set("public", "get_ystar_expectation",
         function(theta) {

           ## E(Xbeta)
           beta <- self$extract_beta(theta)
           Xbeta <- self$X_expectation%*%beta

           # Multiplier
           M <- .sparseDiagonal(n = self$N)
           if (self$has_temporal_lag) {
             gamma <- self$extract_gamma(theta)
             M <- M - .sparseDiagonal(n = self$N)*gamma
           }
           if (self$has_spatial_lag) {
             rho <- self$extract_rho(theta)
             M <- M - rho*self$W_t
           }

           ystar_expectation <- as.vector(solve(M, Xbeta))

           return(ystar_expectation)
         })



#### Method for computing mu if there's a temporal lag
STModel$set("public", "compute_mu_full",
         function(theta) {

          stopifnot(self$has_temporal_lag)

          ## Get gamma
         gamma <- self$extract_gamma(theta)

          ## Xbeta
          beta <- self$extract_beta(theta)
          Xbeta <- self$X%*%beta

          ## Initial (zero) period
          if (self$zero_period) {
            ystar_0 <- self$get_ystar_expectation(theta)
            # For t=1, ystar_1 = X_1%*%beta + rho*W%*%ystar_1 + gamma*ystar_0 + u
            Xbeta[1:self$N,] <- Xbeta[1:self$N,] + gamma*ystar_0
          }

          ## A (the full weights matrix)
          A <- .sparseDiagonal(self$NT)
          if (self$has_spatial_lag) {
            rho <- self$extract_rho(theta)
            A <- A - self$W*rho
          }
          A <- A - self$TM*gamma

          ## Get mu
          if (self$use_py) {
            mu <- scipy$sparse$linalg$bicgstab(r_to_py(A), Xbeta)[[1]]
          } else {
            mu <- solve(A, Xbeta)
          }
          mu <- as.vector(mu)

          return(mu)
         })


#### Method for computing mu if there's NO temporal lag, but a spatial lag
STModel$set("public", "compute_mu_blockwise",
          function(theta) {

            stopifnot(!self$has_temporal_lag)

            ## Xbeta
            beta <- self$extract_beta(theta)
            Xbeta <- self$X%*%beta

            ## A_t (the period-wise weights matrix)
            rho <- self$extract_rho(theta)
            A_t <- .sparseDiagonal(self$N) - self$W_t*rho

            ## Get each mu_t, concatenate
            mu <- c()
            for (t in 1:self$TT) {
              start <- (t-1)*self$N + 1
              end <- start + self$N - 1
              Xbeta_t <- Xbeta[start:end,,drop=FALSE]
              if (self$use_py) {
                mu_t <- scipy$sparse$linalg$bicgstab(r_to_py(A_t), Xbeta_t)[[1]]
              } else {
                mu_t <- solve(A_t, Xbeta_t)
              }
              mu <- c(mu, as.vector(mu_t))
            }

            return(mu)
          })


#### Method for computing mu if there are no dependencies
STModel$set("public", "compute_mu_independent",
         function(theta) {

           beta <- self$extract_beta(theta)
           Xbeta <- self$X%*%beta
           mu <- as.vector(Xbeta)

           return(mu)
         })


#### Compute d using Neumann approximation
STModel$set("public", "compute_d_neumann",
          function(rho) {

            ## Multiply each diagonal with power of rho scalar
            rho_vec <- rho^(1:self$L)
            rw_list <- lapply(1:self$L, function(i) rho_vec[i] * self$neumann_diag_list[[i]])
            d_t <- rep(1, self$N) + Reduce('+', rw_list)

            ## Full diagonal
            d <- rep(d_t, self$TT)

            return(d)
          })


#### Compute likelihood
STModel$set("public", "compute_llik",
          function(theta, verbose = FALSE, max_weights = 1) {


            ## Print the candidat
            if (verbose) {
              print(theta)
              flush.console()
            }

            ## Return NA if weights exceed max_weights
            auto_params <- c()
            if (self$has_spatial_lag) {
              auto_params <- c(auto_params, self$extract_rho(theta))
            }
            if (self$has_temporal_lag) {
              auto_params <- c(auto_params, self$extract_gamma(theta))
            }

            if (length(auto_params) > 0) {
              if (any(abs(auto_params) > max_weights)) {
                return(NA)
              }
              if (max(abs(auto_params)) > 1) {
                return(NA)
              }
            }

            ## Compute mu
            if (self$has_temporal_lag) {
              mu <- self$compute_mu_full(theta)
            } else if (self$has_spatial_lag) {  # only spatial lag
              mu <- self$compute_mu_blockwise(theta)
            } else {
              mu <- self$compute_mu_independent(theta)
            }


            ## Compute d
            if (self$has_spatial_lag) {
              rho <- self$extract_rho(theta)
              d <- self$compute_d_neumann(rho)
            } else {
              d <- rep(1, self$NT)
            }

            ## Compute likelihood
            v <- mu/d
            log_pr <- pnorm(v, log.p = TRUE, lower.tail = TRUE)
            log_not_pr <- pnorm(v, log.p = TRUE, lower.tail = FALSE)
            llik <- self$y*log_pr + (1-self$y)*log_not_pr

            ## Print the llik if verbose
            if (verbose) {
              print(sum(llik))
              flush.console()
            }

            return(sum(llik))
          })


#### Train model
STModel$set("public", "train",
         function(theta_init, verbose = FALSE, max_weights = 1, method = "BFGS") {

           theta_names <- c(self$beta_names, self$auto_param_names)
           names(theta_init) <- theta_names

           ml.fit <- maxLik(logLik = self$compute_llik, start = theta_init, method = method,
                            verbose = verbose, max_weights = max_weights)
           return(ml.fit)

         })


#### Train model with automatic initial parameter selection
STModel$set("public", "autotrain",
         function(method = "BFGS") {

           glm_fit <- glm.fit(x = self$X, y = self$y, family = binomial(link = 'probit'))
           glm_coef <- glm_fit$coefficients
           theta_init <- c(glm_coef, rep(0.1, self$n_auto_params))

           ml.fit <- self$train(theta_init = theta_init, method = method)
           return(ml.fit)
         })
