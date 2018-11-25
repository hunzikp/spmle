
##############################################################
# Installation function for python dependencies
##############################################################

install_pydep <- function(method = "auto", conda = "auto") {
  reticulate::py_install("numpy", method = method, conda = conda)
  reticulate::py_install("scipy", method = method, conda = conda)
}


##############################################################
# Helper functions
##############################################################

## Conversion from dgCMatrix -> csc_matrix
as_csc <- function(A) {

  if (inherits(A, what = "dgCMatrix")) {

    scipy$sparse$csc_matrix(tuple(A@x, A@i, A@p),
                            dtype = np$float64,
                            shape = tuple(dim(A)[1], dim(A)[2]))

  } else if (inherits(A, what = "ngCMatrix")) {

    x <- rep(1L, length(A@i))
    scipy$sparse$csc_matrix(tuple(x, A@i, A@p),
                            dtype = np$int32,
                            shape = tuple(dim(A)[1], dim(A)[2]))
  }
}


## Make covariate matrix
make_X <- function(X_list, N, TT) {
  ## Creates a covariate matrix for multiple outcomes
  # X_list: list containing one covariate matrix per outcome
  # N: Number of units
  # TT: Number of periods

  G <- length(X_list)  # Number of outcomes
  K_vec <- unlist(lapply(X_list, ncol))  # Number of regressors per outcome
  K <- sum(K_vec)  # Total number of regressors
  X <- matrix(0, N*TT*G, K)
  for (t in 1:TT) {
    # Determine which rows of X_list to pick
    t_start <- (t-1)*N + 1
    t_end <- t_start + N - 1
    for (k in 1:G) {
      # Determine which columns of X to fill in
      if (k == 1) {
        col_start <- 1
        col_end <- K_vec[1]
      } else {
        col_start <- sum(K_vec[1:(k-1)]) + 1
        col_end <- col_start + K_vec[k] - 1
      }
      # Determine which rows of X to fill in
      row_start <- (t-1)*N*G + (k-1)*N + 1
      row_end <- row_start + N - 1
      # Fill in
      X[row_start:row_end, col_start:col_end] <- X_list[[k]][t_start:t_end,]
    }
  }

  return(X)
}

## Make spatial weights matrix from lattice
# N must have natural root
make_W_t <- function(N) {

  ras <- raster(matrix(1, sqrt(N), sqrt(N)))
  spdf <- rasterToPolygons(ras)
  B <- gIntersects(spdf, byid = TRUE)
  diag(B) <- FALSE
  B <- B*1
  W_t <- B / rowSums(B)
  W_t <- Matrix(W_t, sparse = TRUE)

  return(W_t)
}

## Make period-wise spatial weights matrices for multiple outcomes
make_spatial_weights <- function(W_t, G) {

  matrix_list <- vector('list', G)
  N <- nrow(W_t)
  sparseZeroMatrix <- Matrix(matrix(0, N, N), sparse = TRUE)

  # Make a block-digonal matrix with one block per outcome
  for (k in 1:G) {
    bmat_list <- rep(list(sparseZeroMatrix), G)
    bmat_list[[k]] <- W_t
    matrix_list[[k]] <- bdiag(bmat_list)
  }

  return(matrix_list)
}

## Make period-wise outcome-interdependence weights matrices
make_outcome_weights <- function(N, G) {

  matrix_list <- vector('list', G)

  for (k in 1:G) {
    not_k_vec <- (1:G)[-k]
    ivec <- (k-1)*N + rep(1:N, G-1)
    jvec <- c()
    for (not_k in not_k_vec) {
      this_j <- (not_k-1)*N + 1:N
      jvec <- c(jvec, this_j)
    }
    this_M <- sparseMatrix(i = ivec, j = jvec, dims = c(N*G, N*G))
    matrix_list[[k]] <- this_M
  }

  return(matrix_list)
}

## Make temporal-lag weights matrix
make_temporal_lag_matrix <- function(N, G, TT) {

  NGT <- N*G*TT
  index.mat <- cbind(1:NGT, 1:NGT)
  index.mat[,1] <- index.mat[,1]+N*G
  index.mat <- index.mat[index.mat[,1] <= NGT,]
  TM <- sparseMatrix(i = index.mat[,1], j = index.mat[,2], dims = c(NGT, NGT))

  return(TM)
}




##############################################################
# Takahashi inverse
##############################################################

takahashi_inverse <- function(A) {
  # A: Sparse square matrix
  # See https://www.mathworks.com/matlabcentral/fileexchange/33966-sparseinv--sparse-inverse-subset

  # Perform LDU decomposition, such that C = P*A*Q = (L+I)*D*(U+I)
  LU.ls <- expand(lu(A))
  P <- LU.ls$P
  U <- LU.ls$U
  L <- LU.ls$L
  Q <- t(LU.ls$Q)

  d <- diag(U)
  D <- diag(d)
  U <- U/d

  I <- .sparseDiagonal(nrow(A))
  U <- U-I
  L <- L-I

  # C <- P%*%A%*%Q

  # Get Zpattern
  # Zpattern is the symbolic Cholesky factorization of C+C',
  # so it includes all entries in L+U and its transpose.
  Zpattern <- (t(L + U + I)!=0)

  # Apply Takahashi equations using sparseinv implementation
  Z <- sparseinv:::.sparseinv_wrapper(L = L, d = d, U = t(U), Zpattern = Zpattern)

  # Permute Z so that it has the same structure as A
  Z <- Q%*%Z%*%P

  return(Z)
}


##############################################################
# Main class definition
##############################################################

#### Fields & Constructor
MBDM <- R6Class("MBDM",
                 public = list(

                   #### CLASS FIELDS
                   ## Data param fields
                   G = NULL,
                   N = NULL,
                   TT = NULL,
                   NGT = NULL,
                   NG = NULL,
                   K = NULL,

                   ## Data fields
                   y = NULL,
                   X = NULL,

                   ## Weights matrices
                   omega_list = NULL,
                   omega_t_list = NULL,

                   ## Weight parameters
                   eta_names = NULL,
                   eta_t_names = NULL,
                   n_eta = NULL,
                   n_eta_t = NULL,

                   ## Parameter indices
                   beta_indices = NULL,
                   eta_indices = NULL,
                   eta_t_indices = NULL,
                   n_parameters = NULL,
                   beta_names = NULL,

                   ## Switches
                   use_neumann = NULL,
                   has_temporal_lag = NULL,

                   ## Neumann-related objects
                   L = NULL,
                   neumann_diag_list = NULL,

                   ## Python objects
                   I_py = NULL,
                   omega_list_py = NULL,
                   omega_t_list_py = NULL,
                   I_t_py = NULL,


                   #### CONSTRUCTOR
                   initialize = function(y_list, X_list, N, TT,
                                         W_t = NULL, has_temporal_lag = FALSE) {


                     # hack to import compute_A function via main module
                     # reticulate::py_run_string(mbdmb_pycode)
                     # mbdmpy <<- reticulate::import_main()

                     ## Set data parameters
                     G <- length(y_list)  # number of outcomes
                     self$G <- G
                     self$N <- N
                     self$TT <- TT
                     self$NGT <- N*G*TT
                     self$NG <- N*G

                     ## Make y vector
                     y <- c()
                     for (t in 1:TT) {
                       start <- (t-1)*N + 1
                       end <- start + N - 1
                       for (k in 1:G) {
                         y <- c(y, y_list[[k]][start:end])
                       }
                     }
                     self$y <- y


                     ## Make X matrix, get (or make) beta names
                     self$X <- make_X(X_list = X_list, N = N, TT = TT)
                     self$K <- ncol(self$X)
                     beta_names <- unlist(lapply(X_list, colnames))
                     if (is.null(beta_names)) {
                       beta_names <- unlist(lapply(X_list, function(x) {
                         nc <- ncol(x)
                         return(paste0('beta_', 1:nc))
                       }))
                     }
                     self$beta_names <- beta_names

                     ## Create the period-wise weights matrices
                     omega_t_list <- list()
                     eta_t_names <- c() # Vector to store names of eta_t parameters

                     # Spatial weights
                     if (!is.null(W_t)) {
                       spatial_matrices <- make_spatial_weights(W_t, G)
                       omega_t_list <- c(omega_t_list, spatial_matrices)
                       eta_t_names <- c(eta_t_names, paste0('rho_', 1:G))
                     }

                     # Outcome-interdependence weights
                     # We create one weights-matrix per outcome
                     # Each of these maps the non-k outcomes onto the kth outcome
                     if (G > 1) {
                       outcome_matrices <- make_outcome_weights(N, G)
                       omega_t_list <- c(omega_t_list, outcome_matrices)
                       eta_t_names <- c(eta_t_names, paste0('lambda_', 1:G))
                     }

                     self$omega_t_list <- omega_t_list


                     ## Create the full weights matrices
                     omega_list <- list()
                     eta_names <- c()  # Vector to store names of eta parameters

                     # Create the temporal lag weights matrix
                     self$has_temporal_lag <- has_temporal_lag
                     if (has_temporal_lag) {
                       TM <- make_temporal_lag_matrix(N, G, TT)
                       omega_list[[1]] <- TM
                       eta_names <- c(eta_names, 'gamma')
                     }

                     # Create the block-diagonal spatial/outcome weights matrices
                     if (length(omega_t_list) > 0) {
                       for (i in 1:length(omega_t_list)) {
                         omega_list[[length(omega_list)+1]] <- bdiag(rep(list(omega_t_list[[i]]), TT))
                       }
                     }
                     self$omega_list <- omega_list
                     eta_names <- c(eta_names, eta_t_names)

                     ## Store the names of the weight (eta) parameters
                     self$eta_names <- eta_names
                     self$eta_t_names <- eta_t_names

                     ## Determine number of weight parameters (eta)
                     self$n_eta <- length(omega_list)
                     self$n_eta_t <- length(omega_t_list)

                     ## Determine whether we can use the Neumann approximation to compute d
                     ## If so, precompute diagonals
                     self$use_neumann <- (self$n_eta_t == 1)
                     L <- 8L  # We use the 8th order Neumann series to approximate
                     self$L <- L
                     if (self$use_neumann) {
                       neumann_list <- vector('list', L)
                       neumann_diag_list <- vector('list', L)
                       R <- omega_t_list[[1]]
                       neumann_list[[1]] <- R
                       neumann_diag_list[[1]] <- diag(R)
                       for (i in 2:L) {
                         neumann_list[[i]] <- neumann_list[[i-1]]%*%R
                         neumann_diag_list[[i]] <- diag(neumann_list[[i]])
                       }
                       self$neumann_diag_list <- neumann_diag_list
                       rm(neumann_list)
                     }

                     ## Determine location of beta/eta_t/eta in theta parameter vector
                     n_parameters <- (self$K + self$n_eta)
                     self$n_parameters <- n_parameters
                     self$beta_indices <- 1:self$K
                     if (self$n_eta > 0) {
                       self$eta_indices <- (self$K+1):(self$K + self$n_eta)
                     }
                     if (self$n_eta_t > 0) {
                       self$eta_t_indices <- (n_parameters - self$n_eta_t + 1):n_parameters
                     }
                   }

                 ))


#### Methods extracting specific parameter vectors from theta
MBDM$set("public", "extract_beta",
         function(theta) {
           beta <- theta[self$beta_indices]
           return(beta)
           })

MBDM$set("public", "extract_eta",
         function(theta) {
           eta <- theta[self$eta_indices]
           names(eta) <- NULL
           return(eta)
         })

MBDM$set("public", "extract_eta_t",
         function(theta) {
           eta_t <- theta[self$eta_t_indices]
           names(eta_t) <- NULL
           return(eta_t)
         })


#### Method for computing mu if there's a temporal lag
MBDM$set("public", "compute_mu_full",
         function(beta, eta) {

          ## Xbeta
          Xbeta <- self$X%*%beta

          ## A (the full weights matrix)
          A <- .sparseDiagonal(self$NGT)
          for (i in 1:length(eta)) {
            A <- A - self$omega_list[[i]]*eta[i]
          }

          ## Get mu
          mu <- scipy$sparse$linalg$bicgstab(r_to_py(A), Xbeta)[[1]]
          mu <- as.vector(mu)

          return(mu)
         })


#### Method for computing mu if there's NO temporal lag
MBDM$set("public", "compute_mu_blockwise",
          function(beta, eta_t) {

            ## Xbeta
            Xbeta <- self$X%*%beta

            ## A_t (the period-wise weights matrix)
            A_t <- .sparseDiagonal(self$NG)
            for (i in 1:length(eta_t)) {
              A_t <- A_t - self$omega_t_list[[i]]*eta_t[i]
            }

            ## Get each mu_t, concatenate
            mu <- c()
            for (t in 1:self$TT) {
              start <- (t-1)*self$NG + 1
              end <- start + self$NG - 1
              Xbeta_t <- Xbeta[start:end,,drop=FALSE]
              mu_t <- scipy$sparse$linalg$bicgstab(r_to_py(A_t), Xbeta_t)[[1]]
              mu <- c(mu, as.vector(mu_t))
            }

            return(mu)
          })


#### Method for computing mu if there are no dependencies
MBDM$set("public", "compute_mu_independent",
         function(beta) {

           Xbeta <- self$X%*%beta
           mu <- as.vector(Xbeta)

           return(mu)
         })


#### Compute d using Neumann approximation
MBDM$set("public", "compute_d_neumann",
          function(eta_t) {

            ## Multiply each diagonal with power of eta_t scalar
            eta_vec <- eta_t^(1:self$L)
            rw_list <- lapply(1:self$L, function(i) eta_vec[i] * self$neumann_diag_list[[i]])
            d_t <- rep(1, self$NG) + Reduce('+', rw_list)

            ## Full diagonal
            d <- rep(d_t, self$TT)

            return(d)
          })

#### Compute d using Takahashi equations
MBDM$set("public", "compute_d_takahashi",
        function(eta_t) {

          ## A_t (the period-wise weights matrix)
          omega_t_plus_list <- lapply(1:self$n_eta_t, function(k) {self$omega_t_list[[k]]*eta_t[k]})
          omega_t_sum <- Reduce("+", omega_t_plus_list)
          A_t <- .sparseDiagonal(self$NG) - omega_t_sum

          ## Get d_t using takahashi
          d_t <- diag(takahashi_inverse(A_t))

          ## Full diagonal
          d <- rep(d_t, self$TT)

          return(d)
        })


#### Compute likelihood
MBDM$set("public", "compute_llik",
          function(theta, verbose = FALSE, max_weights = 1) {

            ## Extract parameters from theta
            beta <- self$extract_beta(theta)
            eta_t <- self$extract_eta_t(theta)
            eta <- self$extract_eta(theta)

            ## Print the candidate eta if verbose
            if (verbose) {
              print(eta)
              flush.console()
            }

            ## Return NA if weights exceed max_weights
            if (self$n_eta > 0) {
              if (length(max_weights) < length(eta)) {
                max_weights <- rep(max_weights[1], length(eta))
              }
              if (any(abs(eta) > max_weights)) {
                return(NA)
              }
              if (max(abs(eta)) > 1) {
                return(NA)
              }
            }

            ## Compute mu
            if (self$n_eta > 0) {
              if (self$has_temporal_lag) {
                mu <- self$compute_mu_full(beta, eta)
              } else {
                mu <- self$compute_mu_blockwise(beta, eta_t)
              }
            } else {
                mu <- self$compute_mu_independent(beta)
            }

            ## Compute d
            if (self$n_eta_t > 0) {
              if (self$use_neumann) {
                d <- self$compute_d_neumann(eta_t)
              } else {
                d <- self$compute_d_takahashi(eta_t)
              }
            } else {
              d <- rep(1, self$NGT)
            }

            ## Compute likelihood
            v <- mu/d
            log_pr <- pnorm(v, log.p = TRUE, lower.tail = TRUE)
            log_not_pr <- pnorm(v, log.p = TRUE, lower.tail = FALSE)
            llik <- self$y*log_pr + (1-self$y)*log_not_pr

            return(sum(llik))
          })


#### Train model
MBDM$set("public", "train",
         function(theta_init, verbose = FALSE, max_weights = 1, method = "BFGS") {

           theta_names <- c(self$beta_names, self$eta_names)
           names(theta_init) <- theta_names

           ml.fit <- maxLik(logLik = self$compute_llik, start = theta_init, method = method,
                            verbose = verbose, max_weights = max_weights)
           return(ml.fit)

         })


##############################################################
# Data simulation function
##############################################################

simulate_data <- function(G, N, TT, K_outcome,
                          has_spatial_lag = FALSE, has_temporal_lag = FALSE,
                          beta = NULL, rho = NULL, lambda = NULL, gamma = NULL,
                          time_variant_X = TRUE) {

  ## Build X
  X_list <- vector('list', G)
  for (k in 1:G) {
    if (time_variant_X) {
      X_list[[k]] <- cbind(1, matrix(rnorm(N*TT*K_outcome), N*TT, K_outcome))
    } else {
      X_t <-  cbind(1, matrix(rnorm(N*K_outcome), N, K_outcome))
      X_list[[k]] <- do.call('rbind', rep(list(X_t), TT))
    }
  }
  X <- make_X(X_list = X_list, N = N, TT = TT)

  ## Make the period-wise weights matrices
  omega_t_list <- list()

  if (has_spatial_lag) {
    # Make weights matrix for one outcome
    W_t <- make_W_t(N)
    # Make weights matrices for all outcomes
    spatial_weights_matrices <- make_spatial_weights(W_t, G)
    omega_t_list <- c(omega_t_list, spatial_weights_matrices)
  } else {
    W_t <- NULL
  }

  if (G > 1) {
    outcome_interdep_matrices <- make_outcome_weights(N, G)
    omega_t_list <- c(omega_t_list, outcome_interdep_matrices)
  }

  ## Make the period-wise weight parameters
  eta_t <- c()
  if (has_spatial_lag) {
    if (is.null(rho)) {
      rho <- rep(0.25, G)
    }
    eta_t <- c(eta_t, rho)
  }
  if (G > 1) {
    if (is.null(lambda)) {
      lambda <- rep(0.25, G)
    }
    eta_t <- c(eta_t, lambda)
  }

  ## Make the temporal lag parameter
  if (has_temporal_lag) {
    if (is.null(gamma)) {
      gamma <- 0.25
    }
  }

  ## Make the period-wise inverse multiplier (A_t)
  At <- .sparseDiagonal(n = N*G)
  if (length(omega_t_list) > 0) {
    for (k in 1:length(omega_t_list)) {
      At <- At - omega_t_list[[k]]*eta_t[k]
    }
  }

  ## Make the systematic component
  if (is.null(beta)) {
    beta <- rep(c(-0.5, rep(1, K_outcome)), G)
  }
  Xbeta <- X%*%beta

  ## Make the latent outcome
  ystar_list <- vector('list', TT)
  # Period 1
  ystar_list[[1]]  <- solve(At, Xbeta[1:(N*G),] + rnorm(N*G))
  # Periods 2-TT
  if (TT > 1) {
    for (t in 2:TT) {
      start <- (t-1)*N*G + 1
      end <- start + N*G - 1
      if (has_temporal_lag) {
        ystar_list[[t]] <- solve(At, Xbeta[start:end,] + gamma*ystar_list[[t-1]] + rnorm(N*G))
      } else {
        ystar_list[[t]] <- solve(At, Xbeta[start:end,] + rnorm(N*G))
      }
    }
  }
  ystar <- as.vector(do.call(rbind, ystar_list))

  ## Compute discrete outcome
  y <- (ystar>0)*1

  ## Separate y into list of N*TT outcomes
  y_list <- vector('list', G)
  for (t in 1:TT) {
    for (k in 1:G) {
      start <- (t-1)*N*G + (k-1)*N + 1
      end <- start + N - 1
      y_list[[k]] <-  c(y_list[[k]], y[start:end])
    }
  }

  out.ls <- list(y=y, X=X, y_list=y_list, X_list=X_list,
                 W_t = W_t, At = At, N = N, TT = TT, G = G)
}


