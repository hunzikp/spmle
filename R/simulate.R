##############################################################
# Simulation functions
##############################################################

make_W_t <- function(N) {
  # N must have natural root

  ras <- raster::raster(matrix(1, sqrt(N), sqrt(N)))
  spdf <- raster::rasterToPolygons(ras)
  B <- gIntersects(spdf, byid = TRUE)
  diag(B) <- FALSE
  B <- B*1
  W_t <- B / rowSums(B)
  W_t <- Matrix(W_t, sparse = TRUE)

  return(W_t)
}

simulate_data <- function(N, TT, K_outcome,
                          has_spatial_lag = FALSE, has_temporal_lag = FALSE,
                          beta = NULL, rho = NULL, gamma = NULL,
                          time_variant_X = TRUE, zero_period = TRUE) {

  ## Build X
  if (time_variant_X) {
    X <- cbind(1, matrix(rnorm(N*TT*K_outcome), N*TT, K_outcome))
  } else {
    X_t <-  cbind(1, matrix(rnorm(N*K_outcome), N, K_outcome))
    X <- do.call('rbind', rep(list(X_t), TT))
  }

  ## Make spatial components
  if (has_spatial_lag) {
    W_t <- make_W_t(N)
    if (is.null(rho)) {
      rho <- 0.25
    }
  } else {
    W_t <- NULL
  }

  if (has_temporal_lag) {
    if (is.null(gamma)) {
      gamma <- 0.25
    }
  }

  ## Make the period-wise inverse multiplier (A_t)
  At <- .sparseDiagonal(n = N)
  if (has_spatial_lag) {
    At <- At - rho*W_t
  }

  ## Make the systematic component
  if (is.null(beta)) {
    beta <- c(-0.5, rep(1, K_outcome))
  }
  Xbeta <- X%*%beta

  ## Compute latent outcome
  ystar_list <- vector('list', TT)

  ## Make first period
  if (zero_period & has_temporal_lag) {

    # Compute E(ystar | W)
    X_t_list <- vector('list', TT)
    for (t in 1:TT) {
      start <- (t-1)*N + 1
      end <- start + N - 1
      X_t_list[[t]] <- X[start:end,]
    }
    X_expectation <- Reduce('+', X_t_list)/TT

    M <- .sparseDiagonal(n = N) - .sparseDiagonal(n = N)*gamma
    if (has_spatial_lag) {
      M <- M - rho*W_t
    }

    ystar_0 <- as.vector(solve(M, X_expectation%*%beta))  # E(ystar | W)

    # First period outcome
    ystar_list[[1]]  <- solve(At, Xbeta[1:N,] + gamma*ystar_0 + rnorm(N))

  } else {
    ystar_list[[1]]  <- solve(At, Xbeta[1:N,] + rnorm(N))
  }

  # Periods 2-TT
  if (TT > 1) {
    for (t in 2:TT) {
      start <- (t-1)*N + 1
      end <- start + N - 1
      if (has_temporal_lag) {
        ystar_list[[t]] <- solve(At, Xbeta[start:end,] + gamma*ystar_list[[t-1]] + rnorm(N))
      } else {
        ystar_list[[t]] <- solve(At, Xbeta[start:end,] + rnorm(N))
      }
    }
  }
  ystar <- as.vector(do.call(rbind, ystar_list))

  ## Compute discrete outcome
  y <- (ystar>0)*1

  out.ls <- list(y=y, X=X, W_t = W_t, At = At, N = N, TT = TT)
}
