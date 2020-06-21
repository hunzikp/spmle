###############################################################################
##### Calabrese & Elkink's Recursive Importance Sampler Implementation

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
