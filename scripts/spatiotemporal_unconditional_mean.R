
library(spmle)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(42)


# Make some data
N <- 9
TT <- 100

W <- spmle::make_W_t(N)
X_means <- runif(N)*3
X_ls <- lapply(1:TT, function(i) {cbind(1, matrix(rnorm(N, X_means, 1), N, 1))})
X_expectation <- Reduce('+', X_ls)/TT

# Set parameters
beta <- c(-0.5, 1)
rho <- 0.25
lambda <- 0.25

# Compute unconditional outcome mean
# E(y*_t) = (I – hI – pW)^(-1) E(X)beta
Y_unconditional <- solve(diag(N) - diag(N)*lambda - rho*W)%*%X_expectation%*%beta

# Simulate process
ymat <- matrix(NA, TT, N)
ymat[1,] <- rep(0, N)  # initial Y
M <- solve(diag(N) - rho*W)
for (tt in 2:TT) {
  ymat[tt,] <- as.vector(M%*%(X_ls[[tt]]%*%beta + lambda*ymat[tt-1,]))
}

# Plot the mofo
ysim_df <- as.data.frame(ymat)
outcome_names <- paste0('outcome', 1:N)
names(ysim_df) <- outcome_names
ysim_df$period <- 1:TT
plot_tb <- ysim_df %>% pivot_longer(-period, names_to = "outcome")
yuc_tb <- tibble(outcome = outcome_names, value = as.vector(Y_unconditional))

ggplot(plot_tb, aes(period, value)) +
  geom_line() +
  geom_hline(data = yuc_tb, aes(yintercept = value), color = 'red') +
  facet_wrap(~outcome)
