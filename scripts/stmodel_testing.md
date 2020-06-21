spmle Examples
================
Philipp Hunziker
June 21, 2020

# Usage Examples

``` r
library(spmle)
```

## Cross-sectional Probit

``` r
## Set data parameters
N <- 1024
TT <- 1
K_outcome <- 1
beta <- c(-0.5, 1)
has_temporal_lag <- FALSE
has_spatial_lag <- FALSE
zero_period <- FALSE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
                         beta = beta,
                         has_spatial_lag = has_spatial_lag,
                         has_temporal_lag = has_temporal_lag,
                         zero_period = zero_period)

## Make the model object
model <- STModel$new(y = data_ls$y, X = data_ls$X,
                     N = data_ls$N, TT = data_ls$TT,
                     W_t = data_ls$W_t,
                     has_temporal_lag = has_temporal_lag,
                     has_spatial_lag = has_spatial_lag,
                     zero_period = zero_period)

## ML fit
ml.fit <- model$autotrain()

## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 5 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -490.3983 
    ## 2  free parameters
    ## Estimates:
    ##        Estimate Std. error t value Pr(> t)    
    ## beta_1 -0.51169    0.04818  -10.62  <2e-16 ***
    ## beta_2  0.93440    0.05752   16.24  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

## Cross-sectional Spatial Probit

``` r
## Set data parameters
N <- 1024
TT <- 1
K_outcome <- 1
beta <- c(-0.5, 1)
rho <- 0.25
has_temporal_lag <- FALSE
has_spatial_lag <- TRUE
zero_period <- FALSE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
                         beta = beta,
                         rho = rho,
                         has_spatial_lag = has_spatial_lag,
                         has_temporal_lag = has_temporal_lag,
                         zero_period = zero_period)

## Make the model object
model <- STModel$new(y = data_ls$y, X = data_ls$X,
                     N = data_ls$N, TT = data_ls$TT,
                     W_t = data_ls$W_t,
                     has_temporal_lag = has_temporal_lag,
                     has_spatial_lag = has_spatial_lag,
                     zero_period = zero_period)

## ML fit
ml.fit <- model$autotrain()

## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 45 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -432.952 
    ## 3  free parameters
    ## Estimates:
    ##        Estimate Std. error t value  Pr(> t)    
    ## beta_1 -0.55301    0.06783  -8.152 3.57e-16 ***
    ## beta_2  1.15860    0.07163  16.175  < 2e-16 ***
    ## rho     0.32148    0.08484   3.789 0.000151 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

## Temporal Autoregressive Probit

``` r
## Set data parameters
N <- 1024
TT <- 10
K_outcome <- 1
beta <- c(-0.5, 1)
gamma <- 0.25
has_temporal_lag <- TRUE
has_spatial_lag <- FALSE
zero_period <- TRUE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
                         beta = beta, 
                         gamma = gamma,
                         has_spatial_lag = has_spatial_lag,
                         has_temporal_lag = has_temporal_lag,
                         zero_period = zero_period)

## Make the model object
model <- STModel$new(y = data_ls$y, X = data_ls$X,
                     N = data_ls$N, TT = data_ls$TT,
                     W_t = data_ls$W_t,
                     has_temporal_lag = has_temporal_lag,
                     has_spatial_lag = has_spatial_lag,
                     zero_period = zero_period)

## ML fit
ml.fit <- model$autotrain()

## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 38 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -4699.288 
    ## 3  free parameters
    ## Estimates:
    ##        Estimate Std. error t value Pr(> t)    
    ## beta_1 -0.45056    0.01448  -31.12  <2e-16 ***
    ## beta_2  0.93315    0.01909   48.88  <2e-16 ***
    ## gamma   0.27258    0.01436   18.98  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

## Spatiotemporal Probit

Note how here we trigger the Python speed-up by setting `use_py = TRUE`
in the model instantiation.

``` r
## Set data parameters
N <- 1024
TT <- 10
K_outcome <- 1
beta <- c(-0.5, 1)
rho <- 0.25
gamma <- 0.25
has_temporal_lag <- TRUE
has_spatial_lag <- TRUE
zero_period <- TRUE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
                         beta = beta, 
                         gamma = gamma,
                         rho = rho,
                         has_spatial_lag = has_spatial_lag,
                         has_temporal_lag = has_temporal_lag,
                         zero_period = zero_period)

## Make the model object
model <- STModel$new(y = data_ls$y, X = data_ls$X,
                     N = data_ls$N, TT = data_ls$TT,
                     W_t = data_ls$W_t,
                     has_temporal_lag = has_temporal_lag,
                     has_spatial_lag = has_spatial_lag,
                     zero_period = zero_period,
                     use_py = TRUE)  # use python speed-up!

## ML fit
ml.fit <- model$autotrain()

## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 38 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -4074.403 
    ## 4  free parameters
    ## Estimates:
    ##        Estimate Std. error t value  Pr(> t)    
    ## beta_1 -0.52893    0.03060 -17.288  < 2e-16 ***
    ## beta_2  0.97993    0.02109  46.461  < 2e-16 ***
    ## rho     0.22241    0.02992   7.433 1.06e-13 ***
    ## gamma   0.23573    0.01495  15.771  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------
