mbdm Tests
================
Philipp Hunziker
November 25, 2018

Set seed, load package
======================

Tests with Simulated Data
=========================

Cross-sectional probit, no interdependence
------------------------------------------

``` r
## Set data parameters
G <- 1
N <- 10000
TT <- 1
K_outcome <- 1
has_temporal_lag = FALSE
has_spatial_lag = FALSE

## Make the data
data.ls <- simulate_data(G, N, TT, K_outcome, has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                         rho = 0.5, beta = rep(c(0,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
```

    ##    user  system elapsed 
    ##   0.167   0.005   0.171

``` r
## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 34 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -4993.476 
    ## 2  free parameters
    ## Estimates:
    ##        Estimate Std. error t value Pr(> t)    
    ## beta_1  0.00831    0.01444   0.576   0.565    
    ## beta_2  1.00262    0.01925  52.074  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

Cross-sectional spatial probit
------------------------------

``` r
## Set data parameters
G <- 1
N <- 10000
TT <- 1
K_outcome <- 1
has_temporal_lag = FALSE
has_spatial_lag = TRUE

## Make the data
data.ls <- simulate_data(G, N, TT, K_outcome, has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                         rho = 0.5, beta = rep(c(0,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
```

    ##    user  system elapsed 
    ##   4.896   1.473   4.803

``` r
## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 63 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -4932.498 
    ## 3  free parameters
    ## Estimates:
    ##        Estimate Std. error t value Pr(> t)    
    ## beta_1 0.006740   0.006996   0.963   0.335    
    ## beta_2 0.980927   0.019092  51.380  <2e-16 ***
    ## rho_1  0.544023   0.019944  27.278  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

Panel spatial probit without temporal autoregression
----------------------------------------------------

``` r
## Set data parameters
G <- 1
N <- 1024
TT <- 10
K_outcome <- 1
has_temporal_lag = FALSE
has_spatial_lag = TRUE

## Make the data
data.ls <- simulate_data(G, N, TT, K_outcome, has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                         rho = 0.5, beta = rep(c(0,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
```

    ##    user  system elapsed 
    ##  10.833   0.000  10.834

``` r
## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 111 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -5069.155 
    ## 3  free parameters
    ## Estimates:
    ##         Estimate Std. error t value Pr(> t)    
    ## beta_1 -0.003129   0.007433  -0.421   0.674    
    ## beta_2  0.990081   0.018928  52.309  <2e-16 ***
    ## rho_1   0.504099   0.020831  24.199  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

Panel spatio-temporal probit
----------------------------

``` r
## Set data parameters
G <- 1
N <- 1024
TT <- 10
K_outcome <- 1
has_temporal_lag = TRUE
has_spatial_lag = TRUE

## Make the data
data.ls <- simulate_data(G, N, TT, K_outcome, has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                         rho = 0.5, gamma = 0.25, beta = rep(c(0,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
```

    ##    user  system elapsed 
    ##   8.779   0.008   8.792

``` r
## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 111 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -5061.396 
    ## 4  free parameters
    ## Estimates:
    ##        Estimate Std. error t value Pr(> t)    
    ## beta_1 0.001498   0.004352   0.344   0.731    
    ## beta_2 0.952888   0.018486  51.548  <2e-16 ***
    ## gamma  0.232480   0.011735  19.811  <2e-16 ***
    ## rho_1  0.503298   0.017667  28.488  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

Multi-outcome cross-section with spatial-lag
--------------------------------------------

``` r
## Set data parameters
G <- 2
N <- 1024
TT <- 1
K_outcome <- 1
has_temporal_lag = FALSE
has_spatial_lag = TRUE

## Make the data
data.ls <- simulate_data(G, N, TT, K_outcome, has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                         rho = c(0.5, 0.5), lambda = c(0.25,0.25), beta = rep(c(0,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
```

    ##    user  system elapsed 
    ##  94.840   1.181  96.024

``` r
## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 113 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -1030.6 
    ## 8  free parameters
    ## Estimates:
    ##          Estimate Std. error t value  Pr(> t)    
    ## beta_1   -0.05631    0.03538  -1.591    0.112    
    ## beta_2    0.86995    0.05548  15.681  < 2e-16 ***
    ## beta_1   -0.02870    0.02908  -0.987    0.324    
    ## beta_2    0.89506    0.05651  15.838  < 2e-16 ***
    ## rho_1     0.40591    0.06621   6.131 8.75e-10 ***
    ## rho_2     0.50498    0.06003   8.412  < 2e-16 ***
    ## lambda_1  0.29720    0.04941   6.015 1.80e-09 ***
    ## lambda_2  0.22532    0.04774   4.720 2.36e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

Multi-outcome panel with spatio-temporal lag (aka 'the full monty')
-------------------------------------------------------------------

``` r
## Set data parameters
G <- 2
N <- 256
TT <- 10
K_outcome <- 1
has_temporal_lag = TRUE
has_spatial_lag = TRUE

## Make the data
data.ls <- simulate_data(G, N, TT, K_outcome, has_temporal_lag = has_temporal_lag, has_spatial_lag = has_spatial_lag,
                         rho = c(0.25, 0.25), lambda = c(0.25,0.25), gamma = 0.25, beta = rep(c(0,1), G))

## Make the mbdm object
mbdm <- MBDM$new(y_list = data.ls$y_list, X_list = data.ls$X_list, N = N, TT = TT,
                 W_t = data.ls$W_t, has_temporal_lag = has_temporal_lag)

## ML fit
theta_init <- runif(mbdm$n_parameters, min = -0.75, max = 0.75)
system.time(ml.fit <- mbdm$train(theta_init))
```

    ##    user  system elapsed 
    ##  57.944   0.000  57.944

``` r
## Summary
summary(ml.fit)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 135 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -2511.704 
    ## 9  free parameters
    ## Estimates:
    ##           Estimate Std. error t value  Pr(> t)    
    ## beta_1   1.795e-02  1.724e-02   1.041    0.298    
    ## beta_2   9.565e-01  3.699e-02  25.855  < 2e-16 ***
    ## beta_1   2.973e-05  2.072e-02   0.001    0.999    
    ## beta_2   9.156e-01  3.626e-02  25.252  < 2e-16 ***
    ## gamma    2.445e-01  1.640e-02  14.910  < 2e-16 ***
    ## rho_1    2.857e-01  4.207e-02   6.791 1.11e-11 ***
    ## rho_2    1.877e-01  4.721e-02   3.976 7.00e-05 ***
    ## lambda_1 2.298e-01  2.846e-02   8.074 6.80e-16 ***
    ## lambda_2 2.901e-01  2.940e-02   9.865  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------
