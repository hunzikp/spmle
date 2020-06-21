library(spmle)
set.seed(0)


#### Cross-sectional probit, no interdependence

## Set data parameters
N <- 1024
TT <- 1
K_outcome <- 1
has_temporal_lag <- FALSE
has_spatial_lag <- FALSE
zero_period <- FALSE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
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



#### Cross-sectional spatial probit

## Set data parameters
N <- 1024
TT <- 1
K_outcome <- 1
has_temporal_lag <- FALSE
has_spatial_lag <- TRUE
zero_period <- FALSE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
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



#### Temporal autoregressive Probit

## Set data parameters
N <- 1024
TT <- 10
K_outcome <- 1
has_temporal_lag <- TRUE
has_spatial_lag <- FALSE
zero_period <- TRUE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
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


#### Spatio-temporal Probit

## Set data parameters
N <- 1024
TT <- 10
K_outcome <- 1
has_temporal_lag <- TRUE
has_spatial_lag <- TRUE
zero_period <- TRUE

## Make the data
data_ls <- simulate_data(N = N, TT = TT,
                         K_outcome = K_outcome,
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
