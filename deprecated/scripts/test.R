#### Calabrese and Elkink MC setup
### 2.1 Parameter settings
## Path variables
work <- "work/"
data.dir <- paste(work, "data/", sep = "")
libs.dir <- paste(work, "Rlibs/", sep = "")
output.dir <- paste(work, "output/", sep = "")
tmp.dir <- paste(work, "tmp/", sep = "")

## Monte Carlo options
mu.x <- 2
sd.x <- 4
b0 <- 4
b1 <- -2
sd.e <- 1
rho <- c(0, .1, .45, .8)
N <- c(50, 500)
## Prepare parameters matrix
params <- expand.grid(rho, N)
colnames(params) <- c("rho", "N")
params <- as.data.frame(params)
params$mu.x <- mu.x
params$sd.x <- sd.x
params$b0 <- b0
params$b1 <- b1
params$sd.e <- sd.e
params$id <- 1:dim(params)[1]

write.csv(params, file = paste(output.dir, "parameters.csv", sep = ""), row.names = FALSE)


### 2.2 Data generation
## Environment variables
methods <- c(Sys.getenv("MC_METHOD"))
sample.size <- as.integer(Sys.getenv("MC_N"))
subset.start <- as.integer(Sys.getenv("MC_START"))
subset.end <- as.integer(Sys.getenv("MC_END"))

process.id <- as.integer(Sys.getenv("MC_ID"))

## Path variables
work <- "work/"
data.dir <- paste(work, "data/", sep = "")
output.dir <- paste(work, "output/", sep = "")
tmp.dir <- paste(work, "tmp/", sep = "")

## Monte Carlo options
params <- read.csv(paste(output.dir, "parameters.csv", sep = ""))
lbl <- sprintf("%s_%s_%d", Sys.info()["nodename"], format(Sys.time(), "%Y%m%d_%H%M%S"), process.id)

## Libraries required
library(rlecuyer)
## Function to generate random W
## Follows description of algorithm by Beron & Vijverberg (2004)
generate.W <- function(ncases)
{
  x <- runif(ncases)
  y <- runif(ncases)
  if (ncases == 50) d <- .21
  else if (ncases == 500) d <- .06
  else if (ncases == 1500) d <- .036
  else stop("Do not know d to generate W!")
  D <- as.matrix(dist(cbind(x, y)))
  W <- D < d
  diag(W) <- 0
  rs <- rowSums(W)
  W[rs > 0, rs > 0] <- W[rs > 0, rs > 0] / rs[rs > 0]
  W
}

## Main monte carlo loop
streamName <- paste("str", process.id, sep = "")
.lec.CreateStream(streamName)
.lec.SetSeed(streamName, sample(1:60000, 6, replace = TRUE) * process.id)
.lec.CurrentStream(streamName)
p <- params[params$N == sample.size, ]
for (s in subset.start:subset.end) {
  for (i in 1:dim(p)[1]) {
    ## Generate exogenous variables
    x <- rnorm(p[i,"N"], p[i,"mu.x"], p[i,"sd.x"])
    Xb <- cbind(1, x) %*% c(p[i,"b0"], p[i,"b1"])
    e <- rnorm(p[i,"N"], 0, p[i,"sd.e"])
    W <- generate.W(p[i,"N"])
    ## Generate Y
    A <- -p[i,"rho"] * W
    diag(A) <- 1
    A <- solve(A)
    y.star <- A %*% Xb + A %*% e
    y <- ifelse(y.star > 0, 1, 0)
    ## Save the data
    save(y, y.star, x, e, W, p, s, i, file =
           sprintf("%sdata_%05d_%04d.Rdata", data.dir, s, p[i, "id"]))
  }
}
.lec.CurrentStreamEnd()


sample.size=50
i=1
s=1
## Main monte carlo loop
# streamName <- paste("str", process.id, sep = "")
# .lec.CreateStream(streamName)
# .lec.SetSeed(streamName, sample(1:60000, 6, replace = TRUE) * process.id)
# .lec.CurrentStream(streamName)
p <- params[params$N == sample.size, ]
for (s in subset.start:subset.end) {
  for (i in 1:dim(p)[1]) {
    ## Generate exogenous variables
    x <- rnorm(p[i,"N"], p[i,"mu.x"], p[i,"sd.x"])
    Xb <- cbind(1, x) %*% c(p[i,"b0"], p[i,"b1"])
    e <- rnorm(p[i,"N"], 0, p[i,"sd.e"])
    W <- generate.W(p[i,"N"])
    ## Generate Y
    A <- -p[i,"rho"] * W
    diag(A) <- 1
    A <- solve(A)
    y.star <- A %*% Xb + A %*% e
    y <- ifelse(y.star > 0, 1, 0)
    ## Save the data
    save(y, y.star, x, e, W, p, s, i, file =
           sprintf("%sdata_%05d_%04d.Rdata", data.dir, s, p[i, "id"]))
  }
}

y
x
W

rm(sample.size, i,s)
