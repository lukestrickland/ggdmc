% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/ggdmc_density.R,
%   R/ggdmc_random.R
\name{rlnr}
\alias{rlnr}
\alias{n1PDFfixedt0_lnr}
\alias{rlnrR}
\title{Generate Random Choice-RT Data from LNR Model}
\usage{
rlnr(n, meanlog, sdlog, t0, st0 = 0)

n1PDFfixedt0_lnr(dt, meanlog, sdlog)

rlnrR(n, meanlog, sdlog, t0, st0 = 0)
}
\arguments{
\item{n}{numbers of observation}

\item{t0}{a scalar, a vector of length number of accumulators or a
matrix with 1 row per accumulator, when start time differ on each trial}

\item{st0}{range of non-decision time variability, must be a scalar, as the same
variability is assumed in a common encoding/production stage}

\item{mean_log}{a n_acc length vector or a n_acc x n matrix. mean
of the distribution on the log scale without default value}

\item{sd_log}{a n_acc length vector or a n_acc x n matrix. Standard
deviation of the distribution on the log scale without default value.}
}
\value{
a matrix
}
\description{
Race among \code{nacc} accumulators, using log-normal race model
}
\examples{
## A simple demo
pmat <- matrix(c(-1, 0, 1, 1, .2, .2, 0, 0), 2)
set.seed(123)
dat0 <- rlnr(4, pmat[,1], pmat[,2], pmat[,3], pmat[1,4])
##           [,1] [,2]
## [1,] 0.4100361    0
## [2,] 0.4922407    0
## [3,] 1.7855260    1
## [4,] 0.4822220    1
##
## Three accumulators
n <- 1e5
meanlog <- c(.5, .75, 1);
sdlog <- c(1,1,1)
t0 <- c(.2,1,1)
set.seed(123)
dat1 <- rlnr(n, meanlog, sdlog, t0)
table(dat1[,2])
hist(dat1[,1], breaks = "fd")

## t0 has one element only
t0 <- .2
set.seed(123)
dat2 <- rlnr(n, meanlog, sdlog, t0)
table(dat2[,2])
hist(dat2[,1], breaks = "fd", freq = FALSE)
## check t0 noise
st0 <- 1
set.seed(123)
dat3 <- rlnr(n, meanlog, sdlog, t0, st0)
table(dat3[,2])
hist(dat3[,1], breaks = "fd", freq = FALSE)

## Check matrix forms work (not yet)
dat4 <- rlnr(n, meanlog, sdlog, t0, st0)
sim <- data.frame(RT = dat4[,1], R = dat4[,2])
dns <- ggdmc:::plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
ichar <- as.character(1)

meanlog <- matrix(rep(meanlog,length(dns[[ichar]]$x)),nrow=n_acc)
sdlog   <- matrix(rep(sdlog,  length(dns[[ichar]]$x)),nrow=n_acc)
t0      <- matrix(rep(t0, length.out=n_acc*length(dns[[ichar]]$x)), nrow=n_acc)

## Almost 8 times faster
library(rbenchmark)
res <- benchmark(r1 = rlnr(n, c(-1,0), c(1,1), c(.2, .2), 0),
                 r2 = ggdmc:::rlnrDF(n, c(-1,0), c(1,1), c(.2, .2), 0),
                 r3 = ggdmc:::rlnr(n, c(-1,0), c(1,1), c(.2, .2), 0),
                 replications = 1000)
print(res[,1:4])

## test replications elapsed relative
##   r1         1000   1.695    7.995
##   r2         1000   0.380    1.792
##   r3         1000   0.212    1.000

}
