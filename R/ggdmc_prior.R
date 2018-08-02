##' Specifying Parameter Prior Distributions
##'
##' \code{BuildPrior} sets up parameter prior distributions for each model
##' parameter. \code{p1} and \code{p2} refer to the first and second parameters
##' a prior distribution.
##'
##' Four distribution types are implemented:
##' \enumerate{
##' \item Normal and truncated normal, where: p1 = mean, p2 = sd. It specifies
##' a normal distribution when bounds are set -Inf and Inf,
##' \item Beta, where: p1 = shape1 and p2 = shape2 (see \link{pbeta}). Note the
##'       uniform distribution is a special case of the beta with p1 and
##'       p2 = 1),
##' \item Gamma, where p1 = shape and p2 = scale (see \link{pgamma}). Note p2 is
##'       scale, not rate,
##' \item Lognormal, where p1 = meanlog and p2 = sdlog (see \link{plnorm}).
##' }
##'
##' @param dists a vector of character string specifying a distribution.
##' @param p1 the first parameter of a distribution
##' @param p2 the second parameter of a distribution
##' @param lower lower support (boundary)
##' @param upper upper support (boundary)
##' @param untrans whether to do log transformation. Default is not
##' @param types available distribution types
##' @return a list of list
##' @export
BuildPrior <- function(p1, p2,
  lower   = rep(NA, length(p1)),
  upper   = rep(NA, length(p1)),
  dists   = rep("tnorm", length(p1)),
  untrans = rep("identity", length(p1)),
  types   = c("tnorm", "beta", "gamma", "lnorm", "cauchy", "constant")) {

  dist <- c()
  np1 <- length(p1)

  if (length(p2) == 1) p2 <- rep(p2, np1)
  if ( np1 != length(p2) )    stop("p1 and p2 must have the same length")
  if ( np1 != length(lower) ) stop("p1 and lower must have the same length")
  if ( np1 != length(upper) ) stop("p1 and upper must have the same length")
  both.not.na <- !is.na(upper) & !is.na(lower)

  if ( any(upper[both.not.na] <= lower[both.not.na]) )
    stop("All elements of upper must be greater than lower")
  if ( np1 != length(dists) ) stop("p1 and dists must have the same length")
  if ( !all(dists %in% types) )
    stop(paste("Unsupported distribution, allowable types are:",
      paste(types, collapse = ", ")))
  name.untrans <- length(untrans) != np1
  if (name.untrans & (is.null(names(untrans)) | is.null(names(p1))))
    stop("If untrans vector is not the same length as p1 it must have p1 names")
  if (!(all(names(untrans) %in% names(p1)))) stop("untrans vector has names not in p1 names")

  prior <- vector(mode = "list", length = np1)
  names(prior) <- names(p1)
  for (i in 1:np1) {
    prior[[i]] <- switch(dists[i],
      tnorm = {
        if (is.na(lower[i])) lower[i] <- -Inf
        if (is.na(upper[i])) upper[i] <- Inf
        p <- c(p1[i], p2[i], lower[i], upper[i])
        names(p) <- c("mean","sd","lower","upper")
        p <- as.list(p)
        attr(p,"dist") <- "tnorm"
        p
      },
      beta = {
        if (is.na(lower[i])) lower[i] <- 0
        if (is.na(upper[i])) upper[i] <- 1
        p <- c(p1[i],p2[i],lower[i],upper[i])
        names(p) <- c("shape1","shape2", "lower","upper")
        p <- as.list(p)
        attr(p, "dist") <- "beta_lu"
        p
      },
      gamma={
        if (is.na(lower[i])) lower[i] <- 0
        if (is.na(upper[i])) upper[i] <- Inf
        p <- c(p1[i],p2[i],lower[i], upper[i])
        names(p) <- c("shape", "scale","lower", "upper")
        p <- as.list(p)
        attr(p, "dist") <- "gamma_l"
        p
      },
      lnorm={
        if (is.na(lower[i])) lower[i] <- 0
        if (is.na(upper[i])) upper[i] <- Inf
        p <- c(p1[i], p2[i], lower[i], upper[i])
        names(p) <- c("meanlog","sdlog","lower", "upper")
        p <- as.list(p)
        attr(p, "dist") <- "lnorm_l"
        p
      },
      {
        p <- c(p1[i], 0, -Inf, Inf)
        names(p) <- c("constant", "sd", "lower", "upper")
        p <- as.list(p)
        attr(p, "dist") <- "constant"
        p
      }
    )
    prior[[i]]$log <- TRUE
    if (!name.untrans) attr(prior[[i]], "untrans") <- untrans[i] else
      if (is.na(untrans[names(p1)[i]]))
        attr(prior[[i]], "untrans") <- "identity" else
          attr(prior[[i]],"untrans") <- untrans[names(p1)[i]]
  }

  class(prior) <- "prior"
  return(prior)
}

##' Parameter Prior Distributions
##'
##' Probability density functions and random generation for parameter prior
##' distributions. \code{rprior} is a faster compatible function, like
##' \code{rprior.dmc} in \code{DMC}. It is a wrapper function for
##' \code{rprior_scalar} and \code{rprior_mat} functions.
##' \code{rprior_scalar} matches five strings:
##' \code{tnorm}, \code{beta_lu}, \code{gamma_l}, \code{lnorm_l}, and
##' \code{constant} to select a random functions.
##'
##' @param p.prior a list of list usually created by BuildPrior to store the
##' information about parameter prior distributions.
##' @param n number of observations/random draws
##' @param pvec a parameter vector
##' @param dists dists argument in BuildPrior
##' @param p1 p1 argument in prior.p.dmc
##' @param p2 p2 argument in prior.p.dmc
##' @param lower lower argument in prior.p.dmc
##' @param upper upper argument in prior.p.dmc
##' @return \code{rprior_vec} gives a column vector; \code{rprior} gives a
##' matrix; \code{rprior} gives a \code{n x npar} named numeric
##' matrix; \code{rprior_scalar} gives a named numeric vector.
##' \code{sum_log_prior} gives prior likelihood.
##'
##' @examples
##' p.prior <- ggdmc::BuildPrior(
##'  dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta"),
##'  p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
##'  p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
##'  lower = c(0,-5, NA, NA, 0, NA),
##'  upper = c(2, 5, NA, NA, 2, NA))
##'
##' rprior(p.prior, 9)
##' ##               a           v         z         sz        sv         t0
##' ## [1,] 0.97413686  0.78446178 0.9975199 -0.5264946 0.5364492 0.55415052
##' ## [2,] 0.72870190  0.97151662 0.8516604  1.6008591 0.3399731 0.96520848
##' ## [3,] 1.63153685  1.96586939 0.9260939  0.7041254 0.4138329 0.78367440
##' ## [4,] 1.55866180  1.43657110 0.6152371  0.1290078 0.2957604 0.23027759
##' ## [5,] 1.32520281 -0.07328408 0.2051155  2.4040387 0.9663111 0.06127237
##' ## [6,] 0.49628528 -0.19374770 0.5142829  2.1452972 0.4335482 0.38410626
##' ## [7,] 0.03655549  0.77223432 0.1739831  1.4431507 0.6257398 0.63228368
##' ## [8,] 0.71197612 -1.15798082 0.8265523  0.3813370 0.4465184 0.23955415
##' ## [9,] 0.38049166  3.32132034 0.9888108  0.9684292 0.8437480 0.13502154
##'
##' pvec <- c(a=1, v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
##' p.prior  <- BuildPrior(
##'   dists = rep("tnorm", 6),
##'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
##'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05) * 5,
##'   lower = c(0,-5, 0, 0, 0, 0),
##'   upper = c(5, 7, 2, 2, 2, 2))
##'
##' ## summed_log_prior(pvec, p.prior)
##'
##' \dontrun{
##' setwd("/home/yslin/Documents/dmc-amsterdam17/")
##' source("dmc/dmc.R")
##' load_model("LBA", "lba_B.R")
##' require(microbenchmark)
##' res <- microbenchmark(ggdmc::rprior(prior1, 6),
##'                       rprior.dmc(prior1, 6),  times=1e3)
##' }
##'
##' ## Unit: microseconds
##' ##                     expr         min       lq      mean   median
##' ## ggdmc::rprior(prior1, 6)      29.964  33.7350  39.92127  42.4295
##' ## rprior.dmc(prior1, 6)        307.585 343.8335 372.32580 359.2680
##' ##       uq      max neval
##' ##  44.7700   98.897  1000
##' ## 376.3445 2446.987  1000
##' @export
rprior <- function(p.prior, n = 1) {
  if(n == 1) { out <- rprior_scalar(p.prior) }
  else       { out <- rprior_mat(p.prior, n) }
  return(out)
}


# dprior <- function(pvec, prior) {
#
#   dists <- c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta")
#   dists_internal <- c("tnorm", "tnorm", "beta_lu", "tnorm", "beta_lu", "beta_lu")
#   p1    <- c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1)
#   p2    <- c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1)
#   lower <- c(0, -5, NA, NA, 0, NA)
#   upper <- c(2,  5, NA, NA, 2, NA)
#   islog <- rep(1, 6)
#   pvec  <- c(a = 1.15, v = -.10, z = .74, sz = 1.23, sv = .11, t0 = .87)
#   pnames <- names(pvec); pnames
#   p.prior <- ggdmc::BuildPrior(dists = dists, p1 = p1, p2 = p2, lower = lower,
#     upper = upper)
#
# }
