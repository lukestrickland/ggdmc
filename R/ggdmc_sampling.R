### Prior Distributions -------------------------------------------------------
#' @importFrom stats dbeta
dbeta_lu <- function(x,shape1,shape2,lower,upper,log=FALSE) {
  # Used with beta prior
  if (log) {dbeta((x-lower)/(upper-lower),shape1,shape2,log=log)}
  else {dbeta((x-lower)/(upper-lower),shape1,shape2,log=log)/(upper-lower)}
}

#' @importFrom stats rbeta
rbeta_lu <- function(n,shape1,shape2,lower,upper) {
  # Used with beta prior
  lower + rbeta(n,shape1,shape2)*(upper-lower)
}

#' @importFrom stats dgamma
dgamma_l <- function(x,shape,scale,lower,log=FALSE) {
  # Used with gamma prior
  dgamma(x-lower,shape=shape,scale=scale,log=log)
}

#' @importFrom stats rgamma
rgamma_l <- function(n,shape,scale,lower) {
  # Used with gamma prior
  lower + rgamma(n,shape=shape,scale=scale)
}

#' @importFrom stats dlnorm
dlnorm_l <- function(x,meanlog,sdlog,lower,log=FALSE) {
  # Used with lognormal prior
  dlnorm(x-lower,meanlog,sdlog,log=log)
}

#' @importFrom stats rlnorm
rlnorm_l <- function(n,meanlog,sdlog,lower) {
  # Used with lognormal prior
  lower + rlnorm(n,meanlog,sdlog)
}

dconstant <- function(x,constant,log=FALSE) {
  # Used with constant prior
  if (log) rep(0,length(constant)) else
    rep(1,length(constant))
}

rconstant <- function(n, constant, sd, lower) {
  # Used by DMC's constant prior; sd and lower are redundant arguments
  rep(constant, n)
}

### Samplers  ------------------------------------------------------------------

#' DE-MC Crossover and Migrate Samplers
#'
#' \code{crossoverR} and \code{migrateR} are two R prototypes of DE-MC
#' samplers, created by AH. \code{crossover} and \code{migrate} derive from
#' \code{crossoverR} and \code{migrateR} are optimised C++ DE-MC samplers.
#' \code{crossoverR} updates one chain at a time at the data level.
#' \code{crossover} updates all chains at once in a C++ for loop.
#' \code{migrateR} updates all chains together at the data level.
#'
#' @param k the index of a current processed chain. This must be an integer.
#' @param pars an integer vector storing parameter index. For example if
#' \code{pvec <- c(a = 1.15, v.f1 = 1.25, v.f2 = 1.85, z = 0.55, sz = 0.15,
#' sv = 0.32, t0 = 0.25)}, pars is \code{c(1:7)}.
#' @param use.theta a nchain x npar named matrix.
#' @param use.logprior a \code{nchain} length summed and logged prior vector.
#' Each element is the sum-log prior likelihood for a chain.
#' @param use.loglike a \code{nchain} length logged model likelihood. Each
#' element is logged model likelihood for a chain.
#' @param p.prior prior list created by prior.p.dmc
#' @param data data model instance
#' @param rp tuning parameter for DEMC. Default is 0.001.
#' @param gamma.mult turning parameter for DEMC. Default is 2.38
#' @param force use only for PDA method, a likelihood re-calculation interval
#' @return crossover returns a column named vector, representing the crossover
#' result for a chain. migrate returns a column named matrix, representing the
#' migration result for all chains.
#' @seealso prior.p.dmc
#' @examples
#' require(ggdmc)
#' m1 <- model.dmc(
#'   p.map     = list(a="1", v="1", z="1", d="1", sz="1", sv="1", t0="1",
#'                    st0="1"),
#'   constants = c(st0=0, d=0),
#'   match.map = list(M = list(s1="r1", s2="r2")),
#'   factors   = list(S = c("s1", "s2")),
#'   responses = c("r1", "r2"),
#'   type      = "rd")
#'
#' pvec1 <- c(a = 1, v = 1, z = .5, sz = .25, sv = .2,t0 = .15)
#' dat1  <- simulate(m1, pvec1, 512)
#' mdi1  <- BindDataModel(dat1, m1)
#'
#' prior1 <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' system.time(s0 <- run(samples.dmc(1e2, prior1, mdi1), report = 1e2,
#' p.migrate = .05))
#'
#' pars <- 1:s0$n.pars
#' th <- s0$theta[,, s0$start]
#' lp <- s0$summed_log_prior[s0$start,]
#' ll <- s0$log_likelihoods[s0$start,]
#' prior <- s0$p.prior
#' dat <- s0$data
#' rp <- s0$rp
#' force <- FALSE
#' gm <- 2.38
#'
#' res1 <- crossover(1, pars, th, lp, ll, prior, dat, rp, force, gm)
#' res2 <- migrate(th, lp, ll, p.prior, data, rp)
#'
#' print(res1)
#' ##    lp       ll        a        v        z       sz       sv       t0
#' ## -1.60 -9727.18     1.94     1.36     0.43     0.07     0.69     0.26
#'
#' print(res2)
#' ##          lp          ll     a    v    z   sz   sv   t0
#' ##       -0.02    -4765.14  2.33 2.04 0.56 0.43 1.45 0.20
#' ##       -1.58    -9744.38  1.94 1.36 0.43 0.07 0.69 0.26
#' ##        3.50   -12782.53  2.35 2.89 0.56 0.26 0.75 0.31
#' ##        1.18   -12134.66  1.00 2.66 0.58 0.27 1.39 0.33
#' ##        0.20    -8061.58  2.55 2.93 0.71 0.30 0.95 0.22
#' ##        0.20    -7996.31  2.54 2.93 0.71 0.29 0.95 0.22
#' ##        1.97   -12600.06  2.14 2.53 0.35 0.47 1.01 0.33
#' ##        3.30   -12587.68  1.71 2.67 0.50 0.39 0.69 0.33
#' ##        1.96    -9775.30  1.64 2.44 0.43 0.26 1.61 0.28
#' ##        3.30   -12583.47  1.71 2.67 0.50 0.39 0.69 0.33
#' ##        3.98   -11536.49  2.35 2.36 0.55 0.37 0.94 0.29
#' ##        2.98   -11872.35  2.62 2.96 0.51 0.35 0.75 0.29
#' ##        1.97    -8572.99  1.87 3.01 0.36 0.19 0.89 0.25
#' ##        3.91   -10199.54  1.62 2.10 0.55 0.29 1.02 0.29
#' ##        3.91   -10210.39  1.62 2.10 0.55 0.29 1.02 0.29
#' ##        3.61   -13393.09  1.74 2.59 0.58 0.29 0.78 0.34
#' ##        2.68    -6582.30  2.14 2.27 0.41 0.30 0.86 0.22
#' ##       -0.03    -4699.96  2.33 2.03 0.56 0.43 1.45 0.20
#'
#' ## New crossover
#' s0 <- ggdmc:::samples(100, p.prior, data, 4, verbose = TRUE)
#' useTheta <- s0$theta[,,1]
#' useLogPrior <- s0$summed_log_prior[1, ]
#' useLogLike <- s0$log_likelihoods[1, ]
#' pnames <- names(p.prior)
#'
#' dists <- numeric(length(p.prior))
#' p1 <- numeric(length(p.prior))
#' p2 <- numeric(length(p.prior))
#' lower <- numeric(length(p.prior))
#' upper <- numeric(length(p.prior))
#' islog <- numeric(length(p.prior))
#' for(i in 1:length(p.prior)) dists[i] <- attr(p.prior[[i]], "dist")
#' for(i in 1:length(p.prior)) p1[i] <- p.prior[[i]][1][[1]]
#' for(i in 1:length(p.prior)) p2[i] <- p.prior[[i]][2][[1]]
#' for(i in 1:length(p.prior)) lower[i] <- p.prior[[i]][3][[1]]
#' for(i in 1:length(p.prior)) upper[i] <- p.prior[[i]][4][[1]]
#' for(i in 1:length(p.prior)) islog[i] <- p.prior[[i]][5][[1]]
#'
#'
#' model    <- attr(data, "model")
#' ise      <- attr(data, "cell.empty")
#' allpar   <- attr(model, "all.par")
#' parnames <- attr(model, "par.names")
#' type     <- attr(model, "type")
#' dim1 <- dimnames(model)[[1]]
#' dim2 <- dimnames(model)[[2]]
#' dim3 <- dimnames(model)[[3]]
#' n1idx    <- attr(model, "n1.order")
#' isr1     <- ggdmc:::check_rd(type, model)
#' cellidx  <- ggdmc:::cellidxmat(data)
#' matchcell<- attr(model, "match.cell")
#' RT <- data$RT
#'
#'
#' tmp <- ggdmc:::crossover(useTheta, useLogPrior, useLogLike, pnames, dists,
#'                          p1, p2, lower, upper, islog,
#'                          allpar, parnames, model, type, dim1, dim2, dim3,
#'                          n1idx, ise, cellidx, RT,
#'                          matchcell, isr1, .001, 2.38, FALSE, 16384, .01)
#'
#' @export
#' @importFrom stats runif
crossoverR <- function(k, pars, use.theta, use.logprior, use.loglike, p.prior,
  data, rp = .001, gamma.mult=2.38, force=FALSE)
{
  ## step size
  if (is.na(gamma.mult)) gamma <- runif(1, 0.5, 1) else
    gamma <-  gamma.mult/sqrt(2*length(pars))

  ## pick two other chains
  index <- sample(c(1:dim(use.theta)[1])[-k], 2)

  # DE step
  theta <- use.theta[k,]
  names(theta) <- names(attr(attributes(data)$model, "p.vector"))

  theta[pars] <-
    use.theta[k,pars] +
    gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) +
    runif(1, -rp, rp)

  ## If resample flag is TRUE, get a new prior and likelihood and put them into
  ## the output storage. That is, no check if it passes Metropolis step and
  ## accept it straight.
  if (force) {
    use.logprior[k] <- summed.log.prior(theta, p.prior)
    use.loglike[k]  <- sum(log.likelihood (theta, data))
  }
  summed.use.post <-  use.loglike[k] + use.logprior[k]

  ## Get a pair of new prior and likelihood and put them into temporary storage
  log.prior <- summed.log.prior(theta, p.prior)
  loglike  <- sum(log.likelihood (theta, data))
  post <- log.prior + loglike
  if ( is.na(post) ) post <- -Inf
  rho <- exp(post - summed.use.post)
  if(is.na(rho)) rho <- -Inf

  # Metropolis step
  if ( runif(1) < rho ) {
    use.theta[k,]   <- theta
    use.logprior[k] <- log.prior
    use.loglike[k]  <- loglike
  }

  c(use.logprior[k], use.loglike[k], use.theta[k,])
}



### Utilities  ----------------------------------------------------------------

##' @export
GetPNames <- function(model) { return(names(attr(model, "p.vector"))) }

##' Check Data-Model Instance
##'
##' Return a model object extracted from either a data-model instance (DMI) or
##' a 'samples' object. Before that check through (1) in the case of a single
##' subject, whether DMI is a \code{data.frame}; (2) in the case of 'samples',
##' whether samples is a list of many subjects; (3) whether model is
##' successfully created; and (4) whether p.prior is suppled or we can extract
##' it from 'samples'
##'
##' @param data a data-model instance
##' @param p.prior a parameter prior list
##' @param theta1 a user-supplied theta cube
##' @param nchain number of MCMC chains
##' @export
CheckDMI <- function(data = NULL, p.prior = NULL, theta1 = NULL,
  nchain = NULL) {

  if (!is.null(data) && !is.data.frame(data)) stop("Data must be a data frame")
  if (is.null(data)) stop("No data-model instance") else model <- attr(data, "model")
  npar <- length(GetPNames(model))
  if (is.null(nchain)) nchain <- 3*npar

  if (is.null(model)) stop("Must specify a model")
  if (is.null(p.prior)) stop("Must specify a p.prior argument")
  if (!is.null(theta1) && !is.matrix(theta1) || (!all(dim(theta1)==c(nchain, npar))))
    stop("theta1 must be a nchain x npar matrix")
  return(model)
}

##' @export
CheckSamples <- function(samples = NULL, p.prior = NULL, theta1 = NULL) {

  if (is.null(samples)) stop("Must supply samples") else data <- samples$data
  if (!is.null(data) && !is.data.frame(data)) stop("Data must be a data frame")
  if (!is.null(samples) && is.null(samples$theta)) stop("Use StartHypersamples")
  if (is.null(data)) model <- attr(samples$data, "model") else model <- attr(data, "model")
  npar <- length(GetPNames(model))
  nchain <- samples$n.chain

  if (is.null(model)) stop("Must specify a model")
  if (is.null(p.prior) && is.null(samples)) stop("Must specify a p.prior argument")
  if (!is.null(theta1) && !is.matrix(theta1) || (!all(dim(theta1)==c(nchain, npar))))
    stop("theta1 must be a nchain x npar matrix")
  return(model)
}

#' @export
makeforce <- function(samples, force) {
  nsamp <- 1 + (samples$nmc - samples$start) * samples$thin
  if (!force) {
    out <- rep.int(0, nsamp)
  } else {
    out <- rep.int(0, nsamp)
    if (!is.numeric(force)) stop("force must be an integer.")
    for(i in 1:nsamp) { if (i %% force == 0) {out[i] <- 1} }
  }
  return(out)
}

### Sampling  ----------------------------------------------------------------

#' Initialize New Samples
#'
#' @param nmc numbers of MCMC iterations/samples
#' @param p.prior parameter priors
#' @param data a data-model instance from \code{BindDataModel}.
#' @param thin thinning length.
#' @param samples a posterior sample. If data is not supplied, the function
#' will look for \code{samples}.
#' @param theta1 user-supplied theta
#' @param rp DE-MCMC tuning parameter. Default is 0.001.
#' @param nchain numbers of Markov chains. Default is 3 times numbers of
#' parameter.
#' @examples
#' m1 <- BuildModel(
#'     p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'     constants = c(st0=0,d=0),
#'     match.map = list(M=list(s1="r1",s2="r2")),
#'     factors   = list(S=c("s1","s2"), F=c("f1", "f2")),
#'     responses = c("r1","r2"),
#'     type      = "rd")
#'
#' ## m1 is "dmc" class
#' class(m1)
#' ## [1] "dmc"
#'
#' pVec <- c(a=1, v.f1=1, v.f2=1.5, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat  <- simulate(m1, nsim=1e2, p.vector=pVec)
#' str(dat)
#' ## 'data.frame':	400 obs. of  4 variables:
#' ## $ S : Factor w/ 2 levels "s1","s2": 1 1 1 1 1 1 1 1 1 1 ...
#' ## $ F : Factor w/ 2 levels "f1","f2": 1 1 1 1 1 1 1 1 1 1 ...
#' ## $ R : Factor w/ 2 levels "r1","r2": 1 1 1 2 1 1 1 1 2 1 ...
#' ## $ RT: num  0.26 0.255 0.572 0.25 0.518 ...
#'
#' ## mdi1 is dmc as well as data.frame class
#' mdi1 <- BindDataModel(dat, m1)
#'
#' p.prior <- prior.p.dmc(
#'    dists = rep("tnorm", 7),
#'    p1    = c(a=2,  v.f1=2.5, v.f2=1.25, z=.5, sz=.3, sv=1,  t0=.3),
#'    p2    = c(a=.5, v.f1=.5,  v.f2=.35,  z=.1, sz=.1, sv=.3, t0=.05),
#'    lower = c(0,-5, -5, 0, 0, 0, 0),
#'    upper = c(5, 7,  7, 2, 2, 2, 2))
#'
#' ## Set up a new DMC sample with 200 iteration. The default thinning length
#' ## is 1
#' samples0 <- samples(nmc=50, p.prior=p.prior, data=mdi1)
#' samples0$nmc
#' ## [1] 50
#'
#' ## Run a fixed-effect model with 5% chance of using migration sampler
#' ## samples0 <- run(samples0, p.migrate=.05)
#'
#' ## Add 200 more iteration on to sample0
#' samples1 <- samples(nmc=50, p.prior=p.prior, samples=samples0, add=TRUE)
#' ## samples1$nmc
#' ## [1] 100
#' @export
StartNewsamples <- function(nmc, data = NULL, p.prior = NULL, thin = 1,
  theta1 = NULL, rp = .001, nchain = NULL) {

  model  <- CheckDMI(data, p.prior, theta1, nchain)
  pnames <- GetPNames(model)
  npar   <- length(pnames)
  if (is.null(nchain)) nchain <- 3*npar

  if (is.null(data)) stop("No Data-Model instance.")
  if (is.null(p.prior)) stop("No p.prior no new samples")
  if (length(p.prior) != npar) stop("p.vector and p.prior incompatialbe")

  ## Temporary measures
  if (is.null(attr(data, "n.pda"))) attr(data, "n.pda") <- 2^14
  if (is.null(attr(data, "bw")))    attr(data, "bw")    <- .01
  if (is.null(attr(data, "debug"))) attr(data, "debug") <- 0
  if (is.null(attr(data, "gpuid"))) attr(data, "gpuid") <- 0
  ncore <- 1
  debug <- 0

  out <- init_new(nmc, p.prior, data, rp, thin, nchain, ncore, debug)
  return(out)
}

##' @export
RestartSamples <- function(nmc, samples = NULL, p.prior = NULL, thin = NULL,
  theta1 = NULL, rp = .001, add = FALSE) {

  model  <- CheckSamples(samples, p.prior, theta1)
  pnames <- GetPNames(model)
  npar   <- length(pnames)
  if (is.null(samples)) stop("Use StartNewsamples")
  if (is.null(thin)) thin <- samples$thin

  if (is.null(attr(samples$data, "n.pda"))) attr(samples$data, "n.pda") <- 2^14
  if (is.null(attr(samples$data, "bw")))    attr(samples$data, "bw")    <- .01
  if (is.null(attr(samples$data, "debug"))) attr(samples$data, "debug") <- 0
  if (is.null(attr(samples$data, "gpuid"))) attr(samples$data, "gpuid") <- 0

  if (add) {
    out <- init_add(nmc, samples, rp, thin) # add onto existing one
  } else {
    out <- init_old(nmc, samples, rp, thin) # start afresh
  }
  return(out)
}

##' @export
CheckHyperDMI <- function(data = NULL, nchain = NULL) {
  if (is.null(data)) stop("No data")
  if (!is.list(data)) stop ("data must be a list")
  if (is.data.frame(data)) stop("data is a list with each if its elements is data.frame")
  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  if (is.null(nchain)) nchain <- length(pnames)
  return(nchain)
}


# hsam0 <- ggdmc::StartNewHypersamples(32, dmi, p.prior, pp.prior)

##' @export
StartNewHypersamples <- function(nmc, data = NULL, p.prior = NULL,
  pp.prior = NULL, thin = 1, rp = .001, nchain = NULL) {
  nchain <- CheckHyperDMI(data, nchain) ## If nchain=NULL, change it to default
  ## nsub   <- length(data)

  if (is.null(p.prior)) stop("No p.prior")   ## Check priors
  if (is.null(pp.prior)) stop("No pp.prior")
  if (!is.list(pp.prior)) stop("pp.prior must be a list")
  if (length(pp.prior[[1]]) < length(pp.prior[[2]]))
    stop("Location priors must have as many or more elements than scale priors")

  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  isppriorok <- pnames %in% names(p.prior)
  islpriorok <- pnames %in% names(pp.prior[[1]])
  isspriorok <- pnames %in% names(pp.prior[[2]])
  if (!all(isppriorok)) stop("p.prior incompatible with model1")
  if (!all(islpriorok)) stop("location prior incompatible with model1")
  if (!all(isppriorok)) stop("scale prior incompatible with model1")

  out <- init_newhier(nmc, data, p.prior, pp.prior, rp, thin, nchain)
  names(out) <- names(data)
  return(out)
}


#' Set up a DMC Sample with Multiple Participants
#'
#' \code{hsample} initialise a DMC object with each particpant as a list
#' element in a list. A participant is himself/herself a list element, carrying
#' the DMC setting, such as theta, summed_log_likelihood, etc.
#'
#' @param nmc number of Markov Chain Monte Carlo iteration
#' @param p.prior prior distribution setting
#' @param data a model data instance created by \code{data.model.dmc}
#' @param pp.prior prior distribution setting, hyper level
#' @param samples a DMC posterior sample
#' @param thin thinning length. Default 1
#' @param theta1 A user supplied initial theta cube
#' @param phi1 A user supplied initial phi cube
#' @param start.prior A user supplied (different) prior distribution setting
#' @param hstart.prior A user supplied (different) hyper-prior
#' distribution setting
#' @param add whether add new MCMC iteration on top of previous DMC sample
#' @param rp a DEMCMC tuning parameter
#' @param setting a list container to store all DMC setting
#' @param verbose whether to print debugging information
#' @examples
#' ## Set up a DDM Model, rate effect of factor F
#' m1 <- model.dmc(
#'   p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'   match.map = list(M=list(s1="r1",s2="r2")),
#'   factors=list(S=c("s1","s2"), F=c("f1","f2")),
#'   constants = c(st0=0,d=0),
#'   responses = c("r1","r2"),
#'   type = "rd")
#'
#' ## Population distribution
#' pop.mean  <- c(a=2,   v.f1=2.5, v.f2=1.5, z=0.5, sz=0.3, sv=1,  t0=0.3)
#' pop.scale <- c(a=0.5, v.f1=.5,  v.f2=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm",7),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0, 0, 0),
#'   upper = c(5, 7,  7, 2, 2, 2, 2))
#'
#' ##  Check population distributions
#' plot_priors(pop.prior)
#'
#' ## Simulate some data
#' dat <- h.simulate.dmc(m1, nsim=20, ns=4, p.prior=pop.prior)
#' mdi <- BindDataModel(dat, m1)
#' head(dat)
#' ##    S  F  R        RT
#' ## 1 s1 f1 r1 0.9227881
#' ## 2 s1 f1 r1 0.7878554
#' ## 3 s1 f1 r1 0.4814711
#' ## 4 s1 f1 r1 0.6864110
#' ## 5 s1 f1 r1 0.5068179
#' ## 6 s1 f1 r1 0.6356547
#'
#' ## Take a look at true parameters
#' ps <- round( attr(dat, "parameters"), 2)
#' ps
#' ##      a v.f1 v.f2    z   sz   sv   t0
#' ## 1 2.83 2.91 1.41 0.66 0.30 0.65 0.30
#' ## 2 2.37 2.42 2.24 0.48 0.28 1.14 0.31
#' ## 3 1.91 2.49 0.98 0.74 0.33 1.20 0.18
#' ## 4 2.14 2.67 2.34 0.65 0.31 1.74 0.27
#'
#' ## FIT FIXED EFFECTS
#' ## specify a broader prior than the true population distribution
#' p.prior <- prior.p.dmc(
#'   dists= rep("tnorm", length(pop.mean)),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' ## Set up multiple participant DMC sample
#' samples0 <- h.samples.dmc(nmc=100, p.prior=p.prior, data=mdi, thin=1)
#'
#' ## Add 400 more iterations and change thinning length to 2
#' samples1 <- h.samples.dmc(nmc=100, p.prior=p.prior, samples=samples0,
#' thin=2, add=TRUE)
#'
#' samples1[[1]]$nmc
#' ## [1] 200
#' @export
hsamples <- function(nmc, p.prior = NULL, data = NULL, thin = NA, samples = NULL,
                     restart = TRUE, add = FALSE, rp  =.001, n.chains = NULL)
{

   if (is.null(samples)) {
     out <- vector(mode = "list", length = length(data))
     if (is.na(thin)) thin <- 1
     if (is.null(data)) stop("data not found.")

     for(i in 1:length(data)) out[[i]] <- samples(nmc, p.prior, data[[i]],
       thin, samples, restart, add, rp, n.chains)

   } else {
      out <- vector(mode = "list", length = length(samples))
      if (is.na(thin)) thin <- samples[[1]]$thin
      for (i in 1:length(samples)) out[[i]] <- samples(nmc, p.prior, data,
          thin, samples[[i]], restart, add, rp, n.chains)
   }

  return(out)
}

##' Fit a Bayesian Model to a Single Participant
##'
##' Use either DE-MCMC or DGMC sampler to fit Bayesian model to a participant
##'
##' @param samples a initialized sample
##' @param report progress report intervel
##' @param ncore number of CPU cores
##' @param pm probability of migration
##' @param qm probability of mutation
##' @param ngroup number of independent groups
##' @param force PDA re-calculate interval
##' @param sampler a string indicating to use which sampler
##' @return Bayesian samples
##' @export
run_one <- function(samples, report, ncore, pm, qm, gammamult, ngroup, force,
  sampler) {

  force  <- makeforce(samples, force)
  pnames <- GetPNames(attr(samples$data, "model"))

  if (is.null(attr(samples$data, "n.pda"))) attr(samples$data, "n.pda") <- 2^14
  if (is.null(attr(samples$data, "bw")))    attr(samples$data, "bw")    <- .01
  if (is.null(attr(samples$data, "debug"))) attr(samples$data, "debug") <- 0
  if (is.null(attr(samples$data, "gpuid"))) attr(samples$data, "gpuid") <- 0

  if (sampler == "DGMC") {
    out <- run_dgmc(samples, force, report, pm, qm, gammamult, ncore, ngroup)
  } else if (sampler == "DE-MCMC") {
    out <- run_dmc(samples, force, report, pm, gammamult, ncore)
  } else {
    stop ("Sampler yet implemented")
  }

  dimnames(out$theta) <- list(NULL, pnames, NULL)
  return(out)
}

##' Fit a Bayesian Model to multiple Participants
##'
##' Use either DE-MCMC or DGMC sampler to fit independent Bayesian model to
##' many participants.
##'
##' @param samples a initialized samples list. Each element should contain
##' samples for a participant.
##' @param report progress report intervel
##' @param ncore number of CPU cores
##' @param pm probability of migration
##' @param qm probability of mutation
##' @param ngroup number of independent groups
##' @param force PDA re-calculate interval
##' @param sampler a string indicating to use which sampler
##' @return Bayesian samples
##' @export
##' @export
run_many <- function(samples, report, ncore, pm, qm, gammamult, ngroup,
  force, sampler) {

  force <- makeforce(samples[[1]], force)
  for(i in 1:length(samples)) {
    if (is.null(attr(samples[[i]]$data, "n.pda"))) attr(samples[[i]]$data, "n.pda") <- 2^14
    if (is.null(attr(samples[[i]]$data, "bw")))    attr(samples[[i]]$data, "bw")    <- .01
    if (is.null(attr(samples[[i]]$data, "debug"))) attr(samples[[i]]$data, "debug") <- 0
    if (is.null(attr(samples[[i]]$data, "gpuid"))) attr(samples[[i]]$data, "gpuid") <- 0
  }

  if (get_os() == "windows" & ncore > 1) {
    cl  <- parallel::makeCluster(ncore)
    if (sampler == "DGMC") {
      out <- parallel::parLapply(cl, samples, run_dgmc, force, report, pm, qm,
        gammamult, ncore, ngroup)
    } else if (sampler == "DE-MCMC") {

      out <- parallel::parLapply(cl, samples, run_dmc, force, report, pm,
        gammamult, ncore)
    } else {
      stop("Sampler unknown")
    }
    stopCluster(cl)

  } else if (ncore > 1) {
    if (sampler == "DGMC") {
      out <- parallel::mclapply(samples, run_dgmc, force, report, pm, qm,
        gammamult, ncore, ngroup, mc.cores=ncore)
    } else if (sampler == "DE-MCMC") {
      out <- parallel::mclapply(samples, run_dmc, force, report, pm,
        gammamult, ncore, mc.cores=ncore)
    } else {
      stop("Sampler unknown")
    }

  } else {
    if (sampler == "DGMC") {
      out <- lapply(samples, run_dgmc, force, report, pm, qm, gammamult, ncore,
        ngroup)
    } else if (sampler == "DE-MCMC") {
      out <- lapply(samples, run_dmc, force, report, pm, gammamult, ncore)
    } else {
      stop ("Sampler unknown")
    }
  }

  for(i in 1:length(out)) {
    dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
  }
  return(out)
}

#' Fit a fixed-effect or hierarchical Bayesian models
#'
#' This function fit a hierarchical or a fixed-effect model, using Bayeisan
#' sampling.  We use pMCMC, with a suite of DE-MCMC, DGMC, and simply,
#' crossover (i.e., DE-MC), mutation, or migration operators. Note that
#' the latter two operators essentially are random-walk Metroplolis, so they
#' will be very inefficient, if been applied alone, even with our fast C++
#' implementation.
#'
#' @param samples a sample list generated by calling DMC's samples.dmc.
#' @param report how many iterations to return a report
#' @param cores a switch for computing the prob density for each trial in
#' parallel. Turn it on by setting any number > 1.
#' @param p.migrate set it greater than 0 to use migration samplers. For example
#' p.migrate=0.05 will use migration in 5\% chance.
#' @param gamma.mult a DEMC tuning parameter, affecting the size of jump
#' @param farjump No funciton for compatibility reason
#' @param force Set force to FALSE for turning off force resampling. Set it
#' as an integer 1 to 10, forcing to resample a new parameter proposal every,
#' e.g., 1, 2, 3 step.
#' @return a DMC sample with class c("list", "dmc")
#' @examples
#' m1 <- model.dmc(
#'     p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'     constants = c(st0=0, d=0),
#'     match.map = list(M=list(s1="r1", s2="r2")),
#'     factors   = list(S=c("s1", "s2")),
#'     responses = c("r1", "r2"),
#'     type      = "rd")
#'
#' ## Use 6 prior truncated normal distributions
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## parameter vector. These are the trun values simulating data
#' p.vector <- c(a=1,v=1, z=.5, sz=.25, sv=.2,t0=.15)
#' dat1 <- simulate(m1, p.vector, 128)
#' dmi1 <- BindDataModel(dat1, m1)
#'
#' ## Use DMC's plot_cell_density to examine distributions
#' ## Accuracy around 70%
#' par(mfrow=c(1,2))
#' ggdmc:::plot.cell.density(data.cell=dmi1[dmi1$S=="s1", ], C="r1", xlim=c(0,2))
#' ggdmc:::plot.cell.density(data.cell=dmi1[dmi1$S=="s2", ], C="r2", xlim=c(0,2))
#' par(mfrow=c(1,1))
#'
#' ## ---------------------------
#' ## Profiles all 6 parameters
#' par(mfrow=c(2,3));
#' ggdmc:::profile.dmc("a",  .1,  2, p.vector, dmi1)
#' ggdmc:::profile.dmc("v",  .1,  2, p.vector, dmi1)
#' ggdmc:::profile.dmc("z",  .2, .8, p.vector, dmi1)
#' ggdmc:::profile.dmc("sz", .1, .9, p.vector, dmi1)
#' ggdmc:::profile.dmc("sv", .1,  2, p.vector, dmi1)
#' ggdmc:::profile.dmc("t0", .1,  .2, p.vector, dmi1)
#' par(mfrow=c(1,1));
#'
#' ## Initialse a DMC sample
#' ## nthin == 1 (default)
#' ## niter == 100
#' ## prior distributions as listed in p.prior
#' ## data  == model data instance 1
#' ## do not use migrate sampler (default p.migrate=0)
#' m1 <- model.dmc(
#'   p.map     = list(a="1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
#'   constants = c(st0=0, d=0),
#'   match.map = list(M = list(s1="r1", s2="r2")),
#'   factors   = list(S = c("s1", "s2")),
#'   responses = c("r1", "r2"),
#'   type      = "rd")
#'
#' pvec1 <- c(a = 1, v = 1, z = .5, sz = .25, sv = .2,t0 = .15)
#' dat1  <- ggdmc:::simulate.model(m1, pvec1, 512)
#' dmi1  <- BindDataModel(dat1, m1)
#' prior1 <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## Because this will take around 20 to 40 s, set it NOT RUN.
#' ##   user  system elapsed
#' ## 17.552   0.004  17.551
#' \dontrun{
#' s0 <- ggdmc:::samples(512, prior1, dmi1)
#' s0 <- ggdmc:::run(s0, 32, 1, .05, 2.38, FALSE)
#' s1 <- ggdmc:::samples(512, prior1, samples = s0)
#' s1 <- ggdmc:::run(s1, 32, 1)
#'
#' ggdmc:::plot.dmc(s0, pll.chain = TRUE)
#? ggdmc:::plot.dmc(s1, pll.chain = TRUE)
#' ggdmc:::plot.dmc(s1, density   = TRUE)
#' ggdmc:::check.recovery.dmc(s1, pvec1)
#' }
#' @export
run <- function(samples, report = 1e2, ncore = 1, pm = 0, qm = 0,
  gamma.mult = 2.38, ngroup = 6, force = FALSE, sampler = "DE-MCMC",
  debug = FALSE) {

  hyper <- attr(samples, "hyper")

  if (!is.null(hyper)) {   ## hierarchical model
    if (is.null(attr(samples[[1]]$data, "bw"))) stop("No GPU attributes")
    if (is.null(attr(samples[[1]]$data, "gpuid"))) stop("No GPU attributes")

    out <- run_hyper_dmc(samples, report, pm, gamma.mult, ncore, debug)
  } else if (any(names(samples) == "theta")) { ## One subject
    out <- run_one(samples, report, ncore, pm, qm, gamma.mult, ngroup, force,
      sampler)
  } else {  ## Multiple independent subjects
    out <- run_many(samples, report, ncore, pm, qm, gamma.mult, ngroup, force,
      sampler)
  }

  cat("\n")
  return(out)
}


#' @rdname run
#' @export
CheckConverged <- function(samples) {
  stuck <- ggdmc::StuckTests(samples, verbose = FALSE, cut = 10)
  flat  <- ggdmc::FlatTests(samples, p1 = 1/3, p2 = 1/3,
    cut.location = 0.25, cut.scale = Inf, verbose = FALSE)
  mix  <- ggdmc::MixTests(samples, cut = 1.05, verbose = FALSE)
  size <- ggdmc::LengthTests(samples, minN = 512, nfun = "mean", FALSE)
  isstuck <- TRUE
  if (stuck == 0) isstuck <- FALSE

  out <- c(isstuck, flat, mix, size)
  names(out) <- c("Stuck", "Flat", "Mix", "ES")
  return(out)
}

#' @rdname run
#' @export
RunTillConverged <- function(nmc, p.prior, data, thin, samples = NULL, report = 128,
                             ncore = 1, p.migrate = 0, gamma.mult = 2.38,
                             force = FALSE, effective.Size = 512, times = 8)
{
  if (is.null(data) & is.null(samples)) stop("Neither data nor samples was found")
  if (is.null(thin)) thin <- 2
  if (!is.null(data)) {
    sam <- ggdmc:::run(ggdmc:::samples(nmc, p.prior, data = data, thin = thin),
                       report, ncore, p.migrate = .05,
                       gamma.mult = gamma.mult, force = force)
  }

  ispassed <- !any(CheckConverged(sam))
  counter <- 1

  if (!ispassed) {
    thin <- 2^(counter+1)
    repeat {
      sam <- ggdmc:::run(ggdmc:::samples(nmc, p.prior, samples = sam,
                                         thin = thin), report, ncore)
      ispassed <- !any(CheckConverged(sam))
      counter  <- counter + 1
      if (ispassed | counter > times) break
    }
  }

  if (!ispassed) {
    cat("Have repeated ", times, "times, but yet converged.\n")
    cat("Chain(s): ", ggdmc:::pickStuck(sam), "was/were removed")
    sam <- ggdmc::unstick(sam, ggdmc:::pickStuck(sam))
  }
  return(sam)
}


trial_log_likes <- function(samples,thin_pointwise=1,
  chain_dot=TRUE,subject_dot=FALSE)
  # Get pointwise log-likelihoods
{

  n.trials <- dim(samples$data)[1]

  nmc_thin <- seq(thin_pointwise,samples$nmc,by=thin_pointwise)
  trial_log_likes  <-
    array(-Inf,c(length(nmc_thin),samples$n.chains, dim(samples$data)[1]))
  dimnames(trial_log_likes) <- list(nmc_thin,NULL,NULL)
  if (chain_dot) cat("Processing chains: ")
  for (j in 1:samples$n.chains) {
    for (i in nmc_thin)
      trial_log_likes[as.character(i),j,]  <-
        log.likelihood(samples$theta[j,,i],samples$data)
    if (chain_dot) cat(".")
  }
  if (chain_dot) cat("\n")
  if (subject_dot) cat(".")
  trial_log_likes
}

group_trial_log_likes <- function(samples,thin_pointwise=1,max_size=1e8)
  # extracts trial_log_likes from a list of subject fits and concatanates
{

  cat("Processing subjects: ")
  tll <- trial_log_likes(samples[[1]],thin_pointwise=thin_pointwise,
    chain_dot=FALSE,subject_dot=TRUE)
  sdim <- dim(tll)
  size <- sum(length(samples)*prod(sdim))
  if (size>max_size) stop(paste("Output size",size,
    "too large,adjust max_size (",max_size,")"))
  tlls <- lapply(samples[-1],trial_log_likes,thin_pointwise=thin_pointwise,
    chain_dot=FALSE,subject_dot=TRUE)
  tlls[[length(samples)]] <- tll
  sdims <- cbind(matrix(unlist(lapply(tlls,dim)),nrow=3))
  if ( !all(sdims[1,1]==sdims[1,-1]) || !all(sdims[2,1]==sdims[2,-1]) )
    stop("Subjects must have the same number of interations and chains")
  out <- array(dim=c(dim(tlls[[1]])[-3],sum(sdims[3,])))
  start <- 1; end <- sdims[3,1]
  for (i in 1:length(samples)) {
    out[,,start:end] <- tlls[[i]]
    if (i<length(samples)) {
      start <- end+1
      end <- start - 1 + sdims[3,i+1]
    }
  }
  out
}


group_subject_log_likes <- function(samples)
  # extracts summed_log_likes from a list of subject fits and concatanates
{

  tlls <- lapply(samples,function(x){x$log_likelihoods})
  sdims <- matrix(unlist(lapply(tlls,dim)),nrow=2)
  if ( !all(sdims[1,1]==sdims[1,-1]) || !all(sdims[2,1]==sdims[2,-1]) )
    stop("Subjects must have the same number of interations and chains")
  out <- array(dim=c(dim(tlls[[1]]),length(samples)))
  for (i in 1:length(samples))
    out[,,i] <- tlls[[i]]
  out
}


### Posterior predictives ----

#' @export
get.dqp <- function(sim,facs,probs,n.post=NA,ns=NA,bw="nrd0") {

  quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
    out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
    names(out) <- probs*100
    if ( all(is.na(out)) || length(x)==1) NULL else out
  }

  qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)

  # get probabilities
  n <- tapply(sim$RT,sim[,c(facs,"R")],length)
  n[is.na(n)] <- 0 # In case some cells are empty
  nok <- tapply(sim$RT,sim[,c(facs,"R")],function(x){sum(!is.na(x))})
  nok[is.na(nok)] <- 0 # In case some cells are empty
  if (is.null(facs)) np <- sum(n) else
    np <- rep(apply(n,1:length(facs),sum),times=length(levels(sim$R)))
  p <- nok/np

  # For a simulation get probability replicates
  if ( !is.na(n.post) && (n.post>1) ) {
    repfac <- rep(1:n.post,each=sum(ns))
    ps <- tapply(sim$RT,cbind(sim[,c(facs,"R"),drop=FALSE],rep=repfac),
      function(x){sum(!is.na(x))})
    ps[is.na(ps)] <- 0 # In case some cells are empty
    ps <- n.post*ps/np
  } else ps=NULL

  # cell names
  cell.names <- dimnames(qs)[[1]]
  if (!is.null(facs)) {
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
  }
  cell.names <- as.vector(cell.names)
  # Get density and make defective
  dens <- tapply(sim$RT,sim[,c(facs,"R"),drop=FALSE],function(x){
    if (all(is.na(x))) NULL else {
      x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
      if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
    }
  })
  for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
  {
    if (!is.null(qs[i][[1]])) {
      names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
      attr(qs[i][[1]],"cell.name") <- cell.names[i]
    }
    if (!is.null(dens[i][[1]]) ) {
      dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
      attr(dens[i][[1]],"cell.name") <- cell.names[i]
    }
  }
  dnd <- dimnames(dens)
  dnq <- dimnames(qs)
  dens <- apply(dens,1:length(facs),function(x){x})
  qs <- apply(qs,1:length(facs),function(x){x})
  if ( is.null(dim(dens)) ) {
    dens <- array(dens,dim=c(length(dens)))
    dimnames(dens) <- dnd[-length(dnd)]
    qs <- array(qs,dim=c(length(qs)))
    dimnames(qs) <- dnq[-length(dnq)]
  }
  list(pdf=dens,cdf=qs,n=n,p=p,ps=ps)
}


#' @export
#' @importFrom stats quantile
post.predict.dmc <-function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
  bw="nrd0",report=10,save.simulation=FALSE,factors=NA,
  save.simulation.as.attribute=FALSE,ignore.R2=FALSE,censor=c(NA,NA),
  gglist=FALSE, probs.gglist=c(0.1, 0.5, 0.9), CI.gglist=c(0.025, 0.975))
  # make list of posterior preditive density, quantiles and response p(robability)
  # NB: quantiles only calcualted for 2 or more RTs
{

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  cvs <- samples$data[,attr(model,"cvs")]
  attr(cvs,"row.facs") <- apply(apply(
    samples$data[,facs,drop=FALSE],2,as.character),1,paste,collapse=".")
  if ( ignore.R2 & any(names(samples$data)=="R2") )
    samples$data <- samples$data[,names(samples$data)[names(samples$data)!="R2"]]
  if (!is.null(factors) ) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs))
      stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  }
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT",names(cvs),"R2")]
  } else {
    # Assumes last two are SSD and RT! FIX ME. EG WONT WORK IF THERE ARE CVS
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- ggdmc:::simulate.model(model, posts[i,], n=ns)
    if (ignore.R2) tmp <- tmp[,names(tmp)[names(tmp)!="R2"]]
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  reps <- rep(1:n.post,each=dim(samples$data)[1])

  if (!is.na(censor[1])) fast <- sim[,"RT"] < censor[1] else fast <- rep(FALSE,dim(sim)[1])
  if (!is.na(censor[2])) slow <- sim[,"RT"] > censor[2] else slow <- rep(FALSE,dim(sim)[1])
  ok <- !fast & !slow
  sim <- sim[ok,]
  reps <- reps[ok]

  if ( save.simulation ) {
    sim <- cbind(reps,sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,facs=factors,probs,n.post,ns=ns,bw=bw)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    out <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[reps==i,]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(out,"dpqs") <- dpqs
    if (save.simulation.as.attribute)
      attr(out,"sim") <- cbind(reps,sim)
    if (gglist) attr(out, "gglist") <-
      get.fitgglist.dmc(sim=cbind(reps,sim),data=samples$data,factors=factors, noR=FALSE,
        quantiles.to.get= probs.gglist, CI = CI.gglist)
    out
  }
}



ppp.dmc <- function(samples,fun=function(x){mean(x$RT)},n.post=500,
  plot.density=TRUE,main="",bw="nrd0",report=10)
  # posterior predictive (pp) p value for function fun of data (p(observed)>pp)
{

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs])
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  posts <- thetas[sample(c(1:dim(thetas)[1]),n.post,replace=F),]
  sim <- vector(mode="list",length=n.post)
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in 1:n.post) {
    sim[[i]] <- ggdmc:::simulate.model(model, posts[i,], ns)
    if ( (i %% report) == 0) cat(".")
  }
  cat("\n")
  pp <- unlist(lapply(sim,fun))
  obs <- fun(samples$data)
  ppp <- mean(obs>pp)
  if (plot.density) {
    plot(density(pp,bw=bw),main=main)
    abline(v=obs,lty=2)
  }
  ppp
}

#' @export
run.unstuck.dmc <- function(samples,nmc=NA,report=10,cores=1,
  cut=10,nbad=0,max.try=100,p.migrate=0,
  gamma.mult=2.38,verbose=FALSE,
  end.no.migrate=FALSE)
  # Repeats sampling until <= nbad stuck chains as defined by cut or max.try
  # If samples has no content fills it in then repeatedly gets new sets of nmc
  # samples (nmc can be specified or taken from samples). If end.no.migrate
  # runs one final time with migration off.
{
  if (is.null(samples$theta))
    stop("For multiple subjects use h.run.unstuck.dmc")
  n.chain <- dim(samples$theta)[1]
  if (any(is.na(samples$theta[,,2])))
    samples <- run.dmc(samples=samples,report=report,cores=cores,
      gamma.mult=gamma.mult,p.migrate=p.migrate)
  if ( is.na(nmc) ) nmc <- samples$nmc
  try.num <- 1
  repeat {
    cat(paste("\nTry",try.num,"\n"))
    if ( length(pickStuck(samples, verbose = verbose)) <= nbad ) break
    samples <- run.dmc(samples.dmc(nmc=nmc,samples=samples),
      report=report,cores=cores,p.migrate=p.migrate,
      gamma.mult=gamma.mult)
    if (try.num >= max.try) break
    try.num <- try.num + 1
  }
  if (end.no.migrate) samples <- run.dmc(samples.dmc(nmc=nmc,samples=samples),
    report=report,cores=cores,gamma.mult=gamma.mult)
  samples
}


#' @export
run.converge.dmc <- function(samples,nmc,report=10,cores=1,gamma.mult=2.38,
  cut=1.1,max.try=100,minN=NA,meanN=NA,
  transform=TRUE,autoburnin=FALSE,split=TRUE,verbose=FALSE)
  # Adds samples repeatedly, throws away intial samples if gelman.diag better
  # and repeats until gelman.diag Multivariate psrf < cut and once that is
  # fulfilled will continue if necessary to get either min or mean effectiveSize
{
  if (!is.na(minN) & !is.na(meanN)) {
    warning("Both minN and meanN specified, using minN")
    meanN <- NA
  }
  if (!is.na(minN)) nfun <- "min"
  if (!is.na(meanN)) {
    nfun <- "mean"
    minN <- meanN
  }
  if ( is.null(samples$theta) )
    stop("For multiple subjects use h.run.converge.dmc")
  if ( any(is.na(samples$theta[,,2])) ) {
    samples <- run.dmc(samples=samples,
      report=report,cores=cores,gamma.mult=gamma.mult)
    gd <- gelman.diag.mpsrf(theta.as.mcmc.list(samples,split=split),
      autoburnin=autoburnin,transform=transform)
    if (verbose) cat(paste("MPSRF: ",gd,"\n"))
  } else gd <- Inf
  if ( !is.na(minN) & (gd <= cut) )
    okN <- do.call(nfun,list(effectiveSize.dmc(samples))) > minN else
      okN <- TRUE
    if ( (gd > cut) | !okN )
    {  # Do more sampling
      if ( is.na(nmc) ) nmc <- samples$nmc
      try.num <- 0
      effectiveN <- NA
      repeat {
        samples <- run.dmc(samples.dmc(samples=samples,add=TRUE,nmc=nmc),
          report=report,cores=cores,gamma.mult=gamma.mult)
        gd <- gelman.diag.mpsrf(theta.as.mcmc.list(samples,split=split),
          autoburnin=autoburnin,transform=transform)
        if ( try.num>0 ) {
          shorter <- samples.dmc(samples=samples,remove=1:nmc,nmc=0,add=TRUE)
          gd.short <- gelman.diag.mpsrf(theta.as.mcmc.list(shorter,split=split),
            autoburnin=autoburnin,transform=transform)
          if (gd.short < gd) {
            samples <- shorter
            gd <- gd.short
            if (verbose) cat(paste("Discarding initial",nmc,"samples.\n"))
          }
        }
        if ( try.num >= max.try ) break
        try.num <- try.num + 1
        if ( gd <= cut ) {
          if ( is.na(minN) ) break else {
            effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples)))
            if (effectiveN > minN) break
          }
        }
        if (verbose) print(paste("N =",dim(samples$theta)[3],
          "Effective N =",effectiveN,
          "Multivariate psrf achieved =",format(gd,digits=3)))
      }
    }
    if (verbose) {
      print(paste("Final multivariate psrf =",gd))
      cat("Effective sample size\n")
      print(effectiveSize.dmc(samples))
    }
    samples
}


#' @rdname run
#' @export
StuckTests <- function(samples, verbose=FALSE, cut = 10) {
  if (verbose) cat("Stuck chains check\n")
  stucks <- PickStuck(samples, cut = cut, verbose = verbose)
  fail <- length( stucks != 0)
  if (verbose) {
    if (!fail) cat(": OK\n") else
      cat(paste(":",length(stucks),"\n"))
  }
  fail
}

#' @rdname run
#' @export
FlatTests <- function(samples,p1=1/3,p2=1/3,cut.location=0.25,cut.scale=Inf,
  verbose=FALSE) {

  gmcmc <- function(samples) mcmc(matrix(aperm(samples$theta,c(1,3,2)),
    ncol=dim(samples$theta)[2],dimnames=list(NULL,dimnames(samples$theta)[[2]])))

  mat <- gmcmc(samples)
  xlen <- round(dim(mat)[1] * p1)
  ylen <- round(dim(mat)[1] * p2)
  # change in mean relative to robst SD
  m.zs <- apply(mat,2,function(x){
    m1 <- median(x[1:xlen])
    m2 <- median(x[(length(x)-ylen):length(x)])
    abs(m1-m2)/IQR(x)
  })
  names(m.zs) <- paste("m",names(m.zs),sep="_")
  fail <- any(m.zs>cut.location)
  if (!fail) out <- "" else
    out <- paste(names(m.zs)[m.zs==max(m.zs)],"=",round(max(m.zs),2))
  if ( is.finite(cut.scale) ) {
    # Change ini IQR reltiave to overall IQR
    s.zs <- apply(mat,2,function(x){
      m1 <- IQR(x[1:xlen])
      m2 <- IQR(x[(length(x)-ylen):length(x)])
      abs(m1-m2)/IQR(x)
    })
    names(s.zs) <- paste("s",names(s.zs),sep="_")
    if (out != "") out <- paste(out,", ",sep="")
    if (any(s.zs>cut.scale)) out <-
      paste(out,names(s.zs)[s.zs==max(s.zs)],"=",round(max(s.zs),2))
    fail <- fail | any(s.zs>cut.scale)
  }
  if (verbose) {
    cat("Flat check\n")
    print(round(m.zs,2))
    if ( is.finite(cut.scale) )
      print(round(s.zs,2)) else
        if (!fail) cat(": OK\n") else
          cat(paste(":",out,"\n"))
  }
  fail
}

#' @rdname run
#' @export
MixTests <- function(samples, verbose=FALSE, cut = 1.01, split=TRUE) {
  tmp <- gelman.diag.dmc(samples,split=split)
  gds <- c(tmp$mpsrf,tmp$psrf[,1])
  fail <- max(gds) > cut
  if (verbose) {
    cat("Mixing check\n")
    print(round(gds,2))
    if (!fail) cat(": OK\n") else {
      nam <- names(gds)[gds==max(gds)]
      cat(paste(":",nam,"=",round(max(gds),2),"\n"))
    }
  }
  fail
}

#' @rdname run
#' @export
LengthTests <- function(samples, minN, nfun, verbose=FALSE) {
  n <- do.call(nfun,list(effectiveSize.dmc(samples)))
  fail <- n < minN
  if (verbose) {
    cat("Length check")
    if (!fail) cat(": OK\n") else cat(paste(":",n,"\n"))
  }
  fail
}

#' @rdname run
#' @export
GetPLL <- function(fit, pars, start, end) {
  if (missing(pars)) pars <- fit@sim$pars_oi
  if (missing(start)) start <- fit@sim$warmup + 1
  if (missing(end)) end <- fit@sim$iter
  nchain <- fit@sim$chains

  pll <- matrix(numeric(nchain*(end - start + 1)), (end - start + 1))
  for (i in 1:nchain) {
    chaini  <- fit@sim$samples[[i]]$lp__
    pll[,i] <- chaini[start:end]
  }
  return(pll)
}

# StuckTests.stanfit <- function(fit, pars, start, end, cut = 10) {
#   if (missing(pars)) pars <- fit@sim$pars_oi
#   if (missing(start)) start <- fit@sim$warmup + 1
#   if (missing(end)) end <- fit@sim$iter
#   nchain <- fit@sim$chains
#
#   pll <- GetPLL(fit, pars, start, end)
#   mean.ll <- colMeans(pll)
#   names(mean.ll) <- 1:length(mean.ll)
#
#   dev <- -(sort(mean.ll) - median(mean.ll))
#   bad <- as.numeric(names(dev)[dev > 10])
#   return(bad)
# }

# FlatTests.stanfit <- function(fit, p1=1/3, p2=1/3,
#   cut.location=0.25, cut.scale=Inf, verbose=FALSE) {
#
#   gmcmc <- function(fit) {
#     nmc <- fit@sim$iter
#     nchain <- fit@sim$chains
#     npar <- length(fit@sim$fnames_oi) - 2
#
#     a <- array(0, dim=c(nmc, npar, nchain))
#     pnames <- matrix(numeric(nchain * npar), nchain)
#
#     for (i in 1:nchain) {
#       chaini <- fit@sim$samples[[i]]
#       for (j in 1:npar) {
#         a[,j,i] <- chaini[[j]]
#         pnames[i,j] <- names(chaini[j])
#       }
#     }
#
#     tran_theta <- aperm(a, c(3, 1, 2))
#     pnames <- pnames[1,]
#     tran.mat <- matrix(tran_theta, ncol = npar, dimnames = list(NULL, pnames))
#     return(coda::mcmc(tran.mat))
#   }
#   mat <- gmcmc(fit)
#
#   xlen <- round(dim(mat)[1] * p1)
#   ylen <- round(dim(mat)[1] * p2)
#   # change in mean relative to robst SD
#   m.zs <- apply(mat, 2, function(x){
#     m1 <- median(x[1:xlen])
#     m2 <- median(x[(length(x)-ylen):length(x)])
#     abs(m1-m2)/IQR(x)
#   })
#
#   names(m.zs) <- paste("m", names(m.zs),sep="_")
#   fail <- any(m.zs > cut.location)
#   if (!fail) out <- "" else
#     out <- paste(names(m.zs)[m.zs==max(m.zs)], "=", round(max(m.zs),2))
#   if ( is.finite(cut.scale) ) {
#     # Change ini IQR reltiave to overall IQR
#     s.zs <- apply(mat,2,function(x){
#       m1 <- IQR(x[1:xlen])
#       m2 <- IQR(x[(length(x)-ylen):length(x)])
#       abs(m1-m2)/IQR(x)
#     })
#     names(s.zs) <- paste("s",names(s.zs),sep="_")
#     if (out != "") out <- paste(out,", ",sep="")
#     if (any(s.zs>cut.scale)) out <-
#       paste(out,names(s.zs)[s.zs==max(s.zs)],"=",round(max(s.zs),2))
#     fail <- fail | any(s.zs>cut.scale)
#   }
#
#   if (verbose) {
#     cat("Flat check\n")
#     print(round(m.zs,2))
#     if ( is.finite(cut.scale) )
#       print(round(s.zs,2)) else
#         if (!fail) cat(": OK\n") else
#           cat(paste(":",out,"\n"))
#   }
#   fail
# }

# MixTests.stanfit <- function(fit, verbose=FALSE, cut = 1.05, split=TRUE) {
#   fit_summary <- rstan::summary(fit)
#   npar <- length(fit@sim$fnames_oi) - 2
#   gds <- fit_summary$summary[, "Rhat"][1:npar]
#   fail <- max(gds) > cut
#
#   if (verbose) {
#     cat("Mixing check\n")
#     print(round(gds,2))
#     if (!fail) cat(": OK\n") else {
#       nam <- names(gds)[gds==max(gds)]
#       cat(paste(":",nam,"=",round(max(gds),2),"\n"))
#     }
#   }
#   fail
# }

# LengthTests.stanfit <- function(fit, minN = 512, nfun = "mean",
#   verbose = FALSE) {
#   fit_summary <- rstan::summary(fit)
#   npar <- length(fit@sim$fnames_oi) - 2
#   es <- fit_summary$summary[,"n_eff"][1:npar]
#   n <- do.call(nfun, list(es))
#
#   fail <- n < minN
#   if (verbose) {
#     cat("Length check")
#     if (!fail) cat(": OK\n") else cat(paste(":", round(n, 2), "\n"))
#   }
#   fail
# }

