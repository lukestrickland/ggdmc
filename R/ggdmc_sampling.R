### Prior Distributions -------------------------------------------------------
##' @importFrom stats dbeta
dbeta_lu <- function(x, shape1, shape2, lower, upper, log = FALSE) {
  # Used with beta prior
  if (log) {dbeta((x-lower)/(upper-lower),shape1,shape2,log=log)}
  else {dbeta((x-lower)/(upper-lower),shape1,shape2,log=log)/(upper-lower)}
}

##' @importFrom stats rbeta
rbeta_lu <- function(n, shape1, shape2, lower, upper) {
  # Used with beta prior
  lower + rbeta(n,shape1,shape2)*(upper-lower)
}

##' @importFrom stats dgamma
dgamma_l <- function(x, shape, scale, lower, log = FALSE) {
  # Used with gamma prior
  dgamma(x-lower,shape=shape,scale=scale,log=log)
}

##' @importFrom stats rgamma
rgamma_l <- function(n, shape, scale, lower) {
  # Used with gamma prior
  lower + rgamma(n,shape=shape,scale=scale)
}

##' @importFrom stats dlnorm
dlnorm_l <- function(x, meanlog, sdlog, lower, log = FALSE) {
  # Used with lognormal prior
  dlnorm(x-lower,meanlog,sdlog,log=log)
}

##' @importFrom stats rlnorm
rlnorm_l <- function(n, meanlog, sdlog, lower) {
  # Used with lognormal prior
  lower + rlnorm(n,meanlog,sdlog)
}

dconstant <- function(x, constant, log = FALSE) {
  # Used with constant prior
  if (log) rep(0, length(constant)) else
    rep(1, length(constant))
}

rconstant <- function(n, constant, sd, lower) {
  # Used by DMC's constant prior; sd and lower are redundant arguments
  rep(constant, n)
}

### Utilities  ----------------------------------------------------------------
##' @export
GetPNames <- function(model) { return(names(attr(model, "p.vector"))) }

##' Check Data Model Instance
##'
##' Return a model object extracted either from a data model instance or
##' a 'samples' object. The function checks also (1) in the case of a single
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

  if (is.null(data)) {
    stop("No data-model instance")
  } else {
    model <- attr(data, "model")
  }

  npar <- length(GetPNames(model))
  if (is.null(nchain)) nchain <- 3*npar
  if (is.null(model)) stop("Must specify a model")
  if (is.null(p.prior)) stop("Must specify a p.prior argument")
  if (!is.null(theta1) && !is.matrix(theta1) ||
      (!all(dim(theta1)==c(nchain, npar))))
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
MakeForce <- function(samples, force) {
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
##' @export
CheckHyperDMI <- function(data = NULL, nchain = NULL) {
  if (is.null(data)) stop("No data")
  if (!is.list(data)) stop ("data must be a list")
  if (is.data.frame(data)) stop("data is a list with each if its elements is data.frame")
  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  if (is.null(nchain)) {
    nchain <- 3*length(pnames)
    message("nchain is not supplied. Use default ", nchain, " chains")
    message("If you are using DGMC, this nchain may be inappropriate.")
  }
  return(nchain)
}

##' Initialize New Samples
##'
##' These functions use prior distributions, either from \code{p.prior} or joinly
##' from \code{p.prior} and \code{pp.prior} in the case of hierarchical
##' models to generate over-dispersed initial parameter values.
##'
##'
##' @param nmc numbers of Monte Carlo samples / iterations.
##' @param data a data model instance from \code{BindDataModel}.
##' @param p.prior parameter prior distributions from \code{BuildPrior}.
##' @param pp.prior hyper parameter prior distributions from \code{BuildPrior}.
##' This must be a set of location and scale hyper prior distributions.
##' @param thin thinning length.
##' @param nchain numbers of Markov chains. Default is 3 times the numbers of
##' model parameters.
##' @param rp DE-MCMC tuning parameter to generate random noise either from
##' uniform or Gaussian distribution.
##' @param samples a collection fo posterior samples.
##' @param add add more samples onto an existing samples
#' @examples
#' m1 <- BuildModel(
#'     p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'     constants = c(st0=0, d=0),
#'     match.map = list(M = list(s1="r1",s2="r2")),
#'     factors   = list(S = c("s1","s2"), F = c("f1", "f2")),
#'     responses = c("r1","r2"),
#'     type      = "rd")
#'
#' ## m1 is "model" class
#' class(m1)
#' ## [1] "model"
#'
#' pVec <- c(a=1, v.f1=1, v.f2=1.5, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat  <- simulate(m1, nsim=1e2, ps = pVec)
#' str(dat)
#' ## 'data.frame':	400 obs. of  4 variables:
#' ## $ S : Factor w/ 2 levels "s1","s2": 1 1 1 1 1 1 1 1 1 1 ...
#' ## $ F : Factor w/ 2 levels "f1","f2": 1 1 1 1 1 1 1 1 1 1 ...
#' ## $ R : Factor w/ 2 levels "r1","r2": 1 1 1 2 1 1 1 1 2 1 ...
#' ## $ RT: num  0.26 0.255 0.572 0.25 0.518 ...
#'
#' dmi1 <- BindDataModel(dat, m1)
#' npar <- length(GetPNames(m1))
#'
#' p.prior <- BuildPrior(
#'    dists = rep("tnorm", npar),
#'    p1    = c(a=2,  v.f1=2.5, v.f2=1.25, z=.5, sz=.3, sv=1,  t0=.3),
#'    p2    = c(a=.5, v.f1=.5,  v.f2=.35,  z=.1, sz=.1, sv=.3, t0=.05),
#'    lower = c(0,-5, -5, 0, 0, 0, 0),
#'    upper = c(5, 7,  7, 2, 2, 2, 2))
#'
#' ## Set up a new DMC sample with 16 iteration. The default thin is 1
#' samples0 <- StartNewsamples(nmc = 16, data=dmi1, p.prior=p.prior)
#' samples0$nmc
#' ## [1] 16
#'
#' samples1 <- RestartSamples(16, samples0, p.prior)
#' samples1$nmc
#' ## [1] 16
#' samples2 <- RestartSamples(16, samples0, p.prior, add = TRUE)
#' samples2$nmc
#' ## [1] 32
#'
#' #######################
#' ## Hierarchical      ##
#' #######################
#' pop.mean  <- c(a=1.25, v.f1=4,  v.f2=3,  z=.5, sz=.3, sv=2,  t0=.3)
#' pop.scale <- c(a=.5,   v.f1=.5, v.f2=.5, z=.1, sz=.1, sv=.3, t0=.05)
#' pop.prior <- BuildPrior(
#'   dists = rep("tnorm", npar),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,  -5, -5, 0, 0, 0, 0),
#'   upper = c(5,   7,  7, 2, 2, 2, 2))
#'
#' nsubject <- 10
#' ntrial <- 1e2
#' ## You may have errors here, informing you that DDM parameters are invalid.
#' ## That is from MG's rtdists. Simply re-run simulate several times or alter
#' ## pop.prior
#' dat <- simulate(m1, nsim = ntrial, nsub = nsubject, p.prior = pop.prior)
#' dmi2 <- BindDataModel(dat, m1)
#' ps <- attr(dat, "parameters")
#' p.prior <- BuildPrior(
#'    dists = rep("tnorm", npar),
#'    p1    = pop.mean,
#'    p2    = pop.scale*5,
#'    lower = c(0, -5, -5, 0, 0, 0, 0),
#'    upper = c(5,  7,  7, 2, 2, 2, 2))
#' mu.prior <- BuildPrior(
#'    dists = rep("tnorm", npar),
#'    p1    = pop.mean,
#'    p2    = pop.scale*5,
#'    lower = c(0, -5, -5, 0, 0, 0, 0),
#'    upper = c(5,  7,  7, 2, 2, 2, 2))
#' sigma.prior <- BuildPrior(
#'    dists = rep("beta", npar),
#'    p1    = rep(1, npar),
#'    p2    = rep(1, npar),
#'    upper = rep(2, npar))
#' names(sigma.prior) <- GetPNames(m1)
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' thin <- 2
#' nchain <- npar * 3
#' hsam0 <- StartNewHypersamples(32, dmi2, p.prior, pp.prior, thin, 001, nchain)
#' ## hsam0 <- run(hsam0)
#'
#' @export
StartNewsamples <- function(nmc, data = NULL, p.prior = NULL, thin = 1,
  nchain = NULL, rp = .001) {

  model <- CheckDMI(data, p.prior, NULL, nchain)
  npar  <- length(GetPNames(model))
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
  class(out) <- c("model", "list")
  cat("\n")

  return(out)
}

##' @rdname StartNewsamples
##' @export
RestartSamples <- function(nmc, samples = NULL, thin = NULL, rp = .001,
  add = FALSE) {

  hyper <- attr(samples, "hyper")
  # model <- CheckSamples(samples, p.prior, NULL)
  # npar  <- length(GetPNames(model))
  if (is.null(samples)) stop("Use StartNewsamples")
  if (!is.null(hyper)) stop("Use RestartHypersamples")

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

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
StartManynewsamples <- function(nmc, data = NULL, p.prior = NULL, thin = 1,
                                nchain = NULL, rp = .001) {
  npar <- length(p.prior)
  if (is.null(nchain)) nchain <- 3*npar
  if (is.null(p.prior)) stop("No p.prior no new samples")

  if (is.null(attr(data[[1]], "bw"))) {
    message("No GPU attributes. Default bw = .001, using GPU 0.")
    for(i in 1:length(data)) {
      attr(data[[i]], "n.pda") <- 1e4
      attr(data[[i]], "bw") <- .001
      attr(data[[i]], "gpuid") <- 0
    }
  }

  out <- ggdmc::init_newnonhier(nmc, data, p.prior, rp, thin, nchain)

  for(i in 1:length(out)) {
    pnames <- GetPNames(attr(out[[i]]$data, "model"))
    dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
    class(out[[i]]) <- c("model", "list")
  }
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
RestartManysamples <- function(nmc, samples = NULL, thin = NULL,
  rp = .001, add = FALSE) {

  if (is.null(samples)) stop("Use StartManynewsamples")
  if (is.null(thin)) thin <- samples[[1]]$thin

  if (add) {
    out <- ggdmc::init_addnonhier(nmc, samples, rp, thin)
  } else {
    out <- ggdmc::init_oldnonhier(nmc, samples, rp, thin)
  }

  for(i in 1:length(out)) {
    pnames <- GetPNames(attr(out[[i]]$data, "model"))
    dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
    class(out[[i]]) <- c("model", "list")
  }
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
StartNewHypersamples <- function(nmc, data = NULL, p.prior = NULL,
  pp.prior = NULL, thin = 1, rp = .001, nchain = NULL) {

  nchain <- CheckHyperDMI(data, nchain) ## If nchain=NULL, change it to default
  if (is.null(p.prior)) stop("No p.prior")   ## Check priors
  if (is.null(pp.prior)) stop("No pp.prior")
  if (!is.list(pp.prior)) stop("pp.prior must be a list")
  if (length(pp.prior[[1]]) < length(pp.prior[[2]]))
    stop("Location priors must have as many or more elements than scale priors")

  if (is.null(attr(data[[1]], "bw"))) {
    message("No GPU attributes. Default bw = .001, using GPU 0.")
    for(i in 1:length(data)) {
      attr(data[[i]], "n.pda") <- 1e4
      attr(data[[i]], "bw") <- .001
      attr(data[[i]], "gpuid") <- 0
    }
  }


  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  isppriorok <- pnames %in% names(p.prior)
  islpriorok <- pnames %in% names(pp.prior[[1]])
  isspriorok <- pnames %in% names(pp.prior[[2]])
  if (!all(isppriorok)) stop("p.prior incompatible with model")
  if (!all(islpriorok)) stop("location prior incompatible with model")
  if (!all(isppriorok)) stop("scale prior incompatible with model")

  out <- init_newhier(nmc, data, p.prior, pp.prior, rp, thin, nchain)
  names(out) <- names(data)
  class(out) <- c("model", "list")
  return(out)
}

##' @rdname StartNewsamples
##' @export
RestartHypersamples <- function(nmc, samples = NULL, thin = NULL, rp = .001,
                                add = FALSE) {

  if (is.null(samples)) stop("Use StartNewHypersamples")
  if (is.null(thin)) thin <- samples[[1]]$thin

  if (add) {
    out <- init_addhier(nmc, samples, rp, thin) # add onto existing one
  } else {
    out <- init_oldhier(nmc, samples, rp, thin) # start afresh
  }
  names(out) <- names(samples)
  class(out) <- c("model", "list")
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
run_one <- function(samples, report, pm, qm, gammamult, ngroup, force,
  sampler, debug) {

  ncore <- 1
  force  <- MakeForce(samples, force)
  pnames <- GetPNames(attr(samples$data, "model"))

  if (is.null(attr(samples$data, "n.pda"))) attr(samples$data, "n.pda") <- 2^14
  if (is.null(attr(samples$data, "bw")))    attr(samples$data, "bw")    <- .01
  if (is.null(attr(samples$data, "debug"))) attr(samples$data, "debug") <- 0
  if (is.null(attr(samples$data, "gpuid"))) attr(samples$data, "gpuid") <- 0

  if (sampler == "DGMC") {
    message("Run DGMC")
    out <- run_dgmc(samples, force, report, pm, qm, gammamult, ncore, ngroup)
  } else if (sampler == "DE-MCMC") {
    message("Run DE-MCMC")
    out <- run_dmc(samples, force, report, pm, gammamult, ncore, debug)
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
  force, sampler, debug) {

  force <- MakeForce(samples[[1]], force)
  for(i in 1:length(samples)) {
    if (is.null(attr(samples[[i]]$data, "n.pda"))) attr(samples[[i]]$data, "n.pda") <- 2^14
    if (is.null(attr(samples[[i]]$data, "bw")))    attr(samples[[i]]$data, "bw")    <- .01
    if (is.null(attr(samples[[i]]$data, "debug"))) attr(samples[[i]]$data, "debug") <- 0
    if (is.null(attr(samples[[i]]$data, "gpuid"))) attr(samples[[i]]$data, "gpuid") <- 0
  }

  if (get_os() == "windows" & ncore > 1) {
    cl  <- parallel::makeCluster(ncore)
    if (sampler == "DGMC") {
      message("Run many subjects using DGMC in parallel")

      out <- parallel::parLapply(cl, samples, run_dgmc, force, report, pm, qm,
        gammamult, ncore, ngroup)
    } else if (sampler == "DE-MCMC") {
      message("Run many subjects using DE-MCMC in parallel")

      out <- parallel::parLapply(cl, samples, run_dmc, force, report, pm,
        gammamult, ncore)
    } else {
      stop("Sampler unknown")
    }
    stopCluster(cl)

  } else if (ncore > 1) {
    if (sampler == "DGMC") {
      message("Run many subjects using DGMC in parallel")

      out <- parallel::mclapply(samples, run_dgmc, force, report, pm, qm,
        gammamult, ncore, ngroup, mc.cores=ncore)
    } else if (sampler == "DE-MCMC") {
      message("Run many subjects using DE-MCMC in parallel")

      out <- parallel::mclapply(samples, run_dmc, force, report, pm,
        gammamult, ncore, mc.cores=ncore)
    } else {
      stop("Sampler unknown")
    }

  } else {
    if (sampler == "DGMC") {
      message("Run many subject using DGMC with lapply")

      out <- lapply(samples, run_dgmc, force, report, pm, qm, gammamult, ncore,
        ngroup)
    } else if (sampler == "DE-MCMC") {
      message("Run many subject using DE-MCMC lapply")

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
#' ###############
#' ## Run examples
#' ###############
#' @export
run <- function(samples, report = 1e2, ncore = 1, pm = 0, qm = 0, hpm = 0,
  gamma.mult = 2.38, ngroup = 6, force = FALSE, sampler = "DE-MCMC",
  debug = FALSE) {
  ## passing debug == TRUE to launch old migration operator

  hyper <- attr(samples, "hyper")

  if (!is.null(hyper)) {   ## hierarchical model
    if (is.null(attr(samples[[1]]$data, "bw"))) {
      message("No GPU attributes. Default bw = .001, using GPU 0.")
      for(i in 1:length(data)) {
        attr(samples[[i]]$data, "bw") <- .001
        attr(samples[[i]]$data, "gpuid") <- 0
      }
    }

    if (sampler == "DE-MCMC") {
      message("Run Hierarchical DE-MCMC")
      out <- run_hyper_dmc(samples, report, pm, hpm, gamma.mult, ncore, debug)
    } else if (sampler == "DGMC") {
      message("Run Hierarchical DGMC")
      out <- run_hyper_dgmc(samples, report, pm, hpm, qm, gamma.mult, ngroup,
        ncore)
    } else {
      out <- NULL
      message("Unknown sampler?")
    }

    phi1_tmp <- attr(out, "hyper")$phi[[1]]
    phi2_tmp <- attr(out, "hyper")$phi[[2]]
    pnames <- GetPNames(attr(samples[[1]]$data, "model"))
    dimnames(phi1_tmp) <- list(NULL, pnames, NULL)
    dimnames(phi2_tmp) <- list(NULL, pnames, NULL)
    attr(out, "hyper")$phi[[1]] <- phi1_tmp
    attr(out, "hyper")$phi[[2]] <- phi2_tmp

    for(i in 1:length(out)) {
      pnames <- GetPNames(attr(out[[i]]$data, "model"))
      dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
      class(out[[i]]) <- c("model", "list")
    }

  } else if (any(names(samples) == "theta")) { ## One subject
    out <- run_one(samples, report, pm, qm, gamma.mult, ngroup, force,
      sampler, debug)
    pnames <- GetPNames(attr(out$data, "model"))
    dimnames(out$theta) <- list(NULL, pnames, NULL)

  } else {  ## Multiple independent subjects
    out <- run_many(samples, report, ncore, pm, qm, gamma.mult, ngroup, force,
      sampler, debug)
    for(i in 1:length(out)) {
      pnames <- GetPNames(attr(out[[i]]$data, "model"))
      dimnames(out[[i]]$theta) <- list(NULL, pnames, NULL)
      class(out[[i]]) <- c("model", "list")
    }
  }

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}


##' @rdname run
##' @export
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


##' @rdname run
##' @export
run_unstuck <- function(samples, nmc=NA, report=1e2,
  cut=10, nbad=0, max.try=100, pm=0, qm = 0, gammamult=2.38, verbose=FALSE,
  end.no.migrate=FALSE, sampler = "DE-MCMC", debug = FALSE) {
  # Repeats sampling until <= nbad stuck chains as defined by cut or max.try
  # If samples has no content fills it in then repeatedly gets new sets of nmc
  # samples (nmc can be specified or taken from samples). If end.no.migrate
  # runs one final time with migration off.
  if (is.null(samples$theta)) stop("For multiple subjects use run_many_unstuck")
  nchain <- samples$n.chain

  if (any(is.na(samples$theta[,,2]))) {
     samples <- run_one(samples = samples, report = report, pm=pm, qm = qm,
       gammamult = gammamult, ngroup = 6, force = FALSE,
        sampler = sampler, debug = debug)
  }

  if ( is.na(nmc) ) nmc <- samples$nmc
  try.num <- 1

  repeat {
     cat(paste("\nTry", try.num, "\n"))
     if ( length(PickStuck(samples, verbose = verbose)) <= nbad ) break

     samples <- run_one(RestartSamples(nmc, samples),
       report=report, pm=pm, qm = qm, gammamult = gammamult,
       ngroup = 6, force = FALSE, sampler = sampler, debug = debug)
     if (try.num >= max.try) break
     try.num <- try.num + 1
  }

  if (end.no.migrate) {
     samples <- run_one(RestartSamples(nmc, samples),
       report=report, pm=0, qm = qm, gammamult = gammamult,
       ngroup = 6, force = FALSE, sampler = sampler, debug = debug)
  }

  return(samples)
}

##' @rdname run
##' @export
run_converge <- function(samples, nmc, report=1e2,
  gammamult=2.38, cut=1.1, max.try=100, minN=NA, meanN=NA, transform=TRUE,
  autoburnin=FALSE, split=TRUE, verbose=FALSE, sampler = "DE-MCMC",
  debug = FALSE)
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

  if ( is.null(samples$theta) ) stop("For multiple subjects use run_many_unstuck")

  if ( any(is.na(samples$theta[,,2])) ) {
    samples <- run_one(samples = samples, report = report, pm=0, qm=0,
      gammamult = gammamult, ngroup = 6, force = FALSE,
      sampler = sampler, debug = debug)

    gd <- ggdmc:::gelman.diag.mpsrf(theta.as.mcmc.list(samples, split=split),
      autoburnin=autoburnin, transform=transform)
    if (verbose) cat(paste("MPSRF: ", gd, "\n"))

  } else {
    gd <- Inf
  }

  if ( !is.na(minN) & (gd <= cut) ) {
    okN <- do.call(nfun, list(ggdmc:::effectiveSize.dmc(samples))) > minN
  } else {
    okN <- TRUE
  }

  if ( (gd > cut) | !okN ) {  # Do more sampling
      if (is.na(nmc)) nmc <- samples$nmc
      try.num <- 0
      effectiveN <- NA
      repeat {
        samples <- run_one(samples = samples, report = report, pm=0, qm=0,
          gammamult = gammamult, ngroup = 6, force = FALSE,
          sampler = sampler, debug = debug)

        gd <- gelman.diag.mpsrf(theta.as.mcmc.list(samples,split=split),
          autoburnin=autoburnin,transform=transform)

        if ( try.num > 0 ) {
          # shorter  <- RestartSamples(nmc=0, samples=samples, remove=1:nmc, add=TRUE)
          # gd.short <- gelman.diag.mpsrf(theta.as.mcmc.list(shorter,split=split),
          #   autoburnin=autoburnin,transform=transform)
          # if (gd.short < gd) {
          #   samples <- shorter
          #   gd <- gd.short
          #   if (verbose) cat(paste("Discarding initial", nmc,"samples.\n"))
          # }
        }

        if ( try.num >= max.try ) break
        try.num <- try.num + 1
        if ( gd <= cut ) {

          if ( is.na(minN) ) {
            break
          } else {
            effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples)))
            if (effectiveN > minN) break
          }

        }

        if (verbose) print(paste("N = ", dim(samples$theta)[3],
          "Multivariate psrf achieved = ", format(gd, digits=3)))
      }
  }

  if (verbose) {
      print(paste("Final multivariate psrf =",gd))
      cat("Effective sample size\n")
      print(effectiveSize.dmc(samples))
  }
  return(samples)

}


##' @rdname run
##' @export
run_manyunstuck <- function(samples, nmc=NA, report=1e2, ncore = 1,
  cut = 10, nbad = 0, max.try=20, pm = 0, qm = 0,
  gammamult = 2.38, verbose = TRUE,
  end.no.migrate=FALSE, sampler = "DE-MCMC", debug = FALSE)
{
  if ( !is.null(samples$theta) ) stop("Use run_unstuck")
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc

  if ( any(is.na(samples[[1]]$theta[,,2])) ) {
    cat("Getting initial set of samples\n")
    samples <- run(samples = samples, report = report, ncore = ncore,
      pm = pm, qm = qm, hpm = 0, gamma.mult = gammamult)
  }

  if ( get.os() == "windows" & ncore > 1) {
    cl  <- parallel::makeCluster(ncore)
    samples <- parallel::parLapply(cl, samples, run_unstuck,
      nmc, report, cut, nbad, max.try, pm, qm, gammamult, verbose,
      end.no.migrate, sampler, debug)
    stopCluster(cl)
  } else if ( ncore > 1) {


    samples <- parallel::mclapply(samples, run_unstuck,
      nmc, report, cut, nbad, max.try, pm, qm, gammamult, verbose,
      end.no.migrate, sampler, debug, mc.cores = ncore)
  } else {
    samples <- lapply(samples, run_unstuck,
      nmc, report, cut, nbad, max.try, pm, qm, gammamult, verbose,
      end.no.migrate, sampler, debug)
  }

  for(i in 1:length(samples)) {
    pnames <- GetPNames(attr(samples[[i]]$data, "model"))
    dimnames(samples[[i]]$theta) <- list(NULL, pnames, NULL)
    class(samples[[i]]) <- c("model", "list")
  }

  class(samples) <- c("model", "list")
  return(samples)
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

