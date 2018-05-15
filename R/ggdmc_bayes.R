#' Hierarchical prior, likelihoods and posterior
#'
#' \code{h.log.likelihood.dmc} computes log-likelihood of subject parameters
#' ps under population distribuiton p.prior
#'
#' \code{h.log.posterior.dmc} computes log-likelihood of subject parameters ps
#' under population distribution p.prior,given population parameters pp (phi)
#' with priors pp.prior
#'
#' @examples
#' pp <- list( c(1.80, -0.86, 6.47, 0.64, 0.28, 0.14, 0.79),
#'             c(0.57, 1.40, 1.46, 0.17, 1.65, 0.46, 0.60))
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=0.55, sz=0.15, sv=0.32, t0=0.25)
#' pop.scale <- c(a=0.10, v.f1=.8,   v.f2=.5,   z=0.1,  sz=0.05, sv=0.05, t0=0.05)
#' pop.prior <- ggdmc::prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1=pop.mean,
#'   p2=pop.scale,
#'   lower=c(0,-5, -5, 0, 0,   0, 0),
#'   upper=c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' p.prior <- ggdmc:::prior.p.dmc(
#' dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'      p1=pop.mean,
#'      p2=pop.scale*5,
#'      lower=c(0,-5, -5, 0, 0, 0, 0),
#'      upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' mu.prior <- ggdmc:::prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' sigma.prior <- ggdmc:::prior.p.dmc(
#'   dists = rep("beta", length(p.prior)),
#'   p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),p2=c(1,1,1,1,1,1,1),
#'   upper=c(2,2,2,2,2, 2, 2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' ggdmc:::h.summed.log.prior(pp, pp.prior)
#'
#' ## LBA model
#'
#  ## Population distribution, rate effect on F
#' pop.mean <- c(A=.4, B.r1=.6, B.r2=.8, t0=.3,
#'   mean_v.f1.true=1.5, mean_v.f2.true=1, mean_v.f1.false=0,
#'   mean_v.f2.false=0, sd_v.true = .25)
#' pop.scale <-c(A=.1, B.r1=.1, B.r2=.1, t0=.05,
#'   mean_v.f1.true=.2, mean_v.f2.true=.2, mean_v.f1.false=.2,
#'   mean_v.f2.false=.2, sd_v.true = .1)
#'
#' p.prior <- ggdmc:::prior.p.dmc(
#'   dists = rep("tnorm", 9),
#'   p1 = pop.mean,
#'   p2 = pop.scale,
#'   lower = c(0,0,0, .1, NA,NA,NA,NA,0),
#'   upper = c(NA,NA,NA, 1, NA,NA,NA,NA,NA))
#'
#' dat <- ggdmc:::h.simulate.dmc(model, ns=8, n=250, p.prior=pop.prior)
#' ps <- round( attr(dat, "parameters"), 2)
#'
#' pp <- list( c(1.80, -0.86, 6.47, 0.64, 0.28, 0.14, 0.79, 0.2, 0.3),
#'              c(0.57, 1.40, 1.46, 0.17, 1.65, 0.46, 0.60, 0.3, 0.1))
#' ggdmc:::h.log.likelihood.dmc(ps, pp, p.prior)
#'
#' res <- assign.pp(pp, prior)
#'
#' @export
h.summed.log.prior <- function (pp, pp.prior) {
  p.vector1 <- pp[[1]]
  p.vector2 <- pp[[2]]
  names(p.vector1) <- names(pp.prior[[1]])
  names(p.vector2) <- names(pp.prior[[2]])
  sum(logprior(p.vector1, pp.prior[[1]])) +
    sum(logprior(p.vector2, pp.prior[[2]]))
}

#' Add log-likelihoods across subjects at the hyper level
#'
#' @param ps a nsubject x npar matrix
#' @param pp a temporary pp.prior values extracted from phi. The temporary
#' pp.prior has no attached parameter names. phi is two-element list. First
#' element is a location array of nchain x npar x nmc; second element is a
#' scale array of nchain x npar x nmc.
#' @param p.prior prior distributions at data level
#'
#' @export
SummedHyperloglikelihood <- function(ps, pp, p.prior) {
  term1 <- apply(ps, 1, ggdmc::summedlogpriorNV, ggdmc:::assign.pp(pp,p.prior))
  term2 <- apply(term1, 2, sum)
  return(term2)
}

#' @rdname h.summed.log.prior
#' @export
h.log.posterior.dmc <- function(ps,p.prior,pp,pp.prior)
{
  ## 1. sum over subjects of likelihood of pp (pop pars) given ps (subjects pars)
  ## 2. prior probability of pp
  sum(h.log.likelihood.dmc(ps, pp, p.prior)) + h.summed.log.prior(pp, pp.prior)
}


#' @rdname h.summed.log.prior
#' @export
assign.pp <- function(pp, p.prior)
  # Slot pp values into p.prior
{
  for (i in 1:length(p.prior))
    p.prior[[i]][1:2] <- c(pp[[1]][i], pp[[2]][i])
  p.prior
}


#' @rdname h.summed.log.prior
#' @export
log.posterior.dmc <- function(p.vector, p.prior, data, min.like=1e-10)
  # Summed log posterior likelihood
{
  sum (log.likelihood(p.vector,data, min.like=min.like)) +
    summed.log.prior(p.vector, p.prior)
}

#' @rdname h.summed.log.prior
#' @export
log.likelihood <- function(p.vector, dmi, min.like=1e-10) {
  type <- attr(attr(dmi, "model"), "type")

  if (type == "rd") {
    out <- likelihood_rd(p.vector, dmi, min.like)
  } else if (type == "norm") {
    out <- likelihood_norm(p.vector, dmi)
  } else if (type == "norm_pda") {
    out <- likelihood.norm_pda(p.vector, dmi, min.like)
  } else if (type == "plba1_gpu") {
    out <- likelihood.plba1_gpu(p.vector, dmi, min.like)
  } else {
    stop("log.likelihood type not supported")
  }
  return(log(out))
}


