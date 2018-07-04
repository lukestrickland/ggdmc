#' Which Chains Get Stuck
#'
#' Calculate each chain separately for the mean (across many MCMC iterations)
#' of posterior log-likelihood. If the difference of the means and
#' the median (across chains) of the mean of posterior is greater than the
#' \code{cut}, chains are considered stuck. The default value for \code{cut}
#' is 10. \code{unstick} manually removes stuck chains from posterior samples.
#'
#' @param samples posterior samples
#' @param hyper a boolean switch indicating if \code{samples} has hyper
#' parameters
#' @param cut a criterion deciding if a chain is stuck.
#' @param start start to evaluate from which iteration.
#' @param end end at which iteration for evaeuation.
#' @param verbose a boolean switch to print more information
#' @param digits print how many digits. Default is 2
#' @param bad an index vector indicating chain index. R index
#' @return \code{PickStuck} gives an index vector; \code{unstick} gives a DMC
#' sample.
#' @export
#' @importFrom stats median
#' @examples
#' model <- ggdmc::BuildModel(
#' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
#' match.map = list(M = list(s1 = 1, s2 = 2)),
#' factors   = list(S = c("s1", "s2")),
#' constants = c(st0 = 0, sd_v = 1),
#' responses = c("r1", "r2"),
#' type      = "norm")
#'
#' p.vector <- c(A = .75, B = .25, t0 = .2, mean_v.true = 2.5, mean_v.false = 1.5)
#'
#' p.prior <- ggdmc:::prior.p.dmc(
#'   dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
#'   p1    = c(A = .3, B = .3, t0 = 1, mean_v.true = 1, mean_v.false = 0),
#'   p2    = c(1, 1,   1, 3, 3),
#'   lower = c(0,  0,  0, NA, NA),
#'   upper = c(NA,NA,  1, NA, NA))
#'
#' \dontrun{
#' dat <- ggdmc:::BindDataModel(ggdmc:::simulate.model(model, p.vector, 512), model)
#' s0  <- ggdmc::run(ggdmc::samples(1024, p.prior, dat, thin = 2))
#' s1  <- ggdmc::run(ggdmc::samples(1024, p.prior, samples = s0))
#' bad <- ggdmc::PickStuck(s1)
#' s2   <- ggdmc::unstick(s1, bad)
#' s3   <- ggdmc:::theta.as.mcmc.list(s2)
#' plot(s3)
#' }
PickStuck <- function(samples, hyper = FALSE, cut = 10, start = 1, end = NA,
                      verbose = FALSE, digits = 2)
{
  # hyper = FALSE
  # cut = 10
  # start = 1
  # end = NA
  # verbose = FALSE
  # digits = 2
  if (hyper) { # MG: TODO: Haven't split up hyper$pll yet
    hyper <- attr(samples, "hyper")
    if (is.na(end)) end <- dim(hyper$phi[[1]])[3]
    if (end <= start) stop("End must be greater than start")
    mean.ll <- apply(hyper$h_log_likelihoods[start:end,] +
        hyper$h_summed_log_prior[start:end,], 2, mean)
    names(mean.ll) <- 1:length(mean.ll)
  } else {
    if (is.na(end)) end <- samples$nmc
    if (end <= start) stop("End must be greater than start")

    mean.ll <- apply(samples$log_likelihoods[start:end,] +
                     samples$summed_log_prior[start:end,], 2, mean)
    names(mean.ll) <- 1:length(mean.ll)
  }

  dev <- -(sort(mean.ll) - median(mean.ll))
  bad <- as.numeric(names(dev)[dev > cut])

  if (verbose) {
    cat("Deviation of mean chain log-likelihood from median of means\n")
    print(round(dev,digits))
    cat("Bad chains: ")
    if (length(bad)==0) { cat("None\n") } else { cat(bad, "\n") }
  }
  return(bad)
}

#' @rdname PickStuck
#' @export
unstick <- function(samples, bad)
{
  nchain <- samples$n.chains
  if (length(bad) > 0) {
    if (!all(bad %in% 1:nchain)) stop(paste("Index of bad chains must be in 1 to ", nchain))
    samples$theta            <- samples$theta[-bad,,]
    samples$summed_log_prior <- samples$summed_log_prior[,-bad]
    samples$log_likelihoods  <- samples$log_likelihoods[,-bad]
    samples$n.chains         <- samples$n.chains - length(bad)
  }
  return(samples)
}


#' Convert Theta to a mcmc List
#'
#' \code{theta.as.mcmc.list} extracts data-level parameter array, \code{theta},
#' from a DMC sample and convert it to a \pkg{coda} mcmc.list.
#' \code{phi.as.mcmc.list} extracts hyper-level parameter array, \code{phi}
#' from a DMC sample and convert it to a \pkg{coda} mcmc.list.
#' When \code{split} switch is TRUE, the function doubles number of chains,
#' first and second half.
#'
#' @param x a sample list
#' @param start start iteration
#' @param end end iteraton
#' @param split whether to divide one MCMC sequence into two sequences.
#' @importFrom coda mcmc mcmc.list
#' @export
theta.as.mcmc.list <- function(x, start = 1, end = NA, split = FALSE,
  subchain = FALSE, nsubchain = 3, thin = NA) {

  if (is.na(thin)) thin <- x$thin
  nchain <- x$n.chains

  if (subchain) {
    message("pMCMC diagnosis randomly select a subset of chains: ", appendLF = FALSE)
    chain.idx <- base::sample(1:nchain, nsubchain)
    cat(chain.idx, "\n")
    nchain <- nsubchain
  } else {
    chain.idx <- 1:nchain
  }

  if (is.na(end)) end <- x$nmc

  lst <- vector(mode = "list", length = nchain * ifelse(split, 2, 1))
  idx <- start:end

  if (split) is.in <- !as.logical(idx %% 2) else is.in <- rep(TRUE, length(idx))

  if (split) {
    not.is.in <- !is.in
    if ( sum(not.is.in) > sum(is.in) ) {
      not.is.in[1:2][not.is.in[1:2]] <- FALSE
    } else if ( sum(not.is.in) < sum(is.in) ) {
      is.in[1:2][is.in[1:2]] <- FALSE
    }
  }

  for (i in 1:nchain) {
    lst[[i]] <- coda::mcmc( t(x$theta[chain.idx[i], , idx[is.in]]),
      thin = thin)
  }

  if (split) {
    for (i in 1:nchain) {
      lst[[i+nchain]] <- mcmc( t(x$theta[chain.idx[i], , idx[not.is.in]]),
        thin = thin)
    }
  }

  mcmc.list(lst)
}

##' @rdname theta.as.mcmc.list
##' @importFrom coda mcmc mcmc.list
##' @export
phi.as.mcmc.list <- function(x, start = 1, end = NA, split = FALSE,
  subchain = FALSE, nsubchain = 3) {

  thin   <- x$thin   ## x == hyper
  nchain <- x$n.chains

  if (subchain) {
    message("pMCMC diagnosis randomly select a subset of chains: ", appendLF = FALSE)
    chain.idx <- base::sample(1:nchain, nsubchain)
    cat(chain.idx, "\n")
    nchain <- nsubchain
  } else {
    chain.idx <- 1:nchain
  }

  # Parameters that are not constants
  ok1 <- lapply(x$pp.prior,function(x){
    lapply(x,function(y){attr(y, "dist") != "constant"})})
  ok2 <- paste(names(ok1[[2]])[unlist(ok1[[2]])], "h2", sep=".")
  ok1 <- paste(names(ok1[[1]])[unlist(ok1[[1]])], "h1", sep=".")
  if ( is.na(end) ) end <- x$nmc

  lst <- vector(mode = "list", length = nchain)
  indx <- start:end

  if (split) is.in <- !as.logical(indx %% 2) else is.in <- rep(TRUE, length(indx))
  if (split) {
    not.is.in <- !is.in
    if ( sum(not.is.in) > sum(is.in) ) {
      not.is.in[1:2][not.is.in[1:2]] <- FALSE
    } else if ( sum(not.is.in) < sum(is.in) ) {
      is.in[1:2][is.in[1:2]] <- FALSE
    }
  }

  for (i in 1:nchain) {
    tmp1 <- t(x$phi[[1]][chain.idx[i], , indx[is.in]]) ## nmc x npar matrix
    ## attach parnames with h1
    dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1", sep=".")
    tmp1 <- tmp1[,ok1]  ## exclude constant parameter

    ## same thing on scale
    tmp2 <- t(x$phi[[2]][chain.idx[i], , indx[is.in]])
    dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2", sep=".")
    tmp2 <- tmp2[,ok2]

    ## Remove cases with !has.sigma
    tmp2 <- tmp2[,!apply(tmp2, 2, function(x){all(is.na(x))})]
    lst[[i]] <- coda::mcmc(cbind(tmp1, tmp2), thin = thin)
  }

  if (split) {
    for (i in 1:nchain) {
      tmp1 <- t(x$phi[[1]][chain.idx[i], , indx[not.is.in]])
      dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1", sep=".")
      tmp1 <- tmp1[,ok1]
      tmp2 <- t(x$phi[[2]][chain.idx[i], , indx[not.is.in]])
      dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2", sep=".")
      tmp2 <- tmp2[,ok2]
      # Remove cases with !has.sigma
      tmp2 <- tmp2[,!apply(tmp2,2,function(x){all(is.na(x))})]
      lst[[i + nchain]] <- coda::mcmc(cbind(tmp1,tmp2), thin = thin)
    }
  }

  return(coda::mcmc.list(lst))
}

#' Gelman and Rubin Convergence Diagnostic
#'
#' \code{gelman} calls \pkg{coda} gelman.diag to get R hats for one
#' or a list of subjects. It can calculate at the data or hyper level.
#' R hat for one or list of subjects or hyper level.
#' split doubles the number of chains by spliting into 1st and 2nd halves.
#'
#' @param x a DMC sample
#' @param hyper a switch to extract hyper attribute and calculate it
#' @param digits print out how many digits
#' @param start start iteration
#' @param autoburnin turn on auto burnin
#' @param transform turn on transform
#' @param end end iteraton
#' @param ... arguments passing to \code{coda} gelman.diag.
#' @importFrom coda gelman.diag
#' @examples
#' m1 <- model.dmc(
#'      p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'      match.map = list(M=list(s1="r1",s2="r2")),
#'      factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'      constants = c(st0=0,d=0),
#'      responses = c("r1","r2"),
#'      type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=.55,  sz=.15, sv=.32, t0=.25)
#' pop.scale <- c(a=.10,  v.f1=.8,   v.f2=.5,   z=0.1,  sz=.05, sv=.05, t0=.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=8, p.prior=pop.prior)
#' mdi1 <- BindDataModel(dat, m1)
#' ps   <- attr(dat,  "parameters")
#'
#' ### FIT RANDOM EFFECTS
#' p.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' mu.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", length(p.prior)),
#'   p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),p2=c(1,1,1,1,1,1,1),
#'   upper=c(2,2,2,2,2, 2, 2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' hsamples0 <- h.samples.dmc(nmc=10, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0)
#' gelman.diag.dmc(hsamples0, hyper=TRUE)
#' gelman.diag.dmc(hsamples0)
#' @export
gelman <- function(x, hyper = FALSE, start = 1, end=NA, confidence = 0.95,
  transform=TRUE, autoburnin = FALSE, multivariate = TRUE, split = TRUE,
  subchain = FALSE, nsubchain = 3, digits = 2, verbose = TRUE)
{
  if ( hyper ) {
    message("Diagnosing the hyper parameters, phi")
    hyper <- attr(x, "hyper")
    thin <- hyper$thin
    if (is.na(end)) end <- hyper$nmc

    mcmclist <- ggdmc:::phi.as.mcmc.list(hyper, start, end, split,
      subchain, nsubchain)

    out <- coda::gelman.diag(mcmclist, confidence, transform, autoburnin,
      multivariate)

  } else {
    ## if x is one subject samples, we should found an elemnet called theta
    if ( !is.null(x$theta) ) {
      message("Diagnosing a single participant, theta")
      if (is.na(end)) end <- x$nmc
      mcmclist <- ggdmc:::theta.as.mcmc.list(x, start, end, split, subchain,
        nsubchain)
      out <- coda::gelman.diag(mcmclist, confidence, transform, autoburnin,
        multivariate)
    } else {
      message("Diagnosing theta for many participants separately")
      nsub <- length(x)
      out   <- vector(mode = "list", length = nsub)
      start <- rep(start, length.out = nsub)
      end   <- rep(end, length.out = nsub)
      subchain <- rep(subchain, length.out = nsub)
      nsubchain   <- rep(nsubchain, length.out = nsub)
      names(out) <- names(x)

      for (i in 1:nsub) {
        if ( is.na(end[i]) ) end[i] <- x[[i]]$nmc
        mcmclist <- ggdmc:::theta.as.mcmc.list(x[[i]], start[i], end[i], split,
          subchain[i], nsubchain[i])
        out[[i]] <- gelman.diag(mcmclist, confidence, transform, autoburnin,
          multivariate)
      }

      tmp <- unlist(lapply(out, function(x){x$mpsrf}))

      if (verbose) {
        print(round(sort(tmp), digits))
        cat("Mean\n")
        print(round(mean(tmp), digits))
      }

      invisible(out)
    }
  }

  return(out)
}


##' @rdname gelman
##' @export
hgelman <- function(x, start = 1, end = NA, confidence = 0.95, transform = TRUE,
  autoburnin = FALSE, multivariate = TRUE, split = TRUE, subchain = FALSE,
  nsubchain = 3, digits = 2, verbose=TRUE, ...) {

  step1 <- lapply(gelman(x, start = start, end = end, confidence = confidence,
    transform = transform, autoburnin = autoburnin, multivariate = multivariate,
    split=split, subchain = subchain, nsubchain = nsubchain, verbose=FALSE),
    function(x){x$mpsrf})

  out <- base::sort(unlist(step1)) ## non-hyper

  if ( any(names(attributes(x)) == "hyper") ) {
    hyper.gd <- gelman(x, hyper = TRUE, start = start, end = end,
      confidence = confidence, transform = transform, autoburnin = autoburnin,
      multivariate = multivariate, split=split, subchain = subchain,
      nsubchain = nsubchain, verbose=FALSE)

    out <- c(hyper.gd$mpsrf, out)
    names(out) <- c("hyper", names(step1))
  }


  if (verbose) print(round(out, digits))
  invisible(out)
}

#' Effective Sample Size for Estimating the Mean
#'
#' \code{effectiveSize.dmc} calls \pkg{coda} effectiveSize to effective size
#' for either single or multiple subjects. It can calculate at the data or
#' hyper level, too.
#'
#' @param x a DMC sample
#' @param hyper a switch to extract hyper attribute and calculate it
#' @param digits print out how many digits
#' @param start start iteration
#' @param end end iteraton
#' @importFrom coda effectiveSize
#' @export
#' @examples
#' m1 <- model.dmc(
#'      p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'      match.map = list(M=list(s1="r1",s2="r2")),
#'      factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'      constants = c(st0=0,d=0),
#'      responses = c("r1","r2"),
#'      type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=.55,  sz=.15, sv=.32, t0=.25)
#' pop.scale <- c(a=.10,  v.f1=.8,   v.f2=.5,   z=0.1,  sz=.05, sv=.05, t0=.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=4, p.prior=pop.prior)
#' mdi1 <- BindDataModel(dat, m1)
#' ps   <- attr(dat,  "parameters")
#' ### FIT RANDOM EFFECTS
#' p.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' mu.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", length(p.prior)),
#'   p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),p2=c(1,1,1,1,1,1,1),
#'   upper=c(2,2,2,2,2, 2, 2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' hsamples0 <- h.samples.dmc(nmc=10, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0)
#' es <- effectiveSize.dmc(hsamples0)
#' hes <- effectiveSize.dmc(hsamples0, hyper=TRUE)
effectiveSize.dmc <- function(x, hyper=FALSE, digits=0,start=1,end=NA)
{
  if (hyper) {
    hyper0 <- attr(x, "hyper")
    if (is.na(end)) end <- hyper0$nmc
    round(effectiveSize(phi.as.mcmc.list(attr(x, "hyper"), start=start,
      end=end)), digits)
  } else {
      if (!is.null(x$theta)) {
        if (is.na(end)) end <- x$nmc
        round(effectiveSize(theta.as.mcmc.list(x, start=start,end=end)),digits)
      } else
          lapply(x, function(xx){
            if (is.na(end)) end <- xx$nmc
            round(effectiveSize(theta.as.mcmc.list(xx, start=start, end=end)),
              digits)
          })
  }
}

#' Summarise a DMC Sample with One Participant
#'
#' Call coda package to summarise the model parameters in a DMC samples
#'
#' @param object a model samples
#' @param start summarise from which MCMC iteration. Default uses the first
#' iteration.
#' @param end summarise to the end of MCMC iteration. For example, set
#' \code{start=101} and \code{end=1000}, instructs the function to calculate
#' from 101 to 1000 iteration. Default uses the last iteration.
#' @param ... other arguments
#' @keywords summary
#' @export
#' @examples
#' m1 <- model.dmc(
#' p.map=list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#' constants=c(st0=0,d=0),
#' match.map=list(M=list(s1="r1",s2="r2")),
#' factors=list(S=c("s1","s2")),
#' responses=c("r1","r2"),
#' type="rd")
#'
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1=c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2=c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower=c(0,-5, 0, 0, 0, 0),
#'   upper=c(5, 7, 2, 2, 2, 2))
#'
#' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#'
#' dat1 <- simulate(m1, nsim=30, p.vector=pVec)
#' mdi1 <- BindDataModel(dat1, m1)
#'
#' samples0 <- samples(nmc=100, p.prior=p.prior, data=mdi1)
#' samples0 <- run(samples0, p.migrate=.05)
#' gelman.diag.dmc(samples0)
#' class(samples0)
#' ## [1] "dmc"
#' summary(samples0)
summary.model <- function(object, hyper = FALSE, start = 1, end = NA,
  hmeans = FALSE, hci = FALSE, ..., digits = 2) {

  if (hyper || hmeans || hci) {
    message("Random-effec model with multiple participants")
    hyper <- attr(object, "hyper")
    if (is.na(end)) end <- hyper$nmc
    mcmclist <- phi.as.mcmc.list(hyper, start=start, end=end)
    tmp <- summary(mcmclist)

    if (hmeans) {
      out <- matrix(tmp$statistics[, "Mean"], nrow = 2,
        dimnames = list(c("h1", "h2"), object[[1]]$p.names))
    } else if (hci) {
      tmp1 <- tmp$quantiles[,c(1,3,5)]
      tmp1 <- round(cbind(tmp1[1:(dim(tmp1)[1]/2),],
        tmp1[((dim(tmp1)[1]/2)+1):dim(tmp1)[1],]),2)
      dimnames(tmp1) <- list(unlist(strsplit(
        dimnames(tmp1)[[1]],".h1")),
        paste0(c("MEAN: ","",""," SD:  ","",""),dimnames(tmp1)[[2]]) )
      out <- round(tmp1, digits)
    } else {
      out <- tmp
    }

  } else if (!is.null(object$theta)) {
    message("Single Participant")

    if (is.na(end)) end <- object$nmc
    out <- summary(theta.as.mcmc.list(object, start = start, end = end))
  } else {
    message("Summary each participant separately")

    out   <- vector(mode = "list", length = length(object))
    start <- rep(start, length.out = length(object))
    end <- rep(end, length.out = length(object))
    names(out) <- names(object)
    for (i in 1:length(object)) {
      if ( is.na(end[i]) ) end[i] <- object[[i]]$nmc
      out[[i]] <- summary(theta.as.mcmc.list(object[[i]], start = start[i],
        end = end[i]))
    }
    tmp <- t(data.frame(lapply(out, function(x){x[[1]][, 1]})))
    tmp <- rbind(tmp, apply(tmp, 2, mean))
    row.names(tmp) <- c(names(object), "Mean")
    print(round(tmp, digits))
    invisible(out)
  }

  return(out)

}



### Priors, likelihoods and model selection

pll.dmc <- function(samples,hyper=FALSE,start=1,end=NA,prior=FALSE,
  like=FALSE,subject=NA)
  # extracts posterior log-likelihoods
{
  if ( hyper ) {
    hyper <- attr(samples,"hyper")
    if (is.null(hyper))
      stop("There are no hyper-parameters to plot.")
    if ( is.na(end) ) end <- hyper$nmc
    if ( end <= start )
      stop("End must be greater than start")
    if (prior) out <- hyper$h_summed_log_prior[start:end,] else
    if (like) out <-  hyper$h_log_likelihoods[start:end,] else
    out <- hyper$h_summed_log_prior[start:end,] +
           hyper$h_log_likelihoods[start:end,]
    dimnames(out) <- list(start:end,1:dim(out)[2])
  } else {
    if ( is.null(samples$theta) & !is.na(subject) )
      samples <- samples[[subject]]
    if ( is.null(samples$theta) ) { # multiple subjects
      if ( is.na(end) ) end <- samples[[1]]$nmc
      if ( end <= start ) stop("End must be greater than start")
      out <- lapply(samples,function(x){
      if (prior) out <- x$summed_log_prior[start:end,] else
        if (like) out <- x$log_likelihoods[start:end,] else
          out <- x$summed_log_prior[start:end,] + x$log_likelihoods[start:end,]
        dimnames(out) <- list(start:end,1:dim(out)[2])
        out
     })
    } else { # single subejct
      if ( is.na(end) ) end <- samples$nmc
      if ( end <= start ) stop("End must be greater than start")
      if (prior) out <- samples$summed_log_prior[start:end,] else
      if (like) out <-  samples$log_likelihoods[start:end,] else
      out <- samples$summed_log_prior[start:end,] +
             samples$log_likelihoods[start:end,]
      dimnames(out) <- list(start:end,1:dim(out)[2])
    }
  }
  out
}

#' Calculate Dstats of DDM Density
#'
#' Calculates and returns list with mean (meanD), variance (varD) and min
#' (minD) of Deviance, deviance of mean theta (Dmean) and (optionally) matrix
#' of deviances for each theta.
#'
#' This function is still under tweaking. Please use DMC's \code{Dstat.dmc},
#' instead.
#'
#' @param samples a DMC sample/object
#' @param save a save switch
#' @param fast choose different calculation routine
#' @keywords Dstats
#' @importFrom stats var
#' @export
Dstats.dmc <- function(samples,save=FALSE,fast=TRUE)
  # Calculates and returns list with mean (meanD), variance (varD) and min (minD)
  # of Deviance, deviance of mean theta (Dmean) and (optionally) matrix of
  # deviances for each theta.
{
  if (fast) {
    summed_ll <- apply(samples$log_likelihoods,c(1,2),sum)
    D <- -2*samples$log_likelihoods

  } else D <- apply(samples$theta,c(3,1),function(x){
    -2*sum(log(likelihood.dmc(x,samples$data)))
  })

  mtheta <- apply(samples$theta,2,mean)
  Dmean <- -2*sum(log(likelihood.dmc(mtheta,samples$data)))
  minD <- min(D)
  if (save)
    list(np=dim(samples$theta)[2],meanD=mean(D),varD=var(as.vector(D)),
      minD=min(D),Dmean=Dmean,D=D) else
        list(np=dim(samples$theta)[2],meanD=mean(D),varD=var(as.vector(D)),
          minD=min(D),Dmean=Dmean,D=NA)
}


pd.dmc <- function(ds)
  # effective number of parameters calculated by mean, min and var methods
{
  list(Pmean=ds$meanD-ds$Dmean,Pmin=ds$meanD-ds$minD,Pvar=ds$varD/2)
}

#' @importFrom graphics plot abline
posterior.lr.dmc <- function(D1,D2,main="",
                             plot=FALSE,plot.density=TRUE)
  # Aitkins posterior deviance likelihood ratio test
{
  if (is.list(D1)) D1 <- D1$D
  if (is.list(D2)) D2 <- D2$D
  n <- pmin(length(D1),length(D2))
  dD <- D1[1:n]-D2[1:n]
  if (plot) if (plot.density)
    plot(density(dD),xlab="D1-D2",main=main) else
      hist(dD,breaks="fd",xlab="D1-D2",main=main)
  if (min(dD)<0 & max(dD)>0) abline(v=0)
  c(pD1=mean(dD<0))
}

#' @export
IC.dmc <- function(ds=NULL,samples=NULL,DIC=FALSE,fast=TRUE,use.pd=NA)
  # Calcualte IC1 (i.e., BPIC, default) or DIC
{
  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics") else
      ds <- Dstats.dmc(samples,TRUE,fast)
    pds <- pd.dmc(ds)
    if (is.na(use.pd)) {
      if (ds$minD < ds$Dmean) pd <- pds$Pmin else pd <- pds$Pmean
    } else {
      if (use.pd %in% names(pds))
        pd <- pds[[use.pd]] else
          stop(paste("use.pd must be one of:",paste(names(pds),collapse=",")))
    }
    if (DIC) ds$meanD+pd else ds$meanD+2*pd
}

#' @export
wIC.dmc <- function(ds=NULL,samples=NULL,
  DIC=FALSE,fast=TRUE,use.pd=NA,...)
  # Calculate weights for a set of models
{

  ICs <- function(samples,DIC=FALSE,fast=TRUE,use.pd=NA)
    IC.dmc(samples=samples,DIC=DIC,fast=fast,use.pd=use.pd)

  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics")
  if (is.null(samples))
    ics <- unlist(lapply(ds,IC.dmc,DIC=DIC,fast=fast,use.pd=use.pd)) else
      ics <- unlist(lapply(ds,ICs,DIC=DIC,fast=fast,use.pd=use.pd))
    d = ics-min(ics)
    w = exp(-d/2)/sum(exp(-d/2))
    if (!is.null(names(ics))) mnams=names(ics) else
      mnams <- paste("M",1:length(d),sep="")
    print(t(matrix(c(d,w),nrow=length(d),
      dimnames=list(mnams,c("IC-min","w")))),...)
}

#' @export
waic <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...)
  # Calclaute WAIC
{
  out <- waic(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
  if (mc_se) {
    n.chains <- dim(trial_log_likes)[2]
    waics <- vector(mode="list",length=n.chains)
    for (i in 1:n.chains) waics[[i]] <- waic(trial_log_likes[,i,])
    mc <- sd(unlist(lapply(waics,function(x){x$waic})))/sqrt(n.chains)
    print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic,
      mc_se_waic=mc),...)
  } else
    print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic),...)
  if (save) out
}

#' @importFrom loo loo
#' @importFrom stats sd
looic.dmc <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...)
  # Calcuate looic
{
  out <- loo(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
  if (all(out$pareto_k<.5))
    cat("All Pareto k estimates OK (k < 0.5)\n") else {
    msg <- "See PSIS-LOO description (?'loo-package') for more information"
    if (any(out$pareto_k>1))
      msg1 <- paste(sum(out$pareto_k>1)," (",
                    round(100*sum(out$pareto_k>1)/length(out$pareto_k),1),
                    "%) Pareto k estimates greater than 1\n",sep="") else
      msg1 <- paste(sum(out$pareto_k>0.5)," (",
                    round(100*sum(out$pareto_k>1)/length(out$pareto_k),1),
                    "%) Pareto k estimates estimates between 0.5 and 1\n",sep="")
     warning(msg1,msg)
  }
  if (mc_se) {
    n.chains <- dim(trial_log_likes)[2]
    loos <- vector(mode="list",length=n.chains)
    for (i in 1:n.chains) loos[[i]] <- loo(trial_log_likes[,i,])
    mc <- sd(unlist(lapply(loos,function(x){x$looic})))/sqrt(n.chains)
    print(c(p=out$p_loo,se_p=out$se_p,looic=out$looic,se_looic=out$se_looic,
      mc_se_loo=mc),...)
  } else
    print(c(p=out$p_loo,se_p=out$se_p,looic=out$looic,se_looic=out$se_looic),...)
  if (save) out
}

#' @importFrom loo compare
loocompare.dmc <- function(loo1,loo2=NULL,...)
  # Model comparision of objects produced by waic.dmc or looic.dmc
{
  if ( !is.null(loo2) ) {
    tmp <- compare(loo1,loo2)
    out <- c(waic_diff=-2*tmp[[1]],se=2*tmp[[2]])
    print(out,...)
  } else {
    if ( !is.list(loo1) )
      stop("loo1 argument must be a list if loo2 not given")
    if ( any(names(loo1[[1]])=="looic") ) icnam <- "looic" else
                                          icnam <- "waic"
    ics <- unlist(lapply(loo1,function(x){x[[icnam]]}))
    d = ics-min(ics)
    w = exp(-d/2)/sum(exp(-d/2))
    if (!is.null(names(ics))) mnams=names(ics) else
      mnams <- paste("M",1:length(d),sep="")
    print(t(matrix(c(d,w),nrow=length(d),
      dimnames=list(mnams,c("IC-min","w")))),...)
  }
}

p.fun.dmc <- function(samples,fun,hyper=FALSE,ptype=1)
  # applies fun to samples
  {
   if (!hyper) as.vector(apply(samples$theta,c(1,3),fun)) else
     as.vector(apply(attr(samples,"hyper")$phi[[ptype]],c(1,3),fun))
}

#' @importFrom matrixStats colSds
#' @export
CheckRecovery <- function(samples, p.vector, start, hyper = FALSE,
  digits = 2, verbose = TRUE)
{
  if (hyper) {

    est <- summary(samples)
    tmp <- t(data.frame(lapply(est, function(x){x[[1]][, 1]})))
    mean.est <- colMeans(tmp)
    mean.ps <- colMeans(p.vector)
    sd.est <- matrixStats::colSds(tmp)
    sd.ps <- matrixStats::colSds(p.vector)

    pnames <- colnames(ps)
    loc <- rbind(mean.est, mean.ps, mean.ps - mean.est)
    sca <- rbind(sd.est, sd.ps, sd.ps - sd.est)
    out <- rbind(loc, sca)
    rownames(out) <- c("Mean", "True", "Diff", "Sd", "True", "Diff")
    if (verbose) print(round(out, digits))
    invisible(return(out))

  } else {
    if (missing(p.vector)) stop("Please supply true values.")
    if (missing(start)) start <- samples$start

    qs <- summary.model(samples, start = start)$quantiles
    if (!is.null(p.vector) && (!all(dimnames(qs)[[1]] %in% names(p.vector))))
      stop("Names of p.vector do not match parameter names in samples")

    est  <- qs[names(p.vector), "50%"]

    op.vector <- p.vector[order(names(p.vector))]
    oest <- est[order(names(est))]
    bias <- oest- op.vector

    lo   <- qs[names(p.vector), "2.5%"]
    hi   <- qs[names(p.vector), "97.5%"]
    olo <- lo[order(names(lo))]
    ohi <- hi[order(names(hi))]

    out  <- rbind(
      'True'          = op.vector,
      '2.5% Estimate' = olo,
      '50% Estimate'  = oest,
      '97.5% Estimate'= ohi,
      'Median-True'   = bias)

    if (verbose) print(round(out, digits))
    invisible(return(out))
  }


}

# CheckRecovery.stanfit <- function(fit, p.vector, pars, digits = 2,
#                                   verbose = TRUE)
# {
#   if (missing(p.vector)) stop("Please supply true values.")
#   if (missing(pars)) pars <- fit@sim$pars_oi
#   fit_summary <- rstan::summary(fit)
#   qs <- fit_summary$summary
#   dimnames(qs)[[1]] <- c(pars, "s", "lp__")
#   est <- qs[pars, colnames(qs) %in% "50%"]
#   lo  <- qs[pars, colnames(qs) %in% "2.5%"]
#   hi  <- qs[pars, colnames(qs) %in% "97.5%"]
#
#   op.vector <- p.vector[order(names(p.vector))]
#   oest <- est[order(names(est))]
#   bias <- oest - op.vector
#
#   olo <- lo[order(names(lo))]
#   ohi <- hi[order(names(hi))]
#
#   out <- rbind('2.5% Estimate' = olo, '50% Estimate' = oest,
#                '97.5% Estimate '= ohi, 'Median-True' = bias)
#
#   if (!is.null(p.vector)) out <- rbind('True' = op.vector, out)
#   if (verbose) print(round(out, digits))
#   invisible(out)
# }

#' @export
h.check.recovery.dmc <- function(samples, ps=NULL,hyper=FALSE,verbose=TRUE,
  digits=3, ptype=1)
{
  if (hyper) {
    samples <- list(theta=attr(samples,"hyper")$phi[[ptype]])
    samples$nmc <- dim(samples$theta)[3]
    out <- check.recovery.dmc(samples,ps,digits=digits,verbose=verbose)
  } else {
    out <- samples
    for (i in 1:length(samples)) {
      out[[i]] <- check.recovery.dmc(samples[[i]], ps[i,], verbose=FALSE)
      if (i==1) av <- out[[i]] else av <- out[[i]] + av
    }
    if (verbose) print(round(av/length(samples), digits=digits))
  }
  invisible(out)
}


#' @export
gelman.diag.mpsrf <- function(mcmclist,autoburnin,transform)
  # robust version ONLY USED IN sampling.R IN run.converge.dmc
  # SHOULD BE ROLLED OUT OVER FOLLOWING FUNCITONS TO AVOID CRASHES OF
  # AUTO PROCEDURES.
{
  gd <- try(gelman.diag(mcmclist,
    autoburnin=autoburnin,transform=transform),silent=TRUE)
  if (class(gd)=="try-error") Inf else gd$mpsrf
}

#' @export
get.thin <- function(samples,hyper=FALSE) {
  es <- effectiveSize.dmc(samples,hyper=hyper)
  if (hyper) {
    print(es)
    samples <- attr(samples,"hyper")
    n <- prod(dim(samples$phi[[1]])[-2])
    cat("Thin\n")
    print(round(c(mean=n/mean(es),min=n/min(es)),0))
  } else {
    if (any(names(samples)=="theta")) {
      print(es)
      n <- prod(dim(samples$theta)[-2])
      cat("Thin\n")
      print(round(c(mean=n/mean(es),min=n/min(es)),0))
    } else {
      cat("Minimum Effective Size\n")
      print(unlist(lapply(es,min)))
      cat("Mean Effective Size\n")
      print(round(unlist(lapply(es,mean))))
      n <- unlist(lapply(samples,function(x){prod(dim(x$theta)[-2])}))
      mn <- n/unlist(lapply(es,mean))
      mi <- n/unlist(lapply(es,min))
      cat("Thin\n")
      print(round(rbind(mean=mn,min=mi),0))
    }
  }
}
