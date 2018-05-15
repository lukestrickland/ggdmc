##' Return n-cell matrix
##'
##' Contructs a matrix indicating how many responses to simulate in
##' each design cell. It also check whether the format of \code{n} and \code{ns}
##' conform.
##'
##' \code{n} can be (1) one real number for a balanced design, (2) a matrix for
##' an unbalanced design, where rows are subjects and columns are design cells.
##' If the matrix is a row vector, then all subjects have the same \code{n} in
##' each cell. If it is a column vector, then all cells have the same \code{n}.
##' Otherwise each entry specifies the \code{n} for a particular design
##' subject x design cell combination. See examples below.
##'
##' @param model a model object
##' @param ns number of subjects.
##' @param n number of trials/responses.
##'
##' @examples
##' model <- ggdmc::BuildModel(
##' p.map       = list(A="1",B="R",t0="1",mean_v="M",sd_v="M",st0="1"),
##'   match.map = list(M=list(s1=1,s2=2)),
##'   constants = c(sd_v.false=1,st0=0),
##'   factors   = list(S=c("s1","s2")),
##'   responses = c("r1","r2"),
##'   type      = "norm")
##'
##' ## Example 1
##' ggdmc::GetNsim(model, ns = 2, n = 1)
##' #      [,1] [,2]
##' # [1,]    1    1
##' # [2,]    1    1
##'
##' ## Example 2
##' n <- matrix(c(1:2),ncol=1)
##' #      [,1]
##' # [1,]    1  ## subject 1 has 1 response for each cell
##' # [2,]    2  ## subject 2 has 2 responses for each cell
##'
##' ggdmc::GetNsim(model, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    1
##' # [2,]    2    2
##'
##' ## Example 3
##' n <- matrix(c(1:2), nrow=1)
##' #      [,1] [,2]
##' # [1,]    1    2
##' ggdmc::GetNsim(model, ns = 2, n = n)
##' #     [,1] [,2]
##' # [1,]   1    2 ## subject 1 has 1 response for cell 1 and 2 responses for cell 2
##' # [2,]   1    2 ## subject 2 has 1 response for cell 1 and 2 responses for cell 2
##'
##' ## Example 4
##' n <- matrix(c(1:4), nrow=2)
##' #      [,1] [,2]
##' # [1,]    1    3
##' # [2,]    2    4
##' ggdmc::GetNsim(model, ns = 2, n = n)
##' #      [,1] [,2]
##' # [1,]    1    3 ## subject 1 has 1 response for cell 1 and 3 responses for cell 2
##' # [2,]    2    4 ## subject 2 has 2 responses for cell 1 and 4 responses for cell 2
##
##' @export
GetNsim <- function(model, n, ns) {
  if (ns <= 1) stop("Use simulate.model instead to simulate one participant.")
  if (is.vector(n) & (length(n) != 1)) {
    if (!is.matrix(n)) stop("n must be a scalar, a vector, or a matrix")
  }

  facs  <- attr(model, "factors")
  ## lapply(facs, length) return the numbers of level in each factor
  ## unlist and prod then multiply them to get the number of cell
  ncell <- prod(unlist(lapply(facs, length)))

  if (is.matrix(n)) {     ## Untested
    dim1 <- dim(n)[1]
    dim2 <- dim(n)[2]
    if (dim1 == 1) nmat <- matrix(rep(n, each = ns), ns) ## only cells differ
    if (dim2 == 1) nmat <- matrix(rep.int(n, ncell), ns) ## only subjects differ
  } else {
    nmat <- matrix(rep.int(n, ns*ncell), ns)
  }

  dim1 <- dim(nmat)[1]   ## n has been altered in ifelse
  dim2 <- dim(nmat)[2]

  if ( ns != dim1 ) stop(paste0("The n matrix must have ", ns, " rows"))
  if ( ncell != dim2 ) stop(paste0("The n matrix must have ", ncell, " columns"))
  return(nmat)
}


##' Test whether `model` object has many models
##'
##' @param model a model object
##' @param ns number of subjects.
##'
##' @export
ismanymodels <- function(model, ns = NA) {
  if (is.list(model)) {
    if (length(model) != ns) stop("n subject != n model")
    out <- TRUE
  } else {
    if (is.na(ns)) stop("ns not found")
    out <- FALSE
  }
  return(out)
}


make.hstart <- function(fsamples,
  mu.sd=1,sigma.sd=1,
  lower.mu=NULL,upper.mu=NULL,
  lower.sigma=NULL,upper.sigma=NULL,
  do.plot=FALSE,layout=c(2,5),digits=2)
  # Uses the results of fixed effects fitting to get a start prior for hierarchcial
  # By default sets prior sd to 1% of mean (tight), can specify as a vector
  # Prints prior parameters and can plot priors
  ## Change to ggdmc version
{

  mns <- lapply(fsamples,function(x){apply(x$theta,2,mean)})
  mns <- matrix(unlist(mns),nrow=length(mns[[1]]),dimnames=list(names(mns[[1]]),NULL))

  mn <- apply(mns,1,mean)
  if (length(mu.sd)==1) mu.sd <- abs(mn)*mu.sd/100
  if (length(mu.sd)!=length(mn))
    stop(paste("mu.sd must have length 1 or length of parameters:",length(mn)))
  if (is.null(lower.mu)) lower.mu <- rep(-Inf,length(mn))
  if (length(lower.mu)!=length(mn))
    stop(paste("lower.mu must have length of parameters:",length(mn)))
  if (is.null(upper.mu)) upper.mu <- rep(Inf,length(mn))
  if (length(upper.mu)!=length(mn))
    stop(paste("upper.mu must have length of parameters:",length(mn)))

  sd <- apply(mns,1,sd)
  if (length(sigma.sd)==1) sigma.sd <- sd*sigma.sd/100
  if (length(sigma.sd)!=length(sd))
    stop(paste("sigma.sd must have length 1 or length of parameters:",length(sd)))
  if (is.null(lower.sigma)) lower.sigma <- rep(-Inf,length(mn))
  if (length(lower.sigma)!=length(mn))
    stop(paste("lower.sigma must have length of parameters:",length(mn)))
  if (is.null(upper.sigma)) upper.sigma <- rep(Inf,length(mn))
  if (length(upper.sigma)!=length(mn))
    stop(paste("upper.sigma must have length of parameters:",length(mn)))

  mu.prior <- prior.p.dmc(p1=mn,p2=mu.sd,lower=lower.mu,upper=upper.mu)
  sigma.prior <- prior.p.dmc(p1=sd,p2=sigma.sd,lower=lower.sigma,upper=upper.sigma)

  cat("Mu prior\n")
  cat("Mean\n"); print(round(mn,digits))
  cat("SD\n"); print(round(mu.sd,digits))
  cat("Sigma prior\n")
  cat("Mean\n"); print(round(sd,digits))
  cat("SD\n"); print(round(sigma.sd,digits))

  if (do.plot) { plot_priors(mu.prior); plot_priors(sigma.prior) }

  list(mu=mu.prior,sigma=sigma.prior)
}


# samples=hsamples1;nmc=NA;thin=NA;max.try=100
# hstart.prior=NULL;phi1=NULL;start.from=NA

add.hyper.dmc <- function(samples, pp.prior, nmc=NA, thin=NA, max.try=100,
                          hstart.prior=NULL, phi1=NULL, p.prior=NULL)
  # Adds hyper level to fixed subject effect samples object
{

  update_constants <- function(samplesi, hyper, mci=1)
  # update summed_log_prior and log_likelihoods with no.sigma hypers
  {
    consts <- hyper$phi[[1]][,,mci][,!hyper$has.sigma,drop=FALSE]
    p.names <- dimnames(consts)[[2]]
    pp.prior <- hyper$pp.prior[[1]][p.names]
    datai <- samplesi$data
    for (i in 1:dim(samplesi$summed_log_prior)[2]) { # chains

      p.vector <- consts[i,]
      names(p.vector) <- p.names

      samplesi$summed_log_prior[mci,i] <-
        samplesi$summed_log_prior[mci,i] + summed.log.prior(p.vector, pp.prior)

      attr(attr(datai,"model"), "all.par")[p.names] <- p.vector
      samplesi$log_likelihoods[mci,i] <-
        sum(log.likelihood(samplesi$theta[i,,mci], datai))
    }

    if (any(is.na(samplesi$summed_log_prior[mci,])))
      stop("prior=NA for hyper parameter with no sigma, narror prior?")
    if (any(is.na(samplesi$log_likelihoods[mci,])))
      stop("likelihood=NA for hyper parameter with no sigma, narror prior?")
    samplesi
  }


  if (names(samples[[1]])[1]!="theta")
    stop("samples must be a list of subject sample objects")
  if (is.null(p.prior)) p.prior <- samples[[1]]$p.prior else {
    for (i in 1:length(samples)) samples[[i]]$p.prior <- p.prior
  }
  # check pp.prior
  if (!is.list(pp.prior)) stop("pp.prior must be a list")
  if ( length(pp.prior[[1]])<length(pp.prior[[2]]) )
    stop("Hyper-prior must have as many or more location than scale parameters")
  has.sigma <- names(pp.prior[[1]]) %in% names(pp.prior[[2]])
  fixed <- names(pp.prior[[1]])[ !has.sigma ]
  bad <- fixed[!(fixed %in% names(attr(attr(samples[[1]]$data,"model"),"constants")))]
  if ( length(bad)>0 )
    stop(paste("Fixed hyper parameters not constants in data level:",bad))
  pp.names <- lapply(pp.prior,names)
  has.hyper <- names(p.prior) %in% pp.names[[2]]
  if ( !all(names(p.prior)[has.hyper] %in% pp.names[[1]][has.sigma]) )
      stop("pp.prior location parameters not compatible with p.prior")
  if ( !all(names(p.prior)[has.hyper]==pp.names[[2]]) )
      stop("pp.prior scale parameters not compatible with p.prior")
  if (!is.null(hstart.prior)) { # check hstart.prior
    if (!is.list(hstart.prior)) stop("hstart.prior must be a list")
    if (!all(sort(names(pp.prior[[1]]))==sort(names(hstart.prior[[1]]))))
        stop("Incompatible location hstart.prior")
    if (!all(sort(names(pp.prior[[2]]))==sort(names(hstart.prior[[2]]))))
        stop("Incompatible scale hstart.prior")
  }
  if (is.na(nmc)) nmc <- samples[[1]]$nmc
  if (is.na(thin)) thin <- samples[[1]]$thin
  samples <- lapply(samples,function(x){
    samples.dmc(nmc, samples=x, thin=thin)
  })
  n.chains <- samples[[1]]$n.chains
  p.names <- names(pp.prior[[1]])
  n.pars <- length(p.names)
  phi <- array(NA,c(n.chains,n.pars,nmc))
  dimnames(phi)[[2]] <- p.names
  phi <- list(phi,phi)  # pairs of sampled hyper-parameters
  if (!is.null(phi1)) { # check phi1
    if (!all(dim(phi[[1]][,,1])==dim(phi1[[1]])))
      stop("phi1 location not compatible")
    if (!all(dim(phi[[1]][,,1])==dim(phi1[[1]])))
      stop("phi1 scale not compatible")
  }
  h_summed_log_prior <- array(-Inf,c(nmc,n.chains)) # hyper log-likelihoods
  h_log_likelihoods <- array(-Inf,c(nmc,n.chains))  # hyper log-likelihoods
  ntry <- 1
  repeat {
    if ( is.null(phi1) ) {
      cat("Generating hyper-start points for each chain: ")
      for( i in 1:n.chains ) {
        cat(".")
        if ( is.null(hstart.prior) ) { # sample from prior
          phi[[1]][i,,1] <- rprior(pp.prior[[1]])[,p.names]
          phi[[2]][i,has.sigma,1] <-
            rprior(pp.prior[[2]])[,p.names[has.sigma]]
        } else {                        # sample from hstart.prior
          phi[[1]][i,,1] <- rprior(hstart.prior[[1]])[,p.names]
          phi[[2]][i,has.sigma,1] <-
            rprior(hstart.prior[[2]])[,p.names[has.sigma]]
        }
      }
      cat("\n")
    } else {
      phi[[1]][,,1] <- phi1[[1]]
      phi[[2]][,,1] <- phi1[[2]]
    }
    for (i in 1:n.chains) {
      h_summed_log_prior[1,i] <-
        sum(log.prior.dmc(phi[[1]][i,,1],pp.prior[[1]])) +
        sum(log.prior.dmc(phi[[2]][i,has.sigma,1],pp.prior[[2]]))
    }

    # Fill in hyper-likelihoods
    cps <- lapply(samples,function(x){x$theta[,,1]}) # subject list: chains x pars
    for( i in 1:n.chains ) {
      h_log_likelihoods[1,i] <- sum(h.log.likelihood.dmc(p.prior=p.prior[has.hyper],
        ps=t(data.frame(lapply(cps,function(x){x[i,has.hyper]}))),
        pp=list(phi[[1]][i,has.sigma,1],phi[[2]][i,has.sigma,1])))
    }
    if ( any(!is.finite(h_log_likelihoods[1,])) ) {
      ntry <- ntry + 1
      if (ntry > max.try)
        stop(paste("For",max.try,"attempts sampling start points was valid for data but not hyper level\n
          Try a tighter pp.prior or hstart.prior"))
    } else break
  }
  cat("\n")

  # Update subject level priors to fit with hyper ...
  for (j in 1:length(samples)) for ( i in 1:n.chains) {
    data <- samples[[j]]$data
    samples[[j]]$p.prior <- assign.pp(
      list(phi[[1]][i, has.sigma,1],
      phi[[2]][i, has.sigma,1]),
      p.prior = samples[[i]]$p.prior[has.hyper])

    samples[[j]]$summed_log_prior[1,i]  <- summed.log.prior(
      samples[[j]]$theta[i,,1], samples[[j]]$p.prior)
  }

  hyper <- list(phi=phi,h_log_likelihoods=h_log_likelihoods,
    h_summed_log_prior=h_summed_log_prior,pp.prior=pp.prior,start=1,thin=thin,
    n.pars=n.pars,p.names=p.names,rp=samples[[1]]$rp,nmc=nmc,n.chains=n.chains,
    has.sigma=has.sigma,has.hyper=has.hyper)
  if (any(!has.sigma)) # Update data level for !has.sigma hyper parameters
    samples <- lapply(samples,update_constants,hyper=hyper)
  attr(samples,"hyper") <- hyper
  samples
}


#' @importFrom stats runif
h.migrate <- function(use.phi,use.logprior,use.loglike,p.prior,ps,rp,pp.prior,
                      has.sigma,has.hyper,is.constant)
  # DEMCMC migrate set, all chains, hyper level
{

  # Which pars are not constants?
  de <- list(!unlist(is.constant[[1]]), !unlist(is.constant[[2]]))

  n.chains <- dim(use.phi[[1]])[1]
  pars <- dim(use.phi[[1]])[2]
  lnum1 <- sample(c(1:n.chains),1)  # determine how many groups to work with
  lnum2 <- sort(sample(c(1:n.chains),lnum1,replace=F))  # get groups
  phiset <- list( # initialize
    matrix(NA,lnum1,pars,dimnames=list(NULL,dimnames(use.phi[[1]])[[2]])),
    matrix(NA,lnum1,pars,dimnames=list(NULL,dimnames(use.phi[[1]])[[2]]))
  )

  propset.logprior <- propset.loglike <- numeric(lnum1)
  currentset <- propset <- propw <- currw <- numeric(lnum1)

  index <- numeric(lnum1)
  for (i in 1:lnum1) {
    index[i] <- sample(1:n.chains,1,replace=F)
    # create a set of these particles to swap
    phiset[[1]][i,] <- use.phi[[1]][lnum2[i],]
    phiset[[2]][i,] <- use.phi[[2]][lnum2[i],]

    # perturb non-constant parameters
    phiset[[1]][i,de[[1]]] <- phiset[[1]][i,de[[1]]] + runif(1,-rp,rp)
    phiset[[2]][i,de[[2]]] <- phiset[[2]][i,de[[2]]] + runif(1,-rp,rp)

    propset.logprior[i] <- h.summed.log.prior (
      pp=list(phiset[[1]][i, has.sigma], phiset[[2]][i,has.sigma]), pp.prior=pp.prior)
    propset.loglike[i]  <- sum(h.log.likelihood.dmc(ps=ps[i,,has.hyper],
      pp=list(phiset[[1]][i,has.sigma],phiset[[2]][i,has.sigma]),p.prior=p.prior[has.hyper]))

    propset[i] <- propset.logprior[i] + propset.loglike[i]
    propw[i] <- propset[i]
    currentset[i] <- use.logprior[lnum2[i]] + use.loglike[lnum2[i]]

  }
  currw <- currentset

  mh <- exp(propw[lnum1] - currw[1])
  if ( !is.na(mh) && (runif(1) < mh) ) {
    use.phi[[1]][lnum2[1],] <- phiset[[1]][lnum1,]	# swap the 1st with last
    use.phi[[2]][lnum2[1],] <- phiset[[2]][lnum1,]  #  (creating a circle)
    use.logprior[lnum2[1]] <- propset.logprior[lnum1]
    use.loglike[lnum2[1]] <- propset.loglike[lnum1]
  }
  if ( lnum1!=1 ) {										# make sure we are not done yet
    for(i in 1:(lnum1-1)){
      mh <- exp(propw[i] - currw[i+1])
      if ( !is.na(mh) && (runif(1) < mh) ) {
        use.phi[[1]][lnum2[i+1],] <- phiset[[1]][i,]
        use.phi[[2]][lnum2[i+1],] <- phiset[[2]][i,]
        use.logprior[lnum2[i+1]] <- propset.logprior[i]
        use.loglike[lnum2[i+1]] <- propset.loglike[i]
      }
    }
  }
  cbind(use.logprior,use.loglike,use.phi[[1]],use.phi[[2]])
}

blocked.h.crossover <- function(k,blocks,n.pars,use.phi,use.logprior,use.loglike,
  p.prior,ps,pp.prior,rp,has.sigma,has.hyper,is.constant,
  force = FALSE, gamma.mult=2.38,h.gamma.mult=2.38,random.theta=TRUE)
  # hyper level crossover for hypers with a sigma, as a series of blocks
{

  h.crossover <- function(k, pars, use.phi, use.logprior, use.loglike,
    p.prior, ps, is.constant, pp.prior, rp, has.sigma, has.hyper,
    h.gamma.mult=2.38, force = FALSE, random.theta=TRUE)
    # DEMCMC crossover update of one chain, phi level
  {
    if (random.theta) k.theta <- sample(1:dim(ps)[1],1) else k.theta <- k

    # Which pars are not constants?
    de <- list(!unlist(is.constant[[1]][names(is.constant[[1]])[pars]]),
      !unlist(is.constant[[2]][names(is.constant[[1]])[pars]]))
    if ( all(!unlist(de)) ) # all are constants
      return(c(use.logprior[k], use.loglike[k], use.phi[[1]][k,], use.phi[[2]][k,]))

    # step size
    if ( is.na(h.gamma.mult) ) hgamma <- runif(1,0.5,1) else
      hgamma <-  h.gamma.mult/sqrt(2*length(pars)*2) # extra *2 as p1 and p2

    # Update use.loglike for new ps
    phi <- list(use.phi[[1]][k,],use.phi[[2]][k,])
    use.loglike[k] <- sum(h.log.likelihood.dmc(
      ps      = ps[k.theta,,], # has.hyper
      pp      = list(phi[[1]][has.sigma], phi[[2]][has.sigma]),
      p.prior = p.prior[has.hyper]))

    # Calcualte use.post
    use.post <-  use.logprior[k] + use.loglike[k]

    # DE step
    # pick two other chains
    chains <- c(1:dim(use.phi[[1]])[1])[-k]
    index <- sample(chains, 2)

    # update mu
    phi[[1]][pars[de[[1]]]] <- use.phi[[1]][k,pars[de[[1]]]] + runif(1,-rp,rp) +
      hgamma*(use.phi[[1]][index[1],pars[de[[1]]]]-use.phi[[1]][index[2],pars[de[[1]]]])

    # update sigma
    phi[[2]][pars[de[[2]]]] <- use.phi[[2]][k,pars[de[[2]]]] + runif(1,-rp,rp) +
      hgamma*(use.phi[[2]][index[1],pars[de[[2]]]]-use.phi[[2]][index[2],pars[de[[2]]]])

    # Get new post
    logprior <- h.summed.log.prior(pp.prior=pp.prior,
      pp = list(phi[[1]][has.sigma], phi[[2]][has.sigma])) # only for has.sigma

    loglike <- sum(h.log.likelihood.dmc(
      ps=ps[k.theta,,], #has.hyper
      p.prior=p.prior[has.hyper],
      pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma])))
    post <- logprior + loglike

    # Metropolis step
    epup <- exp(post - use.post)
    if ( force || ( !is.na(epup)  && (runif(1) < epup) ) ) {
      use.phi[[1]][k, pars] <- phi[[1]][pars]
      use.phi[[2]][k, pars] <- phi[[2]][pars]
      use.logprior[k] <- logprior
      use.loglike[k] <- loglike
    }

    c(use.logprior[k],use.loglike[k],use.phi[[1]][k,],use.phi[[2]][k,])
  }


  temp <- h.crossover(k, pars=blocks[[1]],use.phi=use.phi,force=force,
    use.logprior=use.logprior,use.loglike=use.loglike,
    p.prior=p.prior,ps=ps,pp.prior=pp.prior,rp=rp,
    h.gamma.mult=h.gamma.mult,is.constant=is.constant,
    has.sigma=has.sigma,has.hyper=has.hyper,
    random.theta=random.theta)

  if ( length(blocks)>1 ) for ( b in 2:length(blocks) ) {
    use.logprior[k]  <- temp[1]
    use.loglike[k]   <- temp[2]
    use.phi[[1]][k,] <- temp[3:(n.pars+2)]
    use.phi[[2]][k,] <- temp[(n.pars+3):(2*n.pars+2)]
    temp <- h.crossover(k, pars=blocks[[b]],use.phi=use.phi,force=force,
      use.logprior=use.logprior,use.loglike=use.loglike,
      p.prior=p.prior,ps=ps,pp.prior=pp.prior,rp=rp,
      h.gamma.mult=h.gamma.mult,is.constant=is.constant,
      has.sigma=has.sigma,has.hyper=has.hyper,
      random.theta=random.theta)
  }
  temp
}


#' @rdname h.run.dmc
#' @export
checkforce <- function(force = FALSE, nsamp = NA) {
  if( is.na(nsamp) ) stop("nsamp not found")
  if (is.numeric(force)) force <- c(rep(FALSE, force - 1), TRUE)
  force <- c(TRUE, rep(force, length.out = nsamp - 1))

  if (!is.logical(force) || (length(force) != nsamp))
    stop(paste("force argument must be a logical vector of length", nsamp-1))
  return(force)
}

#' @rdname h.run.dmc
#' @export
checkblocks <- function(blocks, hyper) {
  if ( any(is.na(blocks)) ) {
    blocks <- as.list(1:hyper$n.pars)
  } else { # check
    if (any(unlist(lapply(blocks, function(x) {
      length(x)==1 || all(hyper$has.sigma[x][1] == hyper$has.sigma[x][-1])
    }))))
      stop("Cant mix hyper-paramaters with and without sigma in a block")
  }
  return(blocks)
}


# samples,cores=1, report=samples[[1]]$nmc, blocks=NA,
# p.migrate=0,h.p.migrate=0,force.hyper=NA,force.data=NA,
# gamma.mult=2.38,h.gamma.mult=NA,random.phi=TRUE,random.theta=TRUE

#' Fit Fixed-Effect or Random-Effect Model
#'
#' \code{h.run.dmc} calls \code{run_data} in the case of fixed-effect model, or
#' \code{run_hyper} in the case of random-effect model. At the moment, it fits
#' either fixed-effect or random-effect model by looking for the
#' \code{hyper} attribute in a DMC sample/object.
#'
#' \code{hrun.dmc} use AH's h.run.dmc, including all his blocked.h.crossover,
#' h.crossover, migrate.h, crossover.h, but use YSL's dprior and density_rd.
#'
#' @param samples a DMC samples, usually generated by calling
#' \code{h.samples.dmc}.
#' @param report the progress report. Default is every 100 iterations report
#' once
#' @param cores a switch for parallel computation.  In the case of fixed-effect
#' model, \code{h.run.dmc} calls \pkg{snow} (via \pkg{snowfall}) or
#' \pkg{parallel}. When setting for example \code{cores=4}, \code{h.run.dmc}
#' requests 4 parallel R processes. In the case of random-effect model,
#' \code{h.run.dmc} calls OpenMP library, which will automatically decide the
#' number of cores to use, so setting core number higher than 1 has the same
#' effect.  Use OpenMP only when the data set is large (e.g., more than 1000
#' trial per condition).
#' @param blocks a list for DE-MCMC parameter for blocking. Unknown way of
#' implmentation.
#' @param p.migrate the probability of applying migration sampler. Set it
#' greater than will trigger the migration sampler. For example
#' \code{p.migrate=0.05} will use migration sampler about 5\% chance.
#' @param h.p.migrate similar to p.migrate for hyper level
#' @param force.hyper This argument has no funciton. It is here only for
#' compatibility (with DMC) reason.
#' @param force.data This argument has no funciton. It is here only for
#' compatibility (with DMC) reason.
#' @param gamma.mult a DE-MCMC tuning parameter, affecting the size of jump.
#' Default is 2.38.
#' @param h.gamma.mult Similar with gamm.mult.  This argument has no funciton.
#' It is here only for compatibility reason.
#' @param setting a list to store setting arguments, such as p.migrate,
#' gamma.mult, etc, sending to C++ function. The user should not use this
#' argument
#' @param verbose a swtich to see debugging information
#' @keywords h.run.dmc
#' @importFrom parallel mclapply makeCluster parLapply stopCluster
#' @examples
#' m1 <- model.dmc(
#' p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#' match.map = list(M=list(s1="r1", s2="r2")),
#' factors   = list(S=c("s1", "s2")),
#' constants = c(st0=0, d=0),
#' responses = c("r1","r2"),
#' type      = "rd")
#'
#' ## Population distribution
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' dat <- h.simulate.dmc(m1, nsim=30, p.prior=pop.prior, ns=8)
#' mdi <- BindDataModel(dat, m1)
#' p.prior  <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## Fixed-effect model
#' samplesInit <- h.samples.dmc(nmc=20, p.prior=p.prior, data=mdi, thin=1)
#' samples0    <- h.run.dmc(samples=samplesInit, report=10)
#'
#' ## Make a hyper-prior list
#' mu.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", 6),
#'   p1    = c(a=1, v=1, z=1, sz=1, sv=1, t0=1),
#'   p2    = c(1,1,1,1,1,1),
#'   upper = c(2,2,2,2,2,2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' ## Random-effect model
#' \dontrun{
#' hs0 <- h.samples.dmc(3e2, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi, thin=1)
#' hs0 <- h.run.dmc(samples=hs0, report=1e2, p.migrate=.05,
#'   h.p.migrate=.05)
#'
#' hs1 <- ggdmc:::h.samples.dmc(3e3, p.prior, pp.prior = pp.prior,
#'          samples = hs0, thin=1)
#' hs1 <- ggdmc:::h.run.dmc(hs1, 1, 5e2)
#' }
#' @export
h.run.dmc <- function(samples, ncore = 1, report=100, p.migrate=0,
  force = NA, gamma.mult = 2.38, debug=FALSE)
{
  os <- get_os()

  ## Set-up default/user DE-MCMC parameters
  if( is.null(setting) ) {
    setting_in <- list(p.migrate=p.migrate, h.p.migrate=h.p.migrate,
      gamma.mult=gamma.mult, h.gamma.mult=2.38, report=report, ncore = ncore)
    if(debug) {
      cat("DE-MCMC tuning paramters: p.migrate =", p.migrate,
        ", h.p.migrate = 0, gamma.mult =", gamma.mult,
        ", h.gamma.mult =", h.p.migrate)
      cat(". ", ncore,  " CPU(s) is/are launched. I'll report progress every",
        report, "step(s).\n")
    }
  } else {
    isSet <- names(setting) %in% c("p.migrate","h.p.migrate","gamma.mult",
      "h.gamma.mult", "cores", "report")
    setting$p.migrate    <- ifelse(isSet[1], setting$p.migrate,    "0")
    setting$h.p.migrate  <- ifelse(isSet[2], setting$h.p.migrate,  "0")
    setting$gamma.mult   <- ifelse(isSet[3], setting$gamma.mult,   "2.38")
    setting$h.gamma.mult <- ifelse(isSet[4], setting$h.gamma.mult, "2.38")
    setting$ncore        <- ifelse(isSet[5], setting$ncore, "1")
    setting$report       <- ifelse(isSet[6], setting$report, "100")
    if(debug) {
      cat("DE-MCMC tuning parameters: p.migrate = ", setting$p.migrate)
      cat(" h.p.migrate =", setting$h.p.migrate)
      cat(" gamma.mult =", setting$gamma.mult)
      cat(" h.gamma.mult =", setting$h.gamma.mult)
      cat(". run will use ", setting$ncore, "CPU(s) and report every ")
      cat(setting$report, " step.\n");
    }
    setting_in <- setting
  }

  ## If hyper attribute is not set, it must be fixed-effect model
  if ( is.null( attr(samples, ("hyper"))) ) {
    ## If the 1st name is theta, assume we get only one participant
    ## Otherwise, mulitple participants
    if (debug) { cat("Run fixed-effect: ") }
    sNames <- names(samples)      ## samples or subjects' names
    if (sNames[1] == "theta" ) {  ## fixed-effect for one subject
      if (ncore > 1) {
        if (debug) { cat("Run with OMP. ncore = 1 turns off parallel\n") }
        out <- run_data_parallel(samples, setting_in)
      } else {
        out <- run_data(samples, setting_in)
      }

      dimnames(out$theta) <- list(NULL, out$p.names, NULL)
    } else {                    ## fixed-effect for multiple subjects

      nsamp <- 1 + (samples[[1]]$nmc - samples[[1]]$start) * samples[[1]]$thin
      force <- rep.int(0, nsamp)
      if (is.numeric(force.data)) {
        for(i in 1:nsamp) { if (i %% force.data == 0) {force[i] <- 1} }
      }

      if (os == "windows" & ncore > 1) {
        cl  <- makeCluster(ncore)
        out <- parLapply(cl, samples, run_data, setting_in, force)
        stopCluster(cl)
      } else if (ncore > 1) {
        ## system.time(s0 <- ggdmc::run(s0, 50, 1, p.migrate = .05))

        out <- mclapply(samples, run, report, ncore, p.migrate, gamma.mult,
          force.data, mc.cores=ncore)
      } else {
        if (debug) cat("Use lapply\n")
        out <- lapply(samples, run_data, setting_in, force)
      }
      for(i in 1:length(out))
        dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
    }
    out <- dmc.list(out)
    ## random-effect
  } else if (!is.null(attr(samples, "hyper"))) {
    if (debug) cat("Run random-effect fit.\n")

    hyper <- attr(samples, "hyper")
    nsamp <- 1 + (hyper$nmc - hyper$start) * hyper$thin
    force <- rep.int(0, nsamp)
    if (is.numeric(force.data)) {
      for(i in 1:nsamp) { if (i %% force.data == 0) {force[i] <- 1} }
    }

    if (ncore > 1) { ## run OMP
      out <- run_hyper_parallel(samples, setting_in)
    } else {
      out <- run_hyper(samples, setting_in, force)
    }

    for(i in 1:length(out)) {
      dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
    }

    hyper0 <- attr(out, "hyper")
    for(j in 1:length(hyper0$phi)) {
      dimnames(hyper0$phi[[j]]) <- list(NULL, out[[1]]$p.names, NULL)
    }

    attr(out, "hyper") <- hyper0
    names(out) <- 1:length(out)
    out <- hyper(out)

  } else {
    stop("Unexpected condition occurred")
  }

  cat("\n")
  return(out)
}

#' @export
h.runR.dmc <- function(samples,cores=1,report=samples[[1]]$nmc,blocks=NA,
  p.migrate=0,h.p.migrate=0,force.hyper=NA,force.data=NA,
  gamma.mult=2.38,h.gamma.mult=NA,random.phi=TRUE,random.theta=TRUE)
  # Fixed effects: run sampling spreading subjects over cores. For cores=1 report
  #   is made for each subject, so default is just reprot when nmc done for each
  # Hierarchical: chains spread over cores at both subject and hyper level.
  #   Hyper sampling is blocked as specified in blocks (hyper-parameter pairs done
  #   together), by default one pair at a time.
{

  # mci=i;thin=hyper$thin;phi=use.phi;has.sigma=hyper$has.sigma;has.hyper=hyper$has.hyper
  run.chains <- function(samplesi,mci,thin,force=FALSE,
    phi,has.sigma,has.hyper,pp.prior,
    p.migrate=p.migrate,gamma.mult=2.38,random.phi=TRUE)
    # Data level for hierarchical
  {

    p.priors <- vector(mode="list",length=samplesi$n.chains)
    consts <- phi[[1]][,!has.sigma,drop=FALSE]
    p.names <- dimnames(consts)[[2]]
    pp.priors <- pp.prior[[1]][p.names]

    for (k in 1:samplesi$n.chains) {
      if (random.phi) k.phi <- sample(1:samplesi$n.chains,1) else k.phi <- k
      p.priors[[k]] <- assign.pp(list(phi[[1]][k.phi,has.sigma],
        phi[[2]][k.phi,has.sigma]),p.prior=samplesi$p.prior[has.hyper])
      if (any(!has.hyper)) for (i in names(samplesi$p.prior)[!has.hyper])
        p.priors[[k]][[i]] <- samplesi$p.prior[[i]]
      # Update use.logprior based on new priors
      attr(samplesi,"use")$logprior[k] <- summed.log.prior(
        p.vector=attr(samplesi,"use")$theta[k,],p.prior=p.priors[[k]],
        p.names=names(attr(attributes(samplesi$data)$model,"p.vector"))) +
        summed.log.prior(consts[k,],pp.priors,p.names)
    }

    if ( runif(1) < p.migrate ) {        # Do migration
      temp <- migrate.h(use.theta      = attr(samplesi,"use")$theta,
        use.logprior       = attr(samplesi,"use")$logprior,
        use.loglike        = attr(samplesi,"use")$loglike,
        p.priors           = p.priors,
        data               = samplesi$data,
        rp                 = samplesi$rp,
        consts             = consts,
        pp.priors          = pp.priors)
    } else {                             # Do crossover
      temp <- t(sapply(1:samplesi$n.chains,crossover.h,
        pars=1:samplesi$n.pars,
        use.theta    = attr(samplesi,"use")$theta,
        use.logprior = attr(samplesi,"use")$logprior,
        use.loglike  = attr(samplesi,"use")$loglike,
        p.priors     = p.priors,
        data         = samplesi$data,
        rp           = samplesi$rp,
        gamma.mult   = gamma.mult,
        consts       = consts,
        force        = force,
        pp.priors    = pp.priors))
    }

    # Harvest results
    attr(samplesi,"use")$logprior <- temp[,1]
    attr(samplesi,"use")$loglike  <- temp[,2]
    attr(samplesi,"use")$theta   <- temp[,-c(1,2)]
    attr(samplesi,"use")$p.priors <- p.priors

    # store samples
    if ( mci %% thin == 0 ) {
      attr(samplesi,"use")$store_i <- attr(samplesi,"use")$store_i + 1
      samplesi$summed_log_prior[attr(samplesi,"use")$store_i,] <- temp[,1]
      samplesi$log_likelihoods[attr(samplesi,"use")$store_i,]  <- temp[,2]
      samplesi$theta[,,attr(samplesi,"use")$store_i]     <- temp[,-c(1,2)]
    }

    samplesi
  }

  no.sig.crossover <- function(k,samples,pars,rp,hgamma,pp.prior,force=FALSE,
    use.phi,use.logpriors,use.loglikes)
    # one chain hyper crossover for no sigma parameters
  {
    # DE step
    # pick two other chains
    index <- sample(c(1:dim(use.phi[[1]])[1])[-k],2,replace=F)
    use.consts <- use.phi[[1]][k,pars]
    consts <- use.phi[[1]][k,pars] + runif(1,-rp,rp) +
      hgamma*(use.phi[[1]][index[1],pars]-use.phi[[1]][index[2],pars])
    p.names <- names(consts)
    pp.priord <- hyper$pp.prior[[1]][p.names]
    loglikes <- logpriors <- numeric(length(samples))
    for ( s in 1:length(samples) ) {
      use.logpriors[k,s] <- summed.log.prior(
        p.vector=attr(samples[[s]],"use")$theta[k,],
        p.prior=attr(samples[[s]],"use")$p.priors[[k]],
        p.names=names(attr(attributes(samples[[s]]$data)$model,"p.vector"))
      ) + summed.log.prior(use.consts,pp.priord,p.names)
      logpriors[s] <- summed.log.prior(
        p.vector=attr(samples[[s]],"use")$theta[k,],
        p.prior=attr(samples[[s]],"use")$p.priors[[k]],
        p.names=names(attr(attributes(samples[[s]]$data)$model,"p.vector"))
      ) + summed.log.prior(consts,pp.priord,p.names)
      datai <- samples[[s]]$data
      attr(attr(datai,"model"),"all.par")[p.names] <- use.consts
      use.loglikes[k,s] <- sum(log.likelihood(
        attr(samples[[s]],"use")$theta[k,],datai))
      attr(attr(datai,"model"),"all.par")[p.names] <- consts
      loglikes[s] <- sum(log.likelihood(
        attr(samples[[s]],"use")$theta[k,],datai))
    }
    post <- sum(logpriors + loglikes)
    use.post <- sum(use.logpriors[k,] + use.loglikes[k,])
    # Metropolis step
    epup <- exp(post-use.post)
    if ( force || (!is.na(epup)  && (runif(1) < epup)) ) {
      use.phi[[1]][k,pars] <- consts
      use.logpriors[k,] <- logpriors
      use.loglikes[k,] <- loglikes
    }
    c(use.logpriors[k,],use.loglikes[k,],use.phi[[1]][k,pars])
  }

  os <- get.os()
  hyper <- attr(samples,"hyper")
  if ( !is.null(hyper) ) # hierarchical sampling
  {

    # Setup
    if ( any(is.na(blocks)) ) blocks <- as.list(1:hyper$n.pars) else { # check
      if (any(unlist(lapply(blocks,function(x){
        length(x)==1 || all(hyper$has.sigma[x][1]==hyper$has.sigma[x][-1])}))))
        stop("Cant mix hyper-paramaters with and without sigma in a block")
    }
    sigma.block <- unlist(lapply(blocks,function(x){all(hyper$has.sigma[x])}))

    if (is.null(hyper$thin)) hyper$thin <- 1
    nsamp <- 1+(hyper$nmc-hyper$start)*hyper$thin

    if (any(is.na(force.data))) force.data <- rep(FALSE,nsamp-1)
    if (any(is.na(force.hyper))) force.hyper <- rep(FALSE,nsamp-1)
    if (!is.logical(force.hyper) || (length(force.hyper)!=(nsamp-1)))
      stop(paste("force.hyper argument must be a logical vector of length",nsamp-1))
    force.hyper <- c(TRUE,force.hyper)
    if (!is.logical(force.data) || (length(force.data)!=(nsamp-1)))
      stop(paste("force.data argument must be a logical vector of length",nsamp-1))
    force.data <- c(TRUE,force.data)

    # Setup "use" for hyper level
    ps <- aperm(array(dim=c(dim(samples[[1]]$theta)[-3],length(samples)),
      data=unlist(lapply(samples,function(x){x$theta[,,hyper$start]}),use.names=FALSE),
      dimnames=list(NULL,samples[[1]]$p.names,NULL)),c(1,3,2))
    use.phi <- list(hyper$phi[[1]][,,hyper$start],hyper$phi[[2]][,,hyper$start])
    use.logprior <- hyper$h_summed_log_prior[hyper$start,]
    use.loglike <- hyper$h_log_likelihoods[hyper$start,]
    # Setup hyper level
    store_i <- hyper$start
    pp.prior <- hyper$pp.prior
    is.constant <- lapply(pp.prior,function(x){
      lapply(x,function(y){attr(y,"dist")=="constant"})})

    # Setup for data level
    p.priors <- vector(mode="list",length=samples[[1]]$n.chains)
    for (k in 1:samples[[1]]$n.chains) {
      p.priors[[k]] <- assign.pp(
        list(use.phi[[1]][k,hyper$has.sigma],
          use.phi[[2]][k,hyper$has.sigma]),
        p.prior=samples[[1]]$p.prior[hyper$has.hyper])
    }
    samples <- lapply(samples,function(x){attr(x,"use") <- list(
      theta=x$theta[,,x$start],
      logprior=x$summed_log_prior[x$start,],
      loglike=x$log_likelihoods[x$start,],
      store_i = x$start,
      p.priors = p.priors
    ); x })
    p.prior=samples[[1]]$p.prior

    if (os == "windows" & cores > 1) {
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK")
      sfClusterSetupRNG()
      sfLibrary(msm)
      sfLibrary(rtdists)
      sfLibrary(pracma)
      sfLibrary(statmod)
      sfExportAll()
    } else if (cores > 1) {
      require(parallel, quietly=TRUE)
    }

    for (i in 2:nsamp)
    {

      # Update data and hyper level for group parameters (hypers with no sigma)
      if ( any(!sigma.block) ) {
        no.sig.blocks=blocks[!sigma.block]
        for ( j in 1:length(no.sig.blocks) ) {
          pars <- no.sig.blocks[[j]]
          # step size
          if ( is.na(h.gamma.mult) ) hgamma <- runif(1,0.5,1) else
            hgamma <-  h.gamma.mult/sqrt(2*length(pars))
          # subject list of logprior and loglike for each chain
          use.logpriors <- matrix(unlist(lapply(samples,function(x){
            attr(x,"use")$logprior})),nrow=samples[[1]]$n.chains)
          use.loglikes <-  matrix(unlist(lapply(samples,function(x){
            attr(x,"use")$loglike})),nrow=samples[[1]]$n.chains)

          if (cores > 1 & os == "windows") {
            temp <- sfLapply(1:hyper$n.chains,no.sig.crossover,samples=samples,
              pars=pars,hgamma=hgamma,pp.prior=pp.prior,use.phi=use.phi,
              force=force.hyper[i],rp=hyper$rp,
              use.logpriors=use.logpriors,use.loglikes=use.loglikes)
          } else if (cores > 1) {
            temp <- t(mclapply(1:hyper$n.chains, no.sig.crossover,
              samples=samples,pars=pars,hgamma=hgamma,pp.prior=pp.prior,
              use.phi=use.phi,force=force.hyper[i],rp=hyper$rp,
              use.logpriors=use.logpriors,
              use.loglikes=use.loglikes,mc.cores=cores))
          } else {
            temp <- lapply(1:hyper$n.chains,no.sig.crossover,samples=samples,
              pars=pars,force=force.hyper[i],hgamma=hgamma,pp.prior=pp.prior,
              use.phi=use.phi,rp=hyper$rp,
              use.logpriors=use.logpriors,use.loglikes=use.loglikes)
          }

          temp <- t(matrix(unlist(temp),ncol=hyper$n.chains))
          for (s in 1:length(samples)) {
            attr(samples[[s]],"use")$logprior <- temp[,s]
            attr(samples[[s]],"use")$loglike <- temp[,s+length(samples)]
          }
          use.phi[[1]][,pars] <- temp[,-c(1:(2*length(samples))),drop=FALSE]
        }
      }

      # Update hyper level with sigma
      if ( runif(1)<h.p.migrate ) {
        temp <- h.migrate(use.phi=use.phi,
          use.logprior=use.logprior,
          use.loglike=use.loglike,
          p.prior=p.prior,is.constant=is.constant,
          ps=ps,pp.prior=pp.prior,
          rp=hyper$rp,has.hyper=hyper$has.hyper,
          has.sigma=hyper$has.sigma)
      } else if ( cores>1 & os == "windows") {
        temp <- t(sfLapply(1:hyper$n.chains,blocked.h.crossover,
          blocks=blocks[sigma.block],n.pars=hyper$n.pars,
          use.phi=use.phi,force=force.hyper[i],
          use.logprior=use.logprior,
          use.loglike=use.loglike,random.theta=random.theta,
          p.prior=p.prior,ps=ps,rp=hyper$rp,
          pp.prior=pp.prior,is.constant=is.constant,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          h.gamma.mult=h.gamma.mult))
        temp=matrix(unlist(temp,use.names=FALSE),nrow=hyper$n.chains,byrow=TRUE)
      } else if ( cores>1 ) {
        temp <- t(mclapply(1:hyper$n.chains,blocked.h.crossover,
          blocks=blocks[sigma.block],n.pars=hyper$n.pars,
          use.phi=use.phi,force=force.hyper[i],
          use.logprior=use.logprior,
          use.loglike=use.loglike,random.theta=random.theta,
          p.prior=p.prior,ps=ps,rp=hyper$rp,
          pp.prior=pp.prior,is.constant=is.constant,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          h.gamma.mult=h.gamma.mult, mc.cores=cores))
        temp=matrix(unlist(temp,use.names=FALSE),nrow=hyper$n.chains,byrow=TRUE)
      } else {
        temp <- t(sapply(1:hyper$n.chains,blocked.h.crossover,
          blocks=blocks[sigma.block],n.pars=hyper$n.pars,
          use.phi=use.phi,force=force.hyper[i],
          use.logprior=use.logprior,
          use.loglike=use.loglike,random.theta=random.theta,
          p.prior=p.prior,ps=ps,rp=hyper$rp,
          pp.prior=pp.prior,is.constant=is.constant,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          h.gamma.mult=h.gamma.mult))
      }

      use.logprior <- temp[,1]
      use.loglike <- temp[,2]
      use.phi[[1]][,] <- temp[,-c(1,2)][,1:hyper$n.pars]
      use.phi[[2]][,] <- temp[,-c(1,2)][,(hyper$n.pars+1):(2*hyper$n.pars)]

      # Update data level for all but !has.sigma parameters
      if (cores>1 & os=="windows") {
        samples <- sfLapply(samples,run.chains,mci=i,p.migrate=p.migrate,
          thin=hyper$thin,gamma.mult=gamma.mult,phi=use.phi,pp.prior=pp.prior,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          force=force.data[i],random.phi=random.phi)
      } else if (cores>1) {
        samples <- mclapply(samples,run.chains,mci=i,
          p.migrate=p.migrate,thin=hyper$thin,gamma.mult=gamma.mult,
          phi=use.phi,pp.prior=pp.prior,has.sigma=hyper$has.sigma,
          has.hyper=hyper$has.hyper,force=force.data[i],
          random.phi=random.phi,mc.cores=cores)
      } else {
        samples <- lapply(samples,run.chains,mci=i,p.migrate=p.migrate,
          thin=hyper$thin,gamma.mult=gamma.mult,phi=use.phi,pp.prior=pp.prior,
          has.sigma=hyper$has.sigma,has.hyper=hyper$has.hyper,
          force=force.data[i],random.phi=random.phi)
      }

      ps <- aperm(array(dim=c(dim(samples[[1]]$theta)[-3],length(samples)),
        data=unlist(lapply(samples,function(x){attr(x,"use")$theta}),use.names=FALSE),
        dimnames=list(NULL,samples[[1]]$p.names,NULL)),c(1,3,2))

      if ( i %% hyper$thin == 0 ) { # store samples
        store_i <- store_i + 1
        if (store_i %% report == 0) cat(store_i," ")
        hyper$h_summed_log_prior[store_i,] <- use.logprior
        hyper$h_log_likelihoods[store_i,] <- use.loglike
        hyper$phi[[1]][,,store_i] <- use.phi[[1]]
        hyper$phi[[2]][,,store_i] <- use.phi[[2]]
      }


    }
    cat("\n")
    if (cores>1 & os=="windows") { sfStop() }
    attr(samples,"hyper") <- hyper
  } else { # fixed effect sampling
    s.names <- names(samples)
    if (any(is.na(force.data))) force.data <- FALSE
    if ( cores>1 & length(samples)>1 & os=="windows") {
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK")
      sfClusterSetupRNG()
      sfLibrary(msm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfExportAll()
      samples <- sfLapply(samples, run.dmc,p.migrate=p.migrate,report=report,
        force=force.data)
      sfStop()
    } else if (cores>1) {
      require(parallel, quietly=TRUE)
      samples <- mclapply(samples, run.dmc, p.migrate=p.migrate,
        report=report,force=force.data,mc.cores=cores)
    } else {
      samples <- lapply(samples,run.dmc,p.migrate=p.migrate,report=report,
        force=force.data)
    }
    names(samples) <- s.names
  }
  samples
}


#' @rdname h.run.dmc
#' @export
h.run.unstuck.dmc <- function(samples,nmc=NA,report=10,cores=1,
  cut=10,nbad=0,max.try=20,verbose=TRUE,
  p.migrate=0,h.p.migrate=0,end.no.migrate=FALSE,
  h.gamma.mult=NA,gamma.mult=2.38,slaveOutfile=NULL)
  # Like run.unstuck but applied to list of subjects. If hyper present then
  # unsticks that as well as subjects.
{
  if ( !is.null(samples$theta) )
    stop("For a single subject use run.unstuck.dmc")
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc
  if ( any(is.na(samples[[1]]$theta[,,2])) ) {
    cat("Getting initial set of samples\n")
    samples <- h.run.dmc(samples=samples,report=report,cores=cores,
      gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult,
      p.migrate=p.migrate,h.p.migrate=h.p.migrate)
  }
  if ( any(names(attributes(samples))=="hyper") ) {
    try.num <- 0
    repeat {
      stucks <- lapply(samples, pickStuck, cut=cut)
      ok <- c(hyper = length(pickStuck(samples, hyper=TRUE, cut=cut)),
        unlist(lapply(stucks, function(x){length(x)})))

      if (verbose) {
        cat(paste("\nAfter try",try.num,"number of stuck chains:\n"))
        print(sort(ok[ok>0],decreasing=TRUE))
      }

      if (try.num > max.try | all(ok<=nbad)) break
      samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc),
        cores=cores,report=report,p.migrate=p.migrate,h.p.migrate=h.p.migrate,
        gamma.mult=gamma.mult)
      try.num <- try.num + 1
    }
    if (end.no.migrate) h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc),
      cores=cores,report=report,gamma.mult=gamma.mult)
  } else {
    os <- get.os()
    if ( cores>1 & os=="windows") {
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
      sfClusterSetupRNG()
      sfLibrary(tnorm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfExportAll()
      samples <- sfLapply(samples,run.unstuck.dmc,nmc=nmc,report=report,cores=1,
        cut=cut,nbad=nbad,max.try=max.try,p.migrate=p.migrate,
        gamma.mult=gamma.mult,verbose=verbose)
      if (cores>1) sfStop()
    } else if (cores>1) {
      require(parallel, quietly=TRUE)
      samples <- mclapply(samples,run.unstuck.dmc,nmc=nmc,
        report=report,cores=1,cut=cut,nbad=nbad,max.try=max.try,
        p.migrate=p.migrate,gamma.mult=gamma.mult,verbose=verbose,
        mc.cores=cores)
    } else {
      samples <- lapply(samples,run.unstuck.dmc,nmc=nmc,report=report,
        cores=1,gamma.mult=gamma.mult,verbose=verbose,
        cut=cut,nbad=nbad,max.try=max.try,p.migrate=p.migrate)
    }
  }
  samples
}


#' @rdname h.run.dmc
#' @export
h.run.converge.dmc <- function(samples,nmc=NA,report=10,cores=1,gamma.mult=2.38,
  h.gamma.mult=NA,random.phi=TRUE,random.theta=TRUE,cut=1.1,max.try=20,minN=NA,meanN=NA,
  slaveOutfile=NULL,digits=2,transform=TRUE,autoburnin=FALSE,split=TRUE,
  save="",verbose=TRUE,thorough=FALSE,
  finalrun=FALSE,finalI=NA,finalminN=NA,finalmeanN=NA,addtofinal=FALSE,removefromfinal=NA)
  # Like run.unstuck but applied to list of subjects.
  # If !thorough and hyper then removes based on only hyper gelman.diag and
  # effective sample size at hyper level only but criterion on gelman.diag
  # applied to hyper AND subjects.
  # If thorough and hyper will check every individual parameter Rhat for both
  # hyper and subjects, exiting when maximum is less than cut. Effective sample
  # size same as not thorough.
  # If finalrun=TRUE gets final fresh set of samples fulfilling finalI/minN/meanN
  # If addtofinal adds on to previous, and can also removefromfinal=1:N before doing so
  # h.save will save each random effects (heirarchical) try to an .RData file
{

  get.gd <- function(samples,thorough=TRUE,autoburnin=FALSE,transform=TRUE,split=TRUE,verbose=FALSE) {
    if (thorough)
      c(gelman.diag.dmc(samples,hyper=TRUE)$psrf[,1],
        unlist(lapply(gelman.diag.dmc(samples),function(x){x$psrf[,1]}))) else
          h.gelman.diag.dmc(samples,autoburnin=autoburnin,
            transform=transform,split=split,verbose=FALSE)
  }

  if ( !is.null(samples$theta) )
    stop("For a single subject use run.converge.dmc")
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc
  if ( any(is.na(samples[[1]]$theta)) )
    samples <- h.run.dmc(samples=samples,report=report,cores=cores,
      gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult,
      random.phi=random.phi,random.theta=random.theta)
  if ( any(names(attributes(samples))=="hyper") ) {
    try.num <- 0
    if (!is.na(minN) & !is.na(meanN)) {
      warning("Both minN and meanN specified, using minN")
      meanN <- NA
    }
    if ( !is.na(minN) ) nfun <- "min"
    if (!is.na(meanN)) {
      nfun <- "mean"
      minN <- meanN
    }

    if ( finalrun & is.na(finalI) ) {
      if ( is.na(finalminN) & is.na(finalmeanN) )
        stop("Must specify a finalI or a finalminN or a finalmeanN if finalrun=TRUE")
      if ( !is.na(finalminN) & !is.na(finalmeanN) ) {
        warning("Both finalminN and finalmeanN specified, using finalminN")
        finalmeanN <- NA
      }
      if ( !is.na(finalminN) ) finalfun <- "min"
      if (!is.na(finalmeanN)) {
        finalfun <- "mean"
        finalminN <- finalmeanN
      }
    }

    gd <- get.gd(samples,thorough=thorough,autoburnin=autoburnin,
      transform=transform,split=split,verbose=FALSE)
    if ( !all(gd < cut) | is.na(minN) ) effectiveN <- NA else
      effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
    if (verbose) {
      if (thorough)
        cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
          ", Effective N = ",effectiveN,", Maximum psrf = ",
          round(max(gd),digits=digits),sep="")) else
          {
            cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
              ", Effective N = ",effectiveN,", Hyper mpsrf = ",
              round(gd["hyper"],digits=digits),
              "\nSubject mpsrf achieved (sorted):\n",sep=""))
            print(round(sort(gd[names(gd)!="hyper"],decreasing=TRUE),digits=digits))
          }
    }
    if ( !(all(gd < cut)) | (!is.na(effectiveN) && effectiveN < minN) ) {
      repeat {
        try.num <- try.num + 1
        samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc,add=TRUE),
          cores=cores,report=report,gamma.mult=gamma.mult,
          h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
        gd <- get.gd(samples,thorough=thorough,autoburnin=autoburnin,
          transform=transform,split=split,verbose=FALSE)
        shorter <- h.samples.dmc(samples=samples,remove=1:nmc,nmc=0,add=TRUE)
        gd.short <- get.gd(shorter,thorough=thorough,autoburnin=autoburnin,
          transform=transform,split=split,verbose=FALSE)
        if (thorough) shorten <- max(gd.short) < max(gd) else
          shorten <- gd.short["hyper"] < gd["hyper"]
        if ( shorten ) {
          samples <- shorter
          gd <- gd.short
          if (verbose) cat(paste("Discarding initial",nmc,"samples.\n"))
        }
        if (verbose) {
          if (thorough)
            cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
              ", Effective N = ",effectiveN,", Maximum psrf = ",
              round(max(gd),digits=digits),sep="")) else
              {
                cat(paste("Iterations = ",dim(attr(samples,"hyper")$phi[[1]])[3],
                  ", Effective N = ",effectiveN,", Hyper mpsrf = ",
                  round(gd["hyper"],digits=digits),
                  "\nSubject mpsrf achieved (sorted):\n",sep=""))
                print(round(sort(gd[names(gd)!="hyper"],decreasing=TRUE),digits=digits))
              }
        }
        if ( try.num >= max.try ) break
        if ( all(gd <= cut) ) {
          if ( is.na(minN) ) break else {
            effectiveN <- do.call(nfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
            if ( effectiveN > minN ) break
          }
        }
        if (save != "") save(samples,file=paste(save,"RData",sep="."))
      }
    }
    if ( addtofinal && (finalrun & !any(is.na(removefromfinal))) &&
        (max(removefromfinal)>dim(samples[[1]]$theta)[3]) ) {
      finalrun <- FALSE
      warning("Final run aborted as cutfromfinal to large")
    }
    if (finalrun) {
      try.num <- 0
      if ( addtofinal & !any(is.na(removefromfinal)) )
        samples <- h.samples.dmc(samples=samples,nmc=0,add=TRUE,remove=removefromfinal)
      if (verbose) cat("\n\nDoing final run\n")
      if ( !is.na(finalI) ) {
        if ( addtofinal ) {
          finalI <- finalI-dim(samples[[1]]$theta)[3]
          if (finalI > 0) {
            samples <- h.samples.dmc(samples=samples,nmc=finalI,add=TRUE)
            samples <- h.run.dmc(samples,cores=cores,report=report,gamma.mult=gamma.mult,
              h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
          } else {
            warning("Already more than finalI available")
          }
        } else samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=finalI),
          cores=cores,report=report,gamma.mult=gamma.mult,
          h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
      } else {
        samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc,add=addtofinal),
          cores=cores,report=report,gamma.mult=gamma.mult,
          h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
        effectiveN <- do.call(finalfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
        if (verbose) cat(paste("Iterations = ",
          dim(attr(samples,"hyper")$phi[[1]])[3],", Effective N = ",effectiveN,"\n\n",sep=""))
        if ( effectiveN < finalminN ) {
          repeat {
            try.num <- try.num + 1
            samples <- h.run.dmc(samples=h.samples.dmc(samples=samples,nmc=nmc,add=TRUE),
              cores=cores,report=report,gamma.mult=gamma.mult,
              h.gamma.mult=h.gamma.mult,random.phi=random.phi,random.theta=random.theta)
            if ( try.num >= max.try ) break
            effectiveN <- do.call(finalfun,list(effectiveSize.dmc(samples,hyper=TRUE)))
            if (verbose) cat(paste("Iterations = ",
              dim(attr(samples,"hyper")$phi[[1]])[3],", Effective N = ",effectiveN,"\n\n",sep=""))
            if ( effectiveN > finalminN ) break
            if (save != "") save(samples,file=paste(save,"RData",sep="."))
          }
        }
      }
    }
  } else {
    os <- get.os()
    if ( cores>1 & os=="windows") {
      require(snowfall,quietly=TRUE)
      require(rlecuyer,quietly=TRUE)
      sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
      sfClusterSetupRNG()
      sfLibrary(tnorm)
      sfLibrary(rtdists)
      sfLibrary(statmod)
      sfLibrary(pracma)
      sfLibrary(coda)
      sfExportAll()
      samples <- sfLapply(samples,run.converge.dmc,report=report,cores=1,
        gamma.mult=gamma.mult,cut=cut,max.try=max.try,minN=minN,meanN=meanN,
        nmc=nmc,transform=transform,autoburnin=autoburnin,split=split,
        verbose=verbose)
      sfStop()
    } else if (cores>1) {
      require(parallel, quietly=TRUE)
      samples <- mclapply(samples,run.converge.dmc,report=report,
        cores=1,gamma.mult=gamma.mult,cut=cut,max.try=max.try,minN=minN,
        meanN=meanN,nmc=nmc,transform=transform,autoburnin=autoburnin,
        split=split,verbose=verbose,mc.cores=cores)
    } else {
      samples <- lapply(samples,run.converge.dmc,report=report,cores=cores,
        gamma.mult=gamma.mult,cut=cut,max.try=max.try,minN=minN,meanN=meanN,
        nmc=nmc,transform=transform,autoburnin=autoburnin,split=split,verbose=verbose)
    }
  }
  samples
}


### Posterior predictive


#' @export
h.post.predict.dmc <- function(samples,n.post=100,probs=c(1:99)/100,bw="nrd0",
  save.simulation=FALSE,factors=NA,save.subject.posts=FALSE,cores=1,ignore.R2=FALSE,
  probs.gglist =c(0.1, 0.5, 0.9), CI.gglist =c(0.025, 0.975),censor=c(NA,NA))
  # apply post.predict to each subject
{
  os <- get.os()
  if ( cores>1 & length(samples)>1 & os=="windows") {
    cat("No progress indication in multi-core\n")
    require(snowfall,quietly=TRUE)
    require(rlecuyer,quietly=TRUE)
    sfInit(parallel=TRUE, cpus=cores, type="SOCK")
    sfClusterSetupRNG()
    sfLibrary(msm)
    sfLibrary(rtdists)
    sfExportAll()
    out <- sfLapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
      factors=factors,save.simulation=save.simulation,gglist=TRUE,
      save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
    sfStop()
  } else if (cores>1) {
    cat("No progress indication in multi-core\n")
    require(parallel, quietly=TRUE)
    out <- mclapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
      factors=factors,save.simulation=save.simulation,gglist=TRUE,
      save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
  } else {
    out <- lapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
      factors=factors,save.simulation=save.simulation,gglist=TRUE,
      save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
  }

  if ( !save.simulation ) { # Get averages
    sim <- do.call(rbind,lapply(out,function(x){attr(x,"sim")}))
    if ( (any(names(samples[[1]]$data)=="R2")) && !ignore.R2 ) for (i in 1:length(samples)) {
      samples[[i]]$data$R <-
        paste(as.character(samples[[i]]$data$R),as.character(samples[[i]]$data$R2),sep="")
      samples[[i]]$data$R[samples[[i]]$data$R2=="DK"] <- "DK"
      samples[[i]]$data$R <- factor(samples[[i]]$data$R)
    }
    dat <- do.call(rbind,lapply(samples,function(x){x$data}))
    facs <- names(attr(attributes(samples[[1]]$data)$model,"factors"))
    if (!is.null(factors)) {
      if (any(is.na(factors))) factors <- facs
      if (!all(factors %in% facs))
        stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
    }
    sim.dqp <- get.dqp(sim[,-1],factors,probs,n.post=1,bw=bw)
    dat.dqp <- get.dqp(sim=dat,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    av <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[sim$reps==i,-1]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(av,"dpqs") <- dpqs
    # Strip global sim attribute and dpqs for each participant
    out <- lapply(out,function(x){
      attr(x,"sim") <- NULL
      if (!save.subject.posts) attr(x,"dpqs") <- NULL
      x
    })
    attr(av, "gglist") <- get.fitgglist.dmc(sim,dat,factors=factors, noR=FALSE,
      quantiles.to.get = probs.gglist, CI= CI.gglist)
    attr(out,"av") <- av
  }
  out
}


