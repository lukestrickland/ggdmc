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


