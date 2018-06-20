#' Plot Cell Density
#'
#' \code{plot.dist} plots the distributions of correct and error RTs and
#' reports accuracy rate.
#'
#' @param mdi a model data instance, a data frame class created by model.dmc
#' @param xlim the range on the x axis. Default is (0, Inf)
#' @param ylim the range on the y axix. Default is (0, Inf)
#' @param main a string for the main title
#' @param save.dat save plot data, instead of plotting
#' @param ... other arguments
#' @keywords plot_dist
#' @import ggplot2
#' @export
#' @examples
#' model <- model.dmc(
#'   p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'   constants = c(st0=0,d=0),
#'   match.map = list(M=list(s1="r1",s2="r2")),
#'   factors   = list(S=c("s1","s2")),
#'   responses = c("r1","r2"),
#'   type      = "rd")
#'
#' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat1 <- simulate(model, nsim=1e2, p.vector=pVec)
#' mdi1 <- BindDataModel(dat1, model)
#' plot_dist(mdi1)
plot_dist <- function(data, xlim=c(0, Inf), ylim=c(0, Inf), main=NULL,
  save.dat=FALSE, ... ) {

  if (is.data.frame(data)) {
    model    <- attr(data, "model")
    matchMap <- attr(model, "match.map")$M
    x0 <- list()
    for(i in 1:length(matchMap)) {
      stimulus <- names(matchMap)[i]
      response <- matchMap[i]
      tmp <- censor(data[data$S==stimulus, ], xlim = xlim)
      df0 <- density.dmc(tmp, C = response)
      df0$S <- stimulus
      x0[[i]] <- df0
    }
    x1 <- do.call(rbind.data.frame, x0)
  } else {
    stop("Multiple participants yet to be implemented\n")
  }

  if (is.infinite(xlim[2])) xlim <- range(x1$x)
  if (is.infinite(ylim[2])) ylim <- range(x1$y)

  acc.C <- numeric(length(x0))
  for (i in 1:length(x0)) { acc.C[i] <- attr(x0[[i]], "accuracy") }
  acc      <- data.frame(C=acc.C, S=names(matchMap))
  ann_data <- data.frame(
    xpos = xlim[2]*(4/5),
    ypos = c(ylim[2]*3/5, ylim[2]/2),
    lpos = paste0("Accuracy (", acc[,2], ")=", round(acc[,1], 2)),
    S    = acc[,2])

  gpline <- aes_string(colour="C", linetype="S")

  ## facets <- facet_grid( formula(paste(". ~", "S")) )
  p <- ggplot(x1, aes_string(x="x", y="y")) +
    geom_line(mapping=gpline,  size=.7) +
    scale_y_continuous(name = "Density" ) +
    scale_x_continuous(name = "RT") +
    scale_color_grey(start=.3, end=.8) +
    coord_cartesian(ylim= c(0, ylim[2]+.05) ) +
    geom_text(aes_string(x="xpos", y="ypos", label="lpos"), ann_data) +
    labs(title=main)
  if (save.dat==TRUE) {
    return(x1)
  } else {
    print(p)
  }
}


#' Posterior Predictive Plot
#'
#' This is posterior predictive plot ggplot2 version
#'
#' @param x post predictive object created usually by
#' \code{post.predict.ggdmc}.
#' @param y default NULL. No function. Just to make it compatible to
#' \code{plot}
#' @param style to plot probability density, \code{pdf}, cumulative
#' density, \code{cdf}, or \code{both}
#' @param gpvar1 add one more experimental factor for faceting
#' @param mode an argument
#' @param size font size
#' @param legPos where to anchor legend, in percentage of the x and y axes
#' @param legKeyS legend key size
#' @param legKeyW legend key width
#' @param xlim x axis range
#' @param ... other arguments
#' @export
#' @import ggplot2
plot.pp.ggdmc <- function(x, y = NULL, style = "pdf", gpvar1=NULL, mode="mode",
  size = 18, legPos = c(.85,.85), legKeyS=unit(2, "lines"),
  legKeyW=unit(1, "cm"), xlim=c(0, Inf), ... ) {

  nlty <- length(unique(x$mode))
  gp1line1 <- aes_string(color="R", linetype=mode)

  if(style == "both") {
    pcdf.string <- "pcdf"
    p0 <- ggplot(x, aes(x, y)) +
      geom_line (mapping=gp1line1,  size=.7) +
      scale_color_grey(start=.1, end=.5) +
      scale_y_continuous(name = "Density" ) +
      scale_x_continuous(name = "RT") +
      scale_linetype_manual(values=1:nlty) +
      theme(legend.direction="horizontal",
            legend.position =legPos)

    if(!is.null(gpvar1)) {
      facets <- facet_grid( paste0(gpvar1,"+", "S", "~", pcdf.string) )
    } else {
      facets <- facet_grid( paste0("S", "~", pcdf.string) )
    }

    print(p0 + facets)
  } else {
    dtmp <- x[x$pcdf %in% style,]
    nlty <- length(unique(dtmp$mode))
    ## pcdf.string <- style

    if(!is.null(gpvar1)) {
      facets <- facet_grid( paste0(gpvar1,"+", "S", "~ .") )
    } else {
      facets <- facet_grid( paste0(". ~", "S"))
    }

    p0 <- ggplot(dtmp, aes(x, y)) +
      geom_line (mapping=gp1line1,  size=.7) +
      scale_color_grey(start=.1, end=.5) +
      scale_linetype_manual(values=1:nlty) +
      scale_y_continuous(name = "RT" ) +
      scale_x_continuous(name = "Density") +
      theme(legend.direction="horizontal",
            legend.position =legPos)
    print(p0 + facets)
  }
}

#' Plot an Autocorrelation Matrix
#'
#' This is a wrapper function to plot a DMC object, using \pkg{ggmcmc}
#' \code{ggs_autocorrelation}
#'
#' @param x a DMC sample/object
#' @param lag.max to make it compartible to \code{acf} of \code{stats} pacakge.
#' Default is NULL
#' @param start plot from which iteration. Default is from the first iteration.
#' @importFrom ggmcmc ggs ggs_autocorrelation
#' @export
acf.dmc <- function(x, lag.max=NULL, start=1) {
  d <- ggs( mcmc.list.dmc(x, start=start) )
  ggs_autocorrelation(d) +
    theme(legend.position="none") +
    theme(axis.text.x  = element_blank())
}

#' Create a Plot Matrix of Posterior Simulations
#'
#' This is a wrapper function to plot a DMC object, using \pkg{ggmcmc}
#' \code{ggs_pairs}
#'
#' @param x a DMC sample/object
#' @param start plot from which iteration. Default is from the first iteration.
#' @param ... other arguments
#' @export
#' @importFrom ggmcmc ggs ggs_pairs
pairs.dmc <- function(x, start=1, ...) {
  d <- ggs( mcmc.list.dmc(x, start=start) )
  print( ggs_pairs(d, lower=list(continuous="density")) )
}

#' Plot DMC Samples
#'
#' Plot trace and probability desntiy, using a model samples.
#'
#' If a samples with hyper attribute is set, plot.hyper will be called.
#' If pll.chain=TRUE changes samples$pll to an mcmc object and plots
#' posterior log-likelihood of chains.
#'
#' @param x a \code{run.dmc} or \code{samples.dmc} generated model samples
#' @param y default NULL. No function. Just to make it compatible to
#' \code{plot}
#' @param start instruct the function to plot starting from which iteration.
#' This indicates how many burn-in interations one requests. For example,
#' start=101, indicates 100 burn-in interations.
#' @param end instruct the function to plot ending at a certain iteration
#' @param save.ll a boolean switch to tell the function to save the mean
#' log-likelihood. This option does not work in DMC's plot.dmc, too.
#' @param main.pll a string as the title for the boxplot. Default is NULL
#' @param pll.chain a boolean switch to plot posterior log likelihoood
#' @param pll.together a boolean switch to plot the posterior log-likelihood
#' chains all together in one canvar
#' @param pll.barplot a boolean switch to plot the means of posterior
#' log-likelihood of all chains as a barplot. By default, it is off.
#' @param only.prior Default is FALSE
#' @param only.like Default is FALSE. only.prior and only.like two switches
#' to plot only prior density or only log-likelihood probability.
#' @param smooth default FALSE
#' @param density plot probability density together with trace? Default FALSE
#' @param save.dat whether save the internal data table out for polish plots
#' @param p.prior prior distribution setting. necessary for plot.prior to work
#' @param natural additional argument for plot.prior
#' @param trans additional argument for plot.prior
#' @param xlim additional argument for plot.prior
#' @param chain1 plot all chains or just chain1
#' @param ... other arguments
#' @keywords plot.dmc
#' @export
#' @import ggplot2
#' @importFrom coda mcmc mcmc.list
#' @importFrom ggmcmc ggs ggs_density ggs_traceplot
#' @importFrom ggthemes theme_fivethirtyeight theme_wsj
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics barplot
#' @examples
#' m1 <- model.dmc(
#'   p.map=list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'   constants=c(st0=0,d=0),
#'   match.map=list(M=list(s1="r1",s2="r2")),
#'   factors=list(S=c("s1","s2")),
#'   responses=c("r1","r2"),
#'   type="rd")
#'
#' p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#'
#' dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)
#' mdi1 <- BindDataModel(dat1, m1)
#'
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1=c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2=c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower=c(0,-5, 0, 0, 0, 0),
#'   upper=c(5, 7, 2, 2, 2, 2))
#'
#' samples0 <- samples.dmc(nmc=100, p.prior=p.prior, data=mdi1)
#' samples0 <- run.dmc(samples0, p.migrate=.05)
#'
#' ## Windows tests produce grid.Call problems
#' ## Use should use with caution.
#' ## plot(samples0)
#' ## plot(samples0, density=TRUE)
#' ## plot(samples0, start=101)
#' ## plot(samples0, start=101, density=TRUE)
#' ## plot(samples0, pll.chain=TRUE)
plot.dmc <- function(x, y=NULL, start=1, end=NA, save.ll=FALSE, main.pll=NULL,
  pll.chain = TRUE, pll.together=TRUE, pll.barplot=FALSE, only.prior=FALSE,
  only.like=FALSE, smooth=FALSE, density=FALSE, save.dat=FALSE,
  p.prior=NULL, natural=TRUE, trans=NA, xlim=NA, chain1=TRUE, ...)
{

  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  if ( pll.chain | pll.barplot ) {
    if (only.prior) {
      chain.pll <- x$summed_log_prior[start:end,]
    } else if (only.like) {
      chain.pll <- x$log_likelihoods[start:end,]
    } else {
      chain.pll <- x$summed_log_prior[start:end,]+x$log_likelihoods[start:end,]
    }

    colnames(chain.pll) <- 1:dim(chain.pll)[2]
    if (pll.barplot) {
      mean.ll        <- colMeans(chain.pll)
      names(mean.ll) <- 1:length(mean.ll)
      barplot(mean.ll, ylab="Mean Post Like", main=main.pll)
      if (save.ll) return(mean.ll)
    } else {
      if (!pll.together) {
        d <- mcmc(chain.pll)
      } else {
        d <- coda::mcmc.list(lapply(data.frame(chain.pll),
          function(xx){mcmc(as.matrix(xx))}))
      }
    }
  } else {
    d <- mcmc.list.dmc(x, start=start, end=end)
  }


  DT <- ggs(d)
  DT$Chain <- factor(DT$Chain)

  f1 <- ggs_traceplot(DT) + theme_fivethirtyeight() +
    theme(legend.position="none")
  if(!pll.together) {
    f1 <- f1 + facet_wrap(Chain~Parameter) +
      theme(strip.text.x = element_blank())
  }

  if (density) {
    ## cat("When density is TRUE, plotting prior density is disable.")
    f2 <- ggs_density(DT) +
      ylab("Density")+ theme_wsj()+
      theme(legend.position="none")

    if (!pll.together) { f2 <- f2 + facet_wrap(Chain~Parameter) +
      theme(strip.text.x = element_blank())
    }


    ## Output 1
    print( grid.arrange(f1, f2, ncol=2, nrow=1) )
    if(save.dat) {return(DT)}
  } else if (!is.null(p.prior)) {
    df1 <- NULL
    for(i in 1:length(p.prior)) {
      xlim <- NA
      if (attr(p.prior[[i]], "dist") != "constant") {

        if (any(is.na(trans))) {
          if (natural) trans <- attr(p.prior[[i]], "untrans")
        } else { trans <- "identity" }

        if (is.numeric(i)) i <- names(p.prior)[i]
        if (!(i %in% names(p.prior))) stop("Parameter not in prior")
        p <- p.prior[[i]]
        p$log <- FALSE
        niter <- attr(DT, "nIterations")

        if (is.na(xlim)) xlim <- range(DT[DT$Parameter==i, "value"])
        p$x <- seq(xlim[1], xlim[2], length.out=niter)
        priorx <- do.call(trans, list(x=p$x))
        priory <- do.call(paste0("d", attr(p, "dist")), p)
        df0 <- data.frame(Iteration=1:niter, Chain=factor(0),
          Parameter=factor(i), x=priorx, y=priory)
        df1 <- rbind(df1, df0)
      }
    }
    ## DT <- rbind(DT, df1)
    ## Overwrite f2 to add chain 0 (ie prior chain)
    if (chain1) {
      DT1 <- DT[DT$Chain==1,]
      attr(DT1, "nChains") <- attr(DT, "nChains")
      f2 <- ggs_density(DT1) + ylab("Density")+ theme_wsj() +
        geom_line(aes(x, y), df1)
    } else {
      f2 <- ggs_density(DT) + ylab("Density")+ theme_wsj() +
        geom_line(aes(x, y), df1)
    }

    ## Output 2
    invisible(print(f2))

    if(save.dat) {return(list(DT, df1))}
  } else {
    ## Output 3
    invisible(print(f1))
    if(save.dat) {return(DT)}
  }
}

#' Plot a DMC Sample with Multiple Participants
#'
#' Plot trace and probability desntiy, using a model samples. This function
#' taks a sample list.
#'
#' @param x a \code{run.dmc} or \code{samples.dmc} generated model samples
#' @param y default NULL. No function. Just to make it compatible to
#' \code{plot}
#' @param start instruct the function to plot starting from which iteration.
#' This indicates how many burn-in interations one requests. For example,
#' start=101, indicates 100 burn-in interations.
#' @param end instruct the function to plot ending at a certain iteration
#' @param save.ll a boolean switch to tell the function to save the mean
#' log-likelihood. This option does not work in DMC's plot.dmc, too.
#' @param main.pll a string as the title for the boxplot. Default is NULL
#' @param pll.chain a boolean switch to plot posterior log likelihoood
#' @param pll.together a boolean switch to plot the posterior log-likelihood
#' chains all together in one canvar
#' @param pll.barplot a boolean switch to plot the means of posterior
#' log-likelihood of all chains as a barplot. By default, it is off.
#' @param only.prior Default is FALSE
#' @param only.like Default is FALSE. only.prior and only.like two switches
#' to plot only prior density or only log-likelihood probability.
#' @param smooth default FALSE
#' @param density plot probability density together with trace? Default FALSE
#' @param save.dat whether save the internal data table out for polish plots
#' @param p.prior prior distribution setting. necessary for plot.prior to work
#' @param natural additional argument for plot.prior
#' @param trans additional argument for plot.prior
#' @param xlim additional argument for plot.prior
#' @param chain1 plot all chains or just chain1
#' @param subject which subject in the list to plot. Default the 1st one.
#' @param ... other arguments
#' @export
#' @import ggplot2
#' @importFrom coda mcmc mcmc.list
#' @importFrom ggmcmc ggs ggs_density ggs_traceplot
#' @importFrom ggthemes theme_fivethirtyeight theme_wsj
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics barplot
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
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' dat <- h.simulate.dmc(m1, p.prior=pop.prior, n=50, ns=4)
#' mdi <- BindDataModel(dat, m1)
#' p.prior  <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05) * 5,
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## Fixed-effect model
#' samplesInit <- h.samples.dmc(nmc=50, p.prior=p.prior, data=mdi, thin=1)
#' samples0    <- h.run.dmc(samples=samplesInit, report=25)
#' ## Windows tests produce a grid.Call problem. The user should use
#' ## with caution.
#' ## plot(samples0) ## traceplot for the first participant
#' ## plot(samples0, density=TRUE)  ## trace- and density-plot
#' ## plot(samples0, density=TRUE, subject=2) ## Plot second participant
#' ## plot(samples0, density=TRUE, subject=3, start=101) ## From 101 iteration
#' ## Plot iteratoin 201 to 400
#' ## plot(samples0, density=TRUE, subject=4, start=201, end=400)
plot.dmc.list <- function(x, y=NULL, start=1, end=NA, save.ll=FALSE,
  main.pll=NULL, pll.chain=FALSE, pll.together=TRUE, pll.barplot=FALSE,
  only.prior=FALSE, only.like=FALSE, smooth=FALSE, density=FALSE,
  save.dat=FALSE, p.prior=NULL, natural=TRUE, trans=NA,
  xlim=NA, chain1=TRUE, subject=1, ...)
{
  x <- x[[subject]]
  if ( is.na(end) ) end <- x$nmc
  if ( end <= start ) stop("End must be greater than start")
  if ( pll.chain | pll.barplot ) {
    if (only.prior) {
      chain.pll <- x$summed_log_prior[start:end,]
    } else if (only.like) {
      chain.pll <- x$log_likelihoods[start:end,]
    } else {
      chain.pll <- x$summed_log_prior[start:end,] +
        x$log_likelihoods[start:end,]
    }

    colnames(chain.pll) <- 1:dim(chain.pll)[2]
    if (pll.barplot) {
      mean.ll        <- colMeans(chain.pll)
      names(mean.ll) <- 1:length(mean.ll)
      barplot(mean.ll, ylab="Mean Post Like", main=main.pll)
      if (save.ll) return(mean.ll)
    } else {
      if (!pll.together) {
        d <- mcmc(chain.pll)
      } else {
        d <- coda::mcmc.list(lapply(data.frame(chain.pll),
          function(xx){mcmc(as.matrix(xx))}))
      }
    }
  } else {
    d <- mcmc.list.dmc(x, start=start, end=end)
  }


  DT <- ggs(d)
  DT$Chain <- factor(DT$Chain)

  f1 <- ggs_traceplot(DT) + theme_fivethirtyeight() +
    theme(legend.position="none")
  if(!pll.together) {
    f1 <- f1 + facet_wrap(Chain~Parameter) +
      theme(strip.text.x = element_blank())
  }

  if (density) {
    ## cat("When density is TRUE, plotting prior density is disable.")
    f2 <- ggs_density(DT) +
      ylab("Density")+ theme_wsj()+
      theme(legend.position="none")

    if (!pll.together) { f2 <- f2 + facet_wrap(Chain~Parameter) +
      theme(strip.text.x = element_blank())
    }


    ## Output 1
    print( grid.arrange(f1, f2, ncol=2, nrow=1) )
    if(save.dat) {return(DT)}
  } else if (!is.null(p.prior)) {
    df1 <- NULL
    for(i in 1:length(p.prior)) {
      xlim <- NA
      if (attr(p.prior[[i]], "dist") != "constant") {

        if (any(is.na(trans))) {
          if (natural) trans <- attr(p.prior[[i]], "untrans")
        } else { trans <- "identity" }

        if (is.numeric(i)) i <- names(p.prior)[i]
        if (!(i %in% names(p.prior))) stop("Parameter not in prior")
        p <- p.prior[[i]]
        p$log <- FALSE
        niter <- attr(DT, "nIterations")

        if (is.na(xlim)) xlim <- range(DT[DT$Parameter==i, "value"])
        p$x <- seq(xlim[1], xlim[2], length.out=niter)
        priorx <- do.call(trans, list(x=p$x))
        priory <- do.call(paste0("d", attr(p, "dist")), p)
        df0 <- data.frame(Iteration=1:niter, Chain=factor(0),
          Parameter=factor(i), x=priorx, y=priory)
        df1 <- rbind(df1, df0)
      }
    }
    ## DT <- rbind(DT, df1)
    ## Overwrite f2 to add chain 0 (ie prior chain)
    if (chain1) {
      DT1 <- DT[DT$Chain==1,]
      attr(DT1, "nChains") <- attr(DT, "nChains")
      f2 <- ggs_density(DT1) + ylab("Density")+ theme_wsj() +
        geom_line(aes(x, y), df1)
    } else {
      f2 <- ggs_density(DT) + ylab("Density")+ theme_wsj() +
        geom_line(aes(x, y), df1)
    }

    ## Output 2
    print(f2)
    if(save.dat) {return(list(DT, df1))}
  } else {
    ## Output 3
    print(f1)
    if(save.dat) {return(DT)}
  }
}


#' Plot DMC Samples at the Hyper level
#'
#' Plot trace and probability desntiy. This method provides a switch to choose
#' plotting at the data or hyper level. The samples must be an object fit by
#' a hierarchical model.
#'
#' @param x \code{run.dmc} or \code{samples.dmc} generated model samples
#' @param y default NULL. No function. Just to make it compatible to
#' \code{plot}
#' @param hyper a boolean switch to draw hyper-parameter. By default it is off.
#' @param start instruct the function to plot starting from which iteration.
#' This indicates how many burn-in interations one requests. For example,
#' start=101, indicates 100 burn-in interations.
#' @param end instruct the function to plot ending at a certain iteration
#' @param save.ll a boolean switch to tell the function to save the mean
#' log-likelihood. This option does not work in DMC's plot.dmc, too.
#' @param main.pll a string as the title for the boxplot. Default is NULL
#' @param pll.chain a boolean switch to plot posterior log likelihoood
#' @param pll.together a boolean switch to plot the posterior log-likelihood
#' chains all together in one canvar
#' @param pll.barplot a boolean switch to plot the means of posterior
#' log-likelihood of all chains as a barplot. By default, it is off.
#' @param only.prior only.prior and only.like pick out one or other
#' component of pll.
#' @param only.like a boolean switch.
#' @param subject indicate which subject to plot in a multi-subject list
#' samples
#' @param smooth default FALSE
#' @param density plot probability density too?
#' @param save.dat whether save the internal data table out for polish plots
#' @param p.prior prior distribution setting. necessary for plot.prior to work
#' @param natural additional argument for plot.prior
#' @param trans additional argument for plot.prior
#' @param xlim additional argument for plot.prior
#' @param chain1 plot all chains or just chain1 as coda does
#' @param ... other arguments
#' @keywords plot.dmc
#' @export
#' @import ggplot2
#' @importFrom coda mcmc mcmc.list
#' @importFrom ggmcmc ggs ggs_density ggs_traceplot
#' @importFrom ggthemes theme_fivethirtyeight theme_wsj
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics barplot
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
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' dat <- h.simulate.dmc(m1, p.prior=pop.prior, n=50, ns=4)
#' mdi <- BindDataModel(dat, m1)
#' p.prior  <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05) * 5,
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## Fixed-effect model
#' samplesInit <- h.samples.dmc(nmc=100, p.prior=p.prior, data=mdi, thin=1)
#' samples0    <- h.run.dmc(samples=samplesInit, report=20)
#'
#' ## Make a hyper-prior list
#' mu.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05) * 5,
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
#' hsamplesInit <- h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi, thin=1)
#' hsamples0 <- h.run.dmc(samples=hsamplesInit)
#'
#' ## Windows Testing errors
#' ## plot(hsamples0) ## Only first participant
#' ## plot(hsamples0, hyper=TRUE) ## Group-level parameters
#' ## plot(hsamples0, hyper=TRUE, density=TRUE) ## Trace and density plots
#'
#' ## plot(hsamples0, p.prior=p.prior) ## plot prior and posterior
#' ## plot(hsamples0, p.prior=pp.prior, hyper=TRUE) ## hyper-level
#'
#' ## plot(hsamples0, hyper=TRUE, start=101) ## starting from 101th iteration
#' ## plot(hsamples0, hyper=TRUE, density=TRUE, start=101)
#'
#' ## Plot posterior likelihood
#' ## plot(hsamples0, hyper=TRUE, pll.chain=TRUE, start=21)
#'
#' ## Save plot data for modifying
#' ## D <- plot(hsamples0, hyper=TRUE, save.dat=TRUE)
#' ## head(D)
#' ## Source: local data frame [6 x 4]
#' ##   Iteration  Chain Parameter    value
#' ##       (int) (fctr)    (fctr)    (dbl)
#' ## 1         1      1      a.h1 2.184504
#' ## 2         2      1      a.h1 2.184504
#' ## 3         3      1      a.h1 2.184504
#' ## 4         4      1      a.h1 2.501221
#' ## 5         5      1      a.h1 2.501221
#' ## 6         6      1      a.h1 2.501221
plot.hyper <- function(x, y=NULL, hyper=FALSE, start=1, end=NA,
  save.ll=FALSE, main.pll=NULL, pll.chain=FALSE, pll.together=TRUE,
  pll.barplot=FALSE, only.prior=FALSE, only.like=FALSE,
  subject=1, smooth=FALSE, density=FALSE, save.dat=FALSE,
  p.prior=NULL, natural=TRUE, trans=NA,
  xlim=NA, chain1=TRUE, ...)
{
  if ( hyper ) {
    hyper.dat <- attr(x, "hyper")
    if (is.null(hyper.dat)) stop("There are no hyper-parameters to plot.")
    if ( is.na(end) ) end <- hyper.dat$nmc
    if ( end <= start ) stop("End must be greater than start")
    if ( pll.chain | pll.barplot ) {
      chain.pll <- hyper.dat$h_summed_log_prior[start:end,] +
        hyper.dat$h_log_likelihoods[start:end,]
      colnames(chain.pll) <- 1:dim(chain.pll)[2]
      if (pll.barplot) {
        mean.ll <- colMeans(chain.pll)
        names(mean.ll) <- 1:length(mean.ll)
        barplot(mean.ll,ylab="Mean Post Like",main=main.pll)
        if (save.ll) return(mean.ll)
      } else {
        if(!pll.together) {
          d <- coda::mcmc(chain.pll)
        } else {
          d <- coda::mcmc.list(lapply(data.frame(chain.pll),
            function(xx){mcmc(as.matrix(xx))}))
        }
      }
    } else {
      d <- phi.as.mcmc.list(hyper.dat, start=start,end=end)
    }
  } else {
    x <- x[[subject]]
    if ( is.na(end) ) end <- x$nmc
    if ( end <= start ) stop("End must be greater than start")
    if ( pll.chain | pll.barplot ) {
      if (only.prior) {
        chain.pll <- x$summed_log_prior[start:end,]
      } else if (only.like) {
        chain.pll <- x$log_likelihoods[start:end,]
      } else {
        chain.pll <- x$summed_log_prior[start:end,] +
          x$log_likelihoods[start:end,]
      }

      colnames(chain.pll) <- 1:dim(chain.pll)[2]
      if (pll.barplot) {
        mean.ll        <- colMeans(chain.pll)
        names(mean.ll) <- 1:length(mean.ll)
        barplot(mean.ll, ylab="Mean Post Like", main=main.pll)
        if (save.ll) return(mean.ll)
      } else {
        if (!pll.together) {
          d <- mcmc(chain.pll)
        } else {
          d <- coda::mcmc.list(lapply(data.frame(chain.pll),
            function(xx){mcmc(as.matrix(xx))}))
        }
      }
    } else {
      d <- ggdmc::mcmc.list.dmc(x, start=start, end=end)
    }
  }

  DT <- ggs(d)
  DT$Chain <- factor(DT$Chain)

  f1 <- ggs_traceplot(DT) + theme_fivethirtyeight() +
    theme(legend.position="none")

  if(!pll.together) {
    f1 <- f1 + facet_wrap(Chain~Parameter) +
      theme(strip.text.x = element_blank())
  }

  if (hyper) { f1 <- f1 + facet_wrap(~Parameter) }

  if (density) {
    ## cat("When density is TRUE, plotting prior density is disable.\n")
    f2 <- ggs_density(DT) +
      ylab("Density")+ theme_wsj()+
      theme(legend.position="none")

    if (!pll.together) { f2 <- f2 + facet_wrap(Chain~Parameter) +
      theme(strip.text.x = element_blank())
    }

    if (hyper) { f2 <- f2 + facet_wrap(~Parameter) }

    ## Output 1
    print( grid.arrange(f1, f2, ncol=2, nrow=1) )
    if(save.dat) {return(DT)}
  } else if (!is.null(p.prior)) {

    if (length(p.prior) == 2 & hyper) {
      pp.prior <- c(p.prior[[1]], p.prior[[2]])
      if (all(names(p.prior[[1]]) == names(p.prior[[2]])))
        location <- paste(names(p.prior[[1]]),"h1", sep=".")
      scale    <- paste(names(p.prior[[2]]),"h2", sep=".")
      names(pp.prior) <- c(location, scale)
      p.prior <- pp.prior ## overwrite p.prior with pp.prior.new
    }

    df1 <- NULL
    for(i in 1:length(p.prior)) {
      xlim <- NA
      if (attr(p.prior[[i]], "dist") != "constant") {

        if (any(is.na(trans))) {
          if (natural) trans <- attr(p.prior[[i]], "untrans")
        } else { trans <- "identity" }

        if (is.numeric(i)) i <- names(p.prior)[i]
        if (!(i %in% names(p.prior))) stop("Parameter not in prior")
        p <- p.prior[[i]]
        p$log <- FALSE
        niter <- attr(DT, "nIterations")

        if (is.na(xlim)) xlim <- range(DT[DT$Parameter==i, "value"])
        p$x    <- seq(xlim[1], xlim[2], length.out=niter)
        priorx <- do.call(trans, list(x=p$x))
        priory <- do.call(paste0("d", attr(p, "dist")), p)
        df0 <- data.frame(Iteration=1:niter, Chain=factor(0),
          Parameter=factor(i), x=priorx, y=priory)
        df1 <- rbind(df1, df0)
      }
    }

    ## Overwrite f2 to add chain 0 (ie prior chain)
    if (chain1) {
      DT1 <- DT[DT$Chain==1,]
      attr(DT1, "nChains") <- attr(DT, "nChains")
      f2 <- ggs_density(DT1) + ylab("Density")+ theme_wsj() +
        geom_line(aes(x, y), df1)
    } else {
      f2 <- ggs_density(DT) + ylab("Density")+ theme_wsj() +
        geom_line(aes(x, y), df1)
    }

    ## if (hyper) { f2 <- f2 + facet_wrap(~Parameter, scales="free") }

    ## Output 2
    print(f2)
    if(save.dat) {return(list(DT, df1))}
  } else {
    ## Output 3
    print(f1)
    if(save.dat) {return(DT)}
  }
}


#' Profile a DMC Object
#'
#' Profile a DMC model based a data model instance (ie \code{fitted}). For
#' the parameter \code{p.name}, e.g., the boundary separation \emph{a} in a DDM,
#' extract it out from a \code{p.vector} and draws the profile likelihood for
#' data and returns the maximum (on a grid of resolution \code{n.point})
#'
#' Set a range for the x axis, which is the profiled parameter. Initiate a
#' log-likelihodd vector with 0 everywhere and length n.point. keep the other
#' parmaters value fixed and change only the target parameter value to ps[i],
#' and then calculate its sum of log-likelihood. Finally, store it in the i
#' position of ll vector. This function currently works for DDM density only.
#'
#' @param p.name indicate which paramter in theta to plot. For example, in a
#' LBA model with a \code{pVec <- c(A="1",B="1",mean_v="M",sd_v="1",t0="1",
#' st0="1")}, one can assign \code{p.name <- "A"} to ask \code{profile.dmc} to
#' profile \emph{A} parameter
#' @param min.p lower bound for p.name parameter
#' @param max.p upper bound for p.name parameter
#' @param p.vector a parameter vector. Use Lognromal LBA model as an example,
#' \code{pVec <- c(A=.25, B=.35, meanlog_v.true=1, meanlog_v.false=.25,
#' sdlog_v.true=.5,t0=.2)}
#' @param n.point grid for p.name parameter
#' @param digits print out how many digits
#' @param ylim y range
#' @importFrom graphics plot
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
#' dat1 <- simulate(m1, nsim=1e2, p.vector=pVec)
#' mdi1 <- BindDataModel(dat1, m1)
#'
#' ## ---------------------------
#' par(mfrow=c(2,3));
#' profile(mdi1, "a",  .1,   2, pVec)
#' profile(mdi1, "v",  .1,   2, pVec)
#' profile(mdi1, "z",  .2,  .8, pVec)
#' profile(mdi1, "sz", .1,  .9, pVec)
#' profile(mdi1, "sv", .1,   2, pVec)
#' profile(mdi1, "t0", .01, .5, pVec)
#' par(mfrow=c(1,1));
#' @export
profile.dmc <- function(p.name, min.p, max.p, p.vector, data,
  n.point = 100, digits = 2, ylim = NA, nthread = 32) {

  if (!(p.name %in% names(p.vector))) stop("parameter not in p.vector")

  RT        <- data$RT
  model     <- attr(data, "model")
  ise       <- attr(data, "cell.empty")
  allpar    <- attr(model, "all.par")
  parnames  <- attr(model, "par.names")
  type      <- attr(model, "type")
  n1idx     <- attr(model, "n1.order")
  mc        <- attr(model, "match.cell")
  isr1      <- check_rd(type, model)
  cellidx   <- cellIdx2Mat(data)
  pnames    <- names(p.vector)
  matchcell <- attr(model, "match.cell")
  ps        <- seq(min.p, max.p, length.out = n.point)
  nsim      <- attr(data, "n.pda")
  bw        <- attr(data, "bw")
  gpuid     <- attr(data, "gpuid")
  debug     <- attr(data, "debug")

  if (type == "norm") {

    ll <- profile_norm(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
      n1idx, ise, cellidx, RT, matchcell, isr1, p.name, ps)

  } else if (type == "norm_pda") {

    ll <- profile_norm_pda(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, matchcell, isr1, p.name, ps, nsim, bw)

  } else if (type == "norm_pda_gpu") {

    ll <- profile_norm_gpu(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, matchcell, isr1, p.name, ps, nsim, bw)

  } else if (type == "plba0_gpu") {
    ll <- profile_plba0_gpu(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, matchcell, isr1, p.name, ps, nsim, bw, gpuid, nthread,
      debug)

  } else if (type == "plba1") {

    ll <- profile_plba1(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, matchcell, isr1, p.name, ps, nsim, bw)

  } else if (type == "plba1_gpu") {

    ll <- profile_plba1_gpu(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]], n1idx,
      ise, cellidx, RT, matchcell, isr1, p.name, ps, nsim, bw, gpuid, nthread,
      debug)

  } else if (type == "rd") {

    ll <- profile_rd(p.vector, pnames, allpar, parnames, model, type,
      dimnames(model)[[1]], dimnames(model)[[2]], dimnames(model)[[3]],
      n1idx, ise, cellidx, RT, matchcell, isr1, p.name, ps)

  } else {
    stop("Type not defined")
  }

  names(ll) <- round(ps, digits)
  plot(ps, ll, type = "l", xlab = p.name, ylab = "log-likelihood")
  ll[ll==max(ll)]
}


#' Plot Prior Distributions
#'
#' \code{plot_prior} plots the i'th member of the list created by
#' \code{BuildPrior}.  If \code{trans=TRUE}, the function will plot on natural
#' (logarithmic) scale using transform specified in
#' \code{attr(p.prior[[i]], "trans")}. \code{plot.prior} plots all parameters
#' at once.
#'
#' \code{plot_prior} checks if any of the elements in \code{trans} is \code{NA}.
#' It then checks if \code{natural} is set \code{TRUE}
#' (by default it's \code{TRUE}). If \code{natural} is set also
#' \code{TRUE} (i.e., the user wants to do transformation), it then checks what
#' has been set in the "untrans" attribute for the i'th parameter.  Otherwise,
#' its default is set all parameters as \code{identity}.
#'
#' @param i an integer or a character string indicating to plot which parameter
#' @param p.prior a list of list storing prior setting.
#' @param xlim set the range of on x axis. This is usually the range for each
#' parameter. Default is NA.
#' @param natural default TRUE
#' @param n.point default to plot 128
#' @param trans default NA. trans can be a scalar or vector.
#' @param save.dat whether to save the data out
#' @param ... other plotting arguments passing throught dot dot dot.
#' @import ggplot2
#' @export
#' @examples
#' p.prior <- BuildPrior(dists = rep("tnorm", 7),
#'            p1    = c(a=2,   v.f1=4,  v.f2=3,  z=0.5, sv=1,  sz=0.3, t0=0.3),
#'            p2    = c(a=0.5, v.f1=.5, v.f2=.5, z=0.1, sv=.3, sz=0.1, t0=0.05),
#'            lower = c(0,-5, -5, 0, 0, 0, 0),
#'            upper = c(5, 7,  7, 1, 2, 1, 1))
#'
#' plot_prior("a", p.prior)
#' plot_prior(2, p.prior)
#'
#' view(p.prior)
#' plot(p.prior)
#' d <- plot_priors(pop.prior, save.dat=TRUE)
#'
#' ## require(ggplot2)
#' ## p2 <- ggplot(d, aes(x = xpos, y = ypos)) + geom_line() +
#' ##       facet_wrap(~gpvar, scales="free") + theme_bw(base_size =14)
plot_prior <- function(i, p.prior, xlim=NA, natural=TRUE, n.point =128,
                       trans=NA, save.dat=FALSE, ... )
{
  if (any(is.na(trans)))   ### Check 1 ###
  {trans <- ifelse(natural, attr(p.prior[[i]], "untrans"), "identity")}
  if (is.numeric(i)) i <- names(p.prior)[i]   ### Check 2 ###
  ## Stop the function running, if the parameters are out-of-range
  if (!(i %in% names(p.prior))) stop("Parameter not in prior")
  ## Finish checks. Calculate the data frame
  p <- p.prior[[i]]    ## Save all parameter setting to another object called p

  ## Do an educated guess for xlim
  ## xlim can be a scalar or vector. test if any of its elements is NA. If it
  if ( any(is.na(xlim)) ) {
    ## does contain NAs, we set a xlim for the user
    ## check what distribution the user wants to as her prior distribution
    ## 4 options: tnorm (default), beta_lu, gamma_l, and lnorm_l
    ## Basically, we use parameter 1 and 2 to do some (random?) arithmetic
    ## and pick (1) the bigger and smaller number (tnorm); (2) use lower and
    ## upper directly (beta_lu); (3) use lower and a different arithmetic
    ## (gamma_l); (4) use lower and another different arithmetic (lnorm_l)
    xlim <- switch(attr(p, "dist"),
                   tn = c(pmax(p$lower, p[[1]]-3*p[[2]]),
                          pmin(p$upper, p[[1]]+3*p[[2]])),
                   tnorm    = c(pmax(p$lower, p[[1]]-3*p[[2]]),
                                pmin(p$upper, p[[1]]+3*p[[2]])),
                   beta_lu  = c(p$lower, p$upper),
                   gamma_l  = c(p$lower, p[[1]]*p[[2]]+3*sqrt(p[[1]])*p[[2]]),
                   lnorm_l  = c(p$lower, exp(p[[1]]+2*p[[2]])),
                   constant = c(-10,10)
    )
  }

  ## Now we get xlim. Then we want to set all the x points for plotting.
  ## By default we plot 100 point (n.point = 1e2)
  x <- seq(xlim[1], xlim[2], length.out = n.point)
  p$x <- x
  p$log <- FALSE    ## Turn off logarithmic transformation

  ## 1. call a log/identity transform function to change x value
  ## 2. call a density function for a distribution to set y (density) value
  df <- data.frame(
    xpos  = do.call(trans, list(x = x)),
    ypos  = do.call(paste("d", attr(p,"dist"), sep=""), p),
    gpvar = rep( names(p.prior[i]), length(x)))

  if (save.dat==TRUE) {
    return(df)
  } else {
    print(ggplot(df, aes_string(x="xpos", y="ypos")) + geom_line() +
            xlab("Value")+ ylab("Density") + facet_grid(.~gpvar))
  }
}

#' @rdname plot_prior
#' @export
plot.prior <- function(p.prior, save.dat=FALSE, ...) {
  pD <- NULL
  for(j in names(p.prior))
    invisible(pD <- rbind(pD, plot_prior(i=j, p.prior=p.prior, save.dat=TRUE)))
  if (save.dat) {
    return(pD)
  } else {
    p0 <- ggplot(pD, aes_string(x = "xpos", y = "ypos")) +
      geom_line() +
      xlab("")+ ylab("")+
      facet_wrap(~gpvar, scales="free")
    print(p0)
  }
}


# plot.stanfit <- function(fit, pars, start=1, end=NA, pll.chain=FALSE,
#   density = FALSE, save.dat=FALSE)
# {
#   ## plot stan object, using DMC format
#   if (missing(pars)) pars <- fit@sim$pars_oi
#   nmc    <- fit@sim$iter
#   nchain <- fit@sim$chains
#   if (is.na(end)) end <- fit@sim$iter
#   if (end <= start) stop("End must be greater than start")
#
#   if (pll.chain) {
#     pll <- matrix(numeric(nchain*(end - start + 1)), (end - start + 1))
#     for (i in 1:nchain) {
#       chaini  <- fit@sim$samples[[i]]$lp__
#       pll[,i] <- chaini[start:end]
#     }
#     d <- coda::mcmc.list(lapply(data.frame(pll),
#       function(xx){mcmc(as.matrix(xx))}))
#   } else {
#     d <- rstan::As.mcmc.list(fit, pars = pars)
#   }
#
#   DT <- ggmcmc::ggs(d)
#   DT$Chain <- factor(DT$Chain)
#   f1 <- ggs_traceplot(DT) + theme_fivethirtyeight() + theme(legend.position="none")
#
#   if (density) {
#     f2 <- ggs_density(DT) + ylab("Density")+ theme_wsj() + theme(legend.position="none")
#     invisible(print(grid.arrange(f1, f2, ncol=2, nrow=1)))
#     if (save.dat) return(DT)
#   } else {
#     invisible(print(f1))
#     if(save.dat) return(DT)
#   }
#
# }
