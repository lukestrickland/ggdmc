# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to generating graphical output
#    Usually user does not need to edit

#' Plot Distributions for Each Cell
#'
#' If !is.na(C) plots density for correct and error responses for a data
#' frame with columns R (a factor) and RT, adding a boolean score column
#' for R=C. Otherwise plots each response. Can deal with NA in the RT column,
#' in which case it provides a summary of p(NA)
#' @param data.cell a data frame with only onn experimental conditoin
#' @param C a correctness column
#' @param xlim x censor range
#' @param ymax the upper bound for y axis when plotting
#' @param save.density whether to save density data
#' @param digits print how many digits
#' @param main main title for the figure
#' @param show.mean whether to show mean
#' @export
#' @importFrom graphics plot lines legend
#' @examples
#' m1 <- model.dmc(
#'   p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'   constants = c(st0=0,d=0),
#'   match.map = list(M=list(s1="r1",s2="r2")),
#'   factors   = list(S=c("s1","s2")),
#'   responses = c("r1","r2"),
#'   type      = "rd")
#'
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#'
#' dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)
#' mdi1 <- BindDataModel(dat1, m1)
#'
#' ## Accuracy around 70%
#' par(mfrow=c(1,2))
#' plot_cell_density(data.cell=mdi1[mdi1$S=="s1", ], C="r1", xlim=c(0,2))
#' plot_cell_density(data.cell=mdi1[mdi1$S=="s2", ], C="r2", xlim=c(0,2))
#' par(mfrow=c(1,1))
plot.cell.density <- function(data.cell,C=NA,xlim=c(0,Inf),ymax=NA,
  save.density=FALSE,digits=2,main="",show.mean=FALSE)
  # If !is.na(C) plots density for correct and error responses for a data frame
  # with columns R (a factor) and RT, adding a boolean score column for R=C.
  # Otherwise plots each response. Can deal with NA in the RT column, in which
  # case it provides a summary of p(NA)
{
  if (!is.factor(data.cell$R)) data.cell$R <- factor(data.cell$R)
  if (length(C)==1) C <- rep(C,dim(data.cell)[1])
  p.na <- mean(is.na(data.cell$RT))
  is.in <- !is.na(data.cell$RT)
  is.in[is.in] <- data.cell$RT[is.in]>xlim[1] & data.cell$RT[is.in]<xlim[2]
  dat <- data.cell[is.in,]
  if ( !any(is.na(C)) ) {
    if ( is.logical(C) & length(C)==dim(data.cell)[1] )
      dat$C <- C[is.in] else dat$C <- dat$R==C[is.in]
      if (length(dat$RT[dat$C])>2)
        dns.correct <- density(dat$RT[dat$C]) else dns.correct <- NULL
        if (length(dat$RT[!dat$C])>2)
          dns.error <- density(dat$RT[!dat$C]) else dns.error <- NULL
          if (is.null(dns.error) & is.null(dns.correct))
            stop("There are no densities to plot")
          acc <- mean(dat$C)
          if (!is.null(dns.correct))
            dns.correct$y <- dns.correct$y*acc*(1-p.na)
          if (!is.null(dns.error))
            dns.error$y <- dns.error$y*(1-acc)*(1-p.na)
          if (is.na(ymax)) ymax <- max(c(dns.correct$y,dns.error$y))
          if (!is.null(dns.correct)) {
            plot(dns.correct,xlab="RT",ylab="density",ylim=c(0,ymax),main=main)
            if (!is.null(dns.error)) lines(dns.error,col="red")
          } else {
            plot(dns.error,xlab="RT",ylab="density",ylim=c(0,ymax),main=main,col="red")
            if (!is.null(dns.correct)) lines(dns.correct)
          }
          nams <- "Accuracy ="
          ps <- round(acc,digits)
          if (p.na!=0) {
            nams <- c(nams,"p(NA) =")
            ps <- c(ps,round(p.na,digits))
          }
          legend("topright",paste(nams,ps),bty="n")
          legend("topright",xjust=0, inset=c(0,0.1), c("correct","error"),bty="n", lty=c(1,1), col=c("black","red"))
          if ( save.density ) list(correct=dns.correct,error=dns.error)
  } else {
    rs <- levels(dat$R)
    dns <- vector(mode="list",length=length(rs))
    names(dns) <- rs
    ps <- table(dat$R)/dim(dat)[1]
    for (i in rs) {
      rt <- dat$RT[dat$R==i]
      if (length(rt)>2) {
        dns[[i]] <- density(rt)
        dns[[i]]$y <- ps[i]*dns[[i]]$y*(1-p.na)
      }
    }
    ymax <- suppressWarnings(max(unlist(lapply(dns,function(x){max(x$y)}))))
    no.dns <- unlist(lapply(dns,is.null))
    if (all(no.dns))
      stop("There are no densities to plot!")
    dns1 <- dns[!no.dns]
    ltys <- c(1:length(dns1))
    plot(dns1[[1]],xlab="RT",ylab="density",ylim=c(0,ymax),lty=ltys[1],
      main=main)
    if (length(dns1)>1) for (i in 2:length(dns1)) lines(dns1[[i]],lty=ltys[i])
    nams <- paste("p(",names(dns1),") =",sep="")
    ps <- round(ps[!no.dns]*(1-p.na),digits)
    if ( p.na!=0 ) {
      nams <- c(nams,"p(NA) =")
      ps <- c(ps,round(p.na,2))
      lty <- c(ltys,NA)
    }
    legend("topright",paste(nams,ps),lty=ltys,bty="n")
    if ( save.density ) dns
  }
}





# if (!(p.name %in% names(p.vector)))
#   stop("p.name not in p.vector")
# p <- p.vector
# ps <- seq(min.p,max.p,length.out=n.point)
# ll <- numeric(n.point)
# for (i in 1:n.point)
# {
#   p[p.name] <- ps[i]
#   ll[i] <- sum(log(likelihood.dmc(p,data)))
# }


#' @importFrom graphics par
ppl.barplots.dmc <- function(samples,start=1,end=NA,layout=c(5,2))
  # Grid of barplots of pll for set of subjects
{
  if (!is.null(samples$theta))
    stop("This function cannot be applied to a single subject.")
  par(mfrow=layout)
  for (i in 1:length(samples))
    plot.dmc(samples,subject=i,pll.barplot=TRUE,start=start,end=end,
             main.pll=paste("Subject",i))
}


#' @importFrom graphics par plot lines legend
plot.pp.dmc <- function(pp,style="pdf",layout=NULL,no.layout=FALSE,
  pos="topleft",percentiles=c(10,30,50,70,90),
  dname="Data",mname="Model",model.legend=TRUE,
  show.fit=TRUE,show.fits=TRUE,x.min.max=NA,ylim=NA,
  show.response=NA,
  ltys=NA,lwds=NA)
  # pdf or cdf of data and posterior predictive fits
  # Show response is a vector of integers in 1:n.resp of which resposnes to show,
  # if NA shows all
{

  plot.one <- function(pp,style,data.style,pos,x.min.max=NA,ylim,
    show.fits,dname,mname,model.legend,lgnd,
    n.resp,show.response,ltys,lwds)
  {

    no.facs <- class(try(pp[[data.style]][1][[1]][[2]],silent=TRUE))=="try-error"
    if (no.facs) ni <- 1 else ni <- length(pp[[1]])
    for ( i in 1:ni ) if ( !all(unlist(lapply(pp[[1]][[i]],is.null))) ) {
      if ( no.facs )
        lgnd <- names(lapply(pp[[style]],function(x){out <- attr(x,"cell.name")})) else
          lgnd <- lapply(pp[[style]][[i]],function(x){out <- attr(x,"cell.name")})
        ok.n <- c(1:n.resp)[!unlist(lapply(pp[[1]][[i]],is.null))]
        is.in <- rep(FALSE,n.resp)
        is.in[show.response] <- TRUE
        is.in <- is.in & ((1:n.resp)==ok.n)
        inj <- c(1:n.resp)[is.in]
        if ( style=="pdf" ) {
          if (is.null(pos)) pos <- "topright"
          if (any(is.na(ylim))) ylim <- c(0,max(
            c(unlist(lapply(pp[[style]][[i]],function(x){
              if (length(x$y)==0) NA else max(x$y)})),
              unlist(lapply(pp[[data.style]][[i]],function(x){
                if (length(x$y)==0) NA else max(x$y)}))),na.rm=TRUE))
          xlim <- c(Inf,-Inf)
          for (j in 1:length(ok.n)) {
            if (no.facs)
              xij <- pp[[style]][[ok.n[j]]][[1]]$x else
                xij <- pp[[style]][i][[1]][[ok.n[j]]]$x
              xlim[1] <- ifelse(min(xij)<xlim[1],
                min(xij),xlim[1])
              xlim[2] <- ifelse(max(xij)>xlim[2],
                max(xij),xlim[2])
          }
          if ( !any(is.na(x.min.max)) ) {
            if (x.min.max[1]>xlim[1]) xlim[1] <- x.min.max[1]
            if (x.min.max[2]<xlim[2]) xlim[2] <- x.min.max[2]
          }
          for (j in show.response) {
            if (no.facs) {
              xy <- pp[[style]][[ok.n[j]]][[1]]
              xy.dat <- pp[[data.style]][[ok.n[j]]][[1]]
            } else {
              xy <- pp[[style]][i][[1]][[ok.n[j]]]
              xy.dat <- pp[[data.style]][i][[1]][[ok.n[j]]]
            }
            if ( j == show.response[1] )
              plot(xy,main="",ylim=ylim,xlim=xlim) else
                lines(xy,lty=ltys[j],lwd=lwds[j])
            lines(xy.dat,lty=ltys[j],lwd=lwds[j],col="red")
          }
        } else {
          if (is.null(pos)) pos <- "topleft"
          if (any(is.na(ylim))) ylim <- c(0,1)
          # Get limits
          xlim <- c(Inf,-Inf)
          for ( j in inj ) {
            if (no.facs) ppj <- pp[[data.style]][[ j ]][[1]] else
              ppj <- pp[[data.style]][i][[1]][[ j ]]
            xlim[1] <- ifelse(suppressWarnings(min(ppj))<xlim[1],min(ppj),xlim[1])
            xlim[2] <- ifelse(suppressWarnings(max(ppj))>xlim[2],max(ppj),xlim[2])
            if ( show.fit ) {
              if (no.facs) ppj <- pp[[style]][[ j ]][[1]] else
                ppj <- pp[[style]][i][[1]][[ j ]]
              xlim[1] <- ifelse(suppressWarnings(min(ppj))<xlim[1],min(ppj),xlim[1])
              xlim[2] <- ifelse(suppressWarnings(max(ppj))>xlim[2],max(ppj),xlim[2])
              if (!is.null(attr(pp,"dpqs"))) for ( k in 1:length(attr(pp,"dpqs")) ) {
                if (no.facs) ppk <- attr(pp,"dpqs")[[k]][[style]][[ j ]][[1]] else
                  ppk <- attr(pp,"dpqs")[[k]][[style]][i][[1]][[ j ]]
                xlim[1] <- ifelse(suppressWarnings(min(ppk))<xlim[1],min(ppk),xlim[1])
                xlim[2] <- ifelse(suppressWarnings(max(ppk))>xlim[2],max(ppk),xlim[2])
              }
            }
          }
          if (!any(is.na(x.min.max))) {
            if (x.min.max[1]>xlim[1]) xlim[1] <- x.min.max[1]
            if (x.min.max[2]<xlim[2]) xlim[2] <- x.min.max[2]
          }
          for (j in inj) {
            if (j == inj[1]) plot(NA,NA,main="",type="l",ylim=ylim,xlim=xlim)
            # Draw uncertianty
            if ( show.fits && !is.null(attr(pp,"dpqs")) ) for (k in 1:length(attr(pp,"dpqs"))) {
              ppk <- attr(pp,"dpqs")[[k]]
              if (no.facs) {
                y <- as.numeric(names(ppk[[style]][[ j ]][[1]]))
                x <- ppk[[style]][[ j ]][[1]]
              } else {
                y <- as.numeric(names(ppk[[style]][i][[1]][[ j ]]))
                x <- ppk[[style]][i][[1]][[ j ]]
              }
              lines(x,y,col="grey",lty=ltys[j])
            }
            if (no.facs) ppj <- pp[[data.style]][[ j ]][[1]] else
              ppj <- pp[[data.style]][i][[1]][[ j ]]
            lines(ppj,as.numeric(names(ppj)),col="red",lty=ltys[j],lwd=lwds[j])
            if ( !any(is.na(percentiles)) & !is.null(ppj) )
              points(ppj[percentiles],as.numeric(names(ppj))[percentiles],col="red",cex=1.25)
            if ( show.fit ) {
              if (no.facs) ppj <- pp[[style]][[ j ]][[1]] else
                ppj <- pp[[style]][i][[1]][[ j ]] # ok.n[j]
              lines(ppj,as.numeric(names(ppj)),lty=ltys[j],lwd=lwds[j])
              if ( !any(is.na(percentiles)) & !is.null(ppj) )
                points(ppj[percentiles],as.numeric(names(ppj))[percentiles],pch=16,cex=.75)
            }
          }
        }
        if ( !is.na(pos) ) {
          if ( model.legend ) legend(pos,c(paste(dname,lgnd)[is.in],paste(mname,lgnd)[is.in]),
            lty=c(ltys[show.response],ltys[show.response]),
            lwd=c(lwds[show.response],lwds[show.response]),bty="n",
            col=c(rep("red", length(show.response)),rep("black",length(show.response))) ) else
              legend(pos,paste(dname,lgnd)[is.in],lty=ltys[show.response],lwd=lwds[show.response],bty="n",
                col=rep("red", length(show.response)) )

        }
    }
  }

  if (!show.fit) {
    model.legend <- FALSE
    show.fits <- FALSE
  }

  if ( !all(names(pp)[1:2]==c("pdf","cdf")) ) {
    style <- "cdf"
    dname <- paste("Average",dname)
    mname <- paste("Average",mname)
    pp <- attr(pp,"av")
  }

  if (!any(style %in% c("pdf","cdf")))
    stop("style must be either \"pdf\" or \"cdf\"")
  facs <- dimnames(pp[[1]])
  n.facs <- length(facs)
  resp <- names(pp[[1]][1][[1]])
  if (is.null(resp)) resp <- names(pp[[1]])
  n.resp <- length(resp)
  n.levs <- lapply(facs,length)
  if ( !no.layout ) {
    if ( is.null(layout) ) {
      if (class(try(pp[[data.style]][1][[1]][[2]],silent=TRUE))=="try-error")
        layout <- 1 else layout <- unlist(n.levs)
    }
    if (length(layout)==1) layout <- c(1,layout) else
      layout <- layout[1:2]
    par(mfrow=layout,mar=c(2,2,0,0))
  }
  lgnd <- NULL
  data.style <- paste("data",style,sep=".")
  if (any(is.na(show.response))) show.response <- 1:n.resp
  if (any(is.na(ltys))) ltys <- c(1:n.resp)
  if (any(is.na(lwds))) lwds <- rep(1,n.resp)
  if ( !(length(ltys)==n.resp) ) stop(paste("Must provide",n.resp,"ltys"))
  if (!(length(lwds)==n.resp)) stop(paste("Must provide",n.resp,"lwds"))

  plot.one(pp,style,data.style,pos,x.min.max,ylim,show.fits,
    dname,mname,model.legend,lgnd,n.resp,show.response,ltys,lwds)
}

#' @importFrom graphics hist abline legend
plot.deviance.dmc <- function(ds=NULL,samples=NULL,digits=2,fast=TRUE,
                              main="",xlim=NA)
  ## Posterior deviance histogram
  ## TODO correct Dstats.dmc class problem
{
  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics") else
      ds <- Dstats.ddm(samples, TRUE, fast)
    if (any(is.na(xlim)))
      hist(ds$D,breaks="fd",xlab="Posterior Deviance",main=main) else
        hist(ds$D,breaks="fd",xlab="Posterior Deviance",main=main,xlim=xlim)
    abline(v=ds$Dmean,lty=2)
    abline(v=ds$meanD,lty=3)
    pds <- pd.dmc(ds)
    legend("topright",c(paste("NP =",ds$np),
                        paste("NPmean =",round(pds$Pmean,digits)),
                        paste("NPmin =",round(pds$Pmin,digits)),
                        paste("NPvar =",round(pds$Pvar,digits)),
                        "D(mean theta)","mean(D)"),
           bty="n",lty=c(NA,NA,NA,NA,2:3))
}



### Fixed Effects

#' @importFrom graphics arrows
plotSpar.dmc <- function(est,p.name, a.len=.05)
  # plots ordered cis for p.name, est produced by summary.dmc
{
  lo <- unlist(lapply(est,function(x){x$quantiles[p.name,c("2.5%")]}))
  med <- unlist(lapply(est,function(x){x$quantiles[p.name,c("50%")]}))
  hi <- unlist(lapply(est,function(x){x$quantiles[p.name,c("97.5%")]}))
  ord <- order(med)
  n <- length(med)
  plot(1:n,med[ord],ylim=c(min(lo),max(hi)),ylab=p.name,xlab="Subject",pch=16,
    main="Medians and 95% credible intervals")
  arrows(1:n,med[ord],1:n,lo[ord],angle=90, length=a.len)
  arrows(1:n,med[ord],1:n,hi[ord],angle=90, length=a.len)
}


### Hierachical

#' @importFrom graphics plot
h.profile.dmc <- function(p.name,p.num,min.p,max.p,ps,p.prior,n.point=100,
                          digits=3,ylim=NA)
  # for parameter p.name at position p.num (1 or 2) in p.prior given subject
  # parameters ps draws a likelihood profile and returns the maximum (on a grid
  # of resolution n.point)
{
  if (!(p.name %in% dimnames(ps)[[2]]))
    stop("p.name not in ps")
  if (!(p.num %in% 1:2))
    stop("p.num must be either 1 or 2")
  pps <- seq(min.p,max.p,length.out=n.point)
  ll <- numeric(n.point)
  pop <- p.prior[[p.name]]
  pop$x <- ps[,p.name]
  for (i in 1:n.point)
  {
    pop[[p.num]] <- pps[i]
    ll[i] <- sum(do.call(paste("d",attr(pop,"dist"),sep=""),pop))
  }
  names(ll) <- round(pps,digits)
  p.name <- paste(p.name,names(pop)[p.num],sep=".")
  if (any(is.na(ylim)))
    plot(pps,ll,type="l",xlab=paste(p.name,p.num,sep="."),
         ylab="log-likelihood") else
           plot(pps,ll,type="l",xlab=paste(p.name,p.num,sep="."),
                ylab="log-likelihood",ylim=ylim)
  ll[ll==max(ll)]
}

#' @importFrom stats cor
cor.plausible <- function(hsamples,p.name,cv,plot=FALSE,
                          xlab="r",ylab="Density",main=p.name,
                          fun=NULL,...)
  # Correlations for each interation and chain with of p.name with cv$p.name
{
  if ( !is.null(fun) ) subject.pars <- do.call(rbind,lapply(hsamples,function(x){
    pars <- aperm(x$theta,c(1,3,2));
    pars <- matrix(as.vector(pars),ncol=dim(pars)[3]);
    dimnames(pars)[[2]] <- dimnames(x$theta)[[2]];
    apply(pars,1,fun)
  })) else subject.pars <- do.call(rbind,lapply(hsamples,function(x){
    as.vector(x$theta[,p.name,])}))
  out <- apply(subject.pars,2,function(x){cor(x,cv[[p.name]])})
  if (plot) {
    dns <- density(out)
#    dns$y <- dns$y*diff(range(dns$x))/length(dns$x)
    plot(dns,xlab=xlab,ylab=ylab,main=main,...)
  }
  invisible(out)
}

