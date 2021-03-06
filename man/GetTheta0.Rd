% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GetTheta0}
\alias{GetTheta0}
\title{Extract Start Posterior Sample}
\usage{
GetTheta0(samples)
}
\arguments{
\item{samples}{a MCMC sample}
}
\value{
a nsubject x npar x nchain array
}
\description{
Extract the theta's of the first MCMC iteration across chains and
participants. Note that the ps array in DMC is a nchain x nsubject x
nparameter array. Armadillo operates on slice (the third dimension), so
chain dimension has to be on slice.
}
\examples{
m1 <- ggdmc::BuildModel(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
  match.map = list(M=list(s1="r1", s2="r2")),
  factors   = list(S=c("s1", "s2")),
  constants = c(st0=0, d=0),
  responses = c("r1","r2"),
  type      = "rd")

## Population distribution
pop.prior <- ggdmc::BuildPrior(
      dists = rep("tnorm", 6),
      p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
      p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
      lower = c(0,-5, 0, 0, 0, 0),
      upper = c(5, 7, 2, 2, 2, 2))

dat <- ggdmc139:::h.simulate.dmc(m1, ns=4, n=100, p.prior=pop.prior)
dmi <- ggdmc139:::data.model.dmc(dat, m1)
ps <- attr(dat, "parameters")

p.prior  <- ggdmc::BuildPrior(
        dists = rep("tnorm", 6),
        p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
        p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
        lower = c(0,-5, 0, 0, 0, 0),
        upper = c(5, 7, 2, 2, 2, 2))

## Make a hyper-prior list
mu.prior <- ggdmc::BuildPrior(
        dists = rep("tnorm", 6),
        p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
        p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
        lower = c(0,-5, 0, 0, 0, 0),
        upper = c(5, 7, 2, 2, 2, 2))

sigma.prior <- ggdmc::BuildPrior(
          dists = rep("beta", 6),
          p1    = c(a=1, v=1, z=1, sz=1, sv=1, t0=1),
          p2    = c(1,1,1,1,1,1),
          upper = c(2,2,2,2,2,2))

pp.prior <- list(mu.prior, sigma.prior)


## Random-effect model
## hs0 <- ggdmc139:::h.samples.dmc(5e2, p.prior, dmi, pp.prior, thin=1)
## hs0 <- ggdmc139:::h.run.dmc(hs0, 1, 1e2, p.migrate=.05, h.p.migrate=.05)
theta0 <- ggdmc:::GetTheta0(hs0)
}
