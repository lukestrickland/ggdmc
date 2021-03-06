% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_plotting.R
\name{plot.cell.density}
\alias{plot.cell.density}
\title{Plot Distributions for Each Cell}
\usage{
\method{plot}{cell.density}(data.cell, C = NA, xlim = c(0, Inf),
  ymax = NA, save.density = FALSE, digits = 2, main = "",
  show.mean = FALSE)
}
\arguments{
\item{data.cell}{a data frame with only onn experimental conditoin}

\item{C}{a correctness column}

\item{xlim}{x censor range}

\item{ymax}{the upper bound for y axis when plotting}

\item{save.density}{whether to save density data}

\item{digits}{print how many digits}

\item{main}{main title for the figure}

\item{show.mean}{whether to show mean}
}
\description{
If !is.na(C) plots density for correct and error responses for a data
frame with columns R (a factor) and RT, adding a boolean score column
for R=C. Otherwise plots each response. Can deal with NA in the RT column,
in which case it provides a summary of p(NA)
}
\examples{
m1 <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
  constants = c(st0=0,d=0),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors   = list(S=c("s1","s2")),
  responses = c("r1","r2"),
  type      = "rd")

p.prior <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
  p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
  lower = c(0,-5, 0, 0, 0, 0),
  upper = c(5, 7, 2, 2, 2, 2))

p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)

dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)
mdi1 <- BindDataModel(dat1, m1)

## Accuracy around 70\%
par(mfrow=c(1,2))
plot_cell_density(data.cell=mdi1[mdi1$S=="s1", ], C="r1", xlim=c(0,2))
plot_cell_density(data.cell=mdi1[mdi1$S=="s2", ], C="r2", xlim=c(0,2))
par(mfrow=c(1,1))
}
