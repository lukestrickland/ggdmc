% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggdmc_utils.R
\name{fac2df}
\alias{fac2df}
\title{Convert factor levels to a data frame}
\usage{
fac2df(model)
}
\arguments{
\item{model}{a model object created by \code{model.dmc}}
}
\value{
a data frame
}
\description{
\code{fac2df} takes a model object created by model.dmc and returns a
data frame with all combination of factor levels.
}
\examples{
m1 <- model.dmc(
      p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
      match.map = list(M=list(s1="r1",s2="r2")),
      factors   = list(S=c("s1","s2"), F=c("f1","f2")),
      constants = c(st0=0,d=0),
      responses = c("r1","r2"),
      type      = "rd")
df1 <- fac2df(m1)
str(df1)
## 'data.frame':	4 obs. of  2 variables:
## $ S: Factor w/ 2 levels "s1","s2": 1 2 1 2
## $ F: Factor w/ 2 levels "f1","f2": 1 1 2 2

m2 <- model.dmc(
            p.map=list(A="1",B="1",v="M",sv="M",t0="1",st0="1"),
            constants=c(st0=0,sv.false=1),
            match.map=list(M=list(s1=1,s2=2)),
            factors=list(S=c("s1","s2")),
            responses=c("r1","r2"),
            type="lnorm")
fac2df(m2)

m3 <- model.dmc(
  p.map     = list(A="1",B="R",t0="1",mean_v=c("F","M"),sd_v="M",st0="1"),
  match.map = list(M=list(s1=1,s2=2)),
  factors   = list(S=c("s1","s2"),F=c("f1","f2")),
  constants = c(sd_v.false=1,st0=0),
  responses = c("r1","r2"),
  type="norm")

fac2df(m3)
##    S  F
## 1 s1 f1
## 2 s2 f1
## 3 s1 f2
## 4 s2 f2
}
\keyword{fac2df}
