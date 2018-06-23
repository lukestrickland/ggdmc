# Dynamic Models of Choice with Better Graphic Tools and Quicker Computations 

The package, evolving from dynamic model of choice (_DMC_,
Heathcote et al., 2018), is a generic tool for conducting hierarchical 
Bayesian Computations on cognitive models.  

1. Instead of using Gibbs or HMC, _ggdmc_ uses population-based MCMC (pMCMC) 
samplers. A notable Gibbs example is the Python-based 
HDDM (Wiecki, Sofer & Frank, 2013), which does not allow the user to 
conveniently set the variabilities of DDM parameters. 

2. Differing from DMC (Heathcote, et al., 2018), with only the DE-MCMC 
(Turner, Sederberg, Brown, & Steyvers, 2013) sampler, _ggdmc_ provides a number 
of different pMCMC samplers. It is up to the user to 
decide which sampler works best for their models.  DMC may incorporate these 
pMCMC varieties in the future.  

## Getting Started
Below is an example using the LBA Model (Brown & Heathcote, 2008). See
the tutorials in Heathcote et al., (2018) for more. Note to be more explicit,
the functions in _ggdmc_ usually informs the user what they are doing, such as
_BuildModel_ below. The syntax differs slightly from DMC. Also note that the 
sequence of parameters in a parameter vector (i.e., p.vector) must 
strictly follow the sequence in the _p.vector_ reported by _BuildModel_. 

```
require(ggdmc) 
model <- BuildModel(p.map = list(A = "1", B = "R", t0 = "1",
                            mean_v = c("F", "M"), sd_v = "M", st0 = "1"),
          match.map = list(M = list(s1 = 1, s2 = 2)),
          factors   = list(S = c("s1", "s2"), F = c("f1", "f2")),
          constants = c(sd_v.false = 1, st0 = 0), 
          responses = c("r1", "r2"),
          type      = "norm")

## Population distribution, rate effect on F
pop.mean <- c(A = .4, B.r1 = 1.2, B.r2 = 2.8, t0 = .2,
              mean_v.f1.true = 2.5, mean_v.f2.true = 1.5, mean_v.f1.false = .35,
              mean_v.f2.false = .25, sd_v.true = .25)
pop.scale <- c(A = .1, B.r1 = .1, B.r2 = .1, t0 = .05,
               mean_v.f1.true = .2, mean_v.f2.true = .2, mean_v.f1.false = .2,
               mean_v.f2.false = .2, sd_v.true = .1)
pop.prior <- BuildPrior(
  dists = rep("tnorm", 9),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,   0, 0, .1, NA, NA, NA,  0,  0),
  upper = c(NA, NA, 1, NA, NA, NA, NA, NA, NA))

## Draw parameter prior distributions to visually check
plot(pop.prior)

## Simulate 40 participants, each condition has 250 responses.
## The true parameter vectors are drawn from parameter prior distribution
## specified in 'pop.prior' 
dat <- simulate(model, nsim = 250, nsub = 40, p.prior = pop.prior)
dmi <- BindDataModel(dat, model)    ## dmi = data model instance (thanks to MG)

## A DMC quick way to check data quality. That we want to a good balance of 
## correct and error responses 
## Upper = biased to resposne, 
## Lower = biased away response. 
## First column = greater rate, 
## second lesser = rate
## pdf("figs/data.pdf")
## par(mfcol = c(2, 2))
## for (i in 1:40) {
##   ggdmc:::plot.cell.density(data.cell=dmi[[i]][dmi[[i]]$F=="f1" & dmi[[i]]$S=="s1",],C="r1")
##   ggdmc:::plot.cell.density(data.cell=dmi[[i]][dmi[[i]]$F=="f1" & dmi[[i]]$S=="s2",],C="r2")
##   ggdmc:::plot.cell.density(data.cell=dmi[[i]][dmi[[i]]$F=="f2" & dmi[[i]]$S=="s1",],C="r1")
##   ggdmc:::plot.cell.density(data.cell=dmi[[i]][dmi[[i]]$F=="f2" & dmi[[i]]$S=="s2",],C="r2")
## }
## dev.off()

## Extract the mean and variabilities of parameters across the 40 participants
## ps <- attr(dat, "parameters")
##
## Please use matrixStats package, which is even faster than C functions in 
## base package
##
## round(matrixStats::colMeans2(ps), 2)
## round(matrixStats::colSds(ps), 2)
##    A  B.r1  B.r2  mean_v.f1.true  mean_v.f2.true mean_v.f1.false 
## 0.43  1.22  1.00            0.21            2.48            1.45 
## 0.10  0.10  0.01            0.04            0.17            0.19 
## mean_v.f2.false  sd_v.true     t0
##            0.34       0.36   0.23
##            0.16       0.18   0.10
           
## FIT HDDM
p.prior <- BuildPrior(
  dists = rep("tnorm", 9),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,   0, 0, .1, NA, NA, NA,  0,  0),
  upper = c(NA, NA, 1, NA, NA, NA, NA, NA, NA))

## Specify prior distributions at the hyper level
mu.prior <- BuildPrior(
  dists = rep("tnorm", 9),
  p1    = pop.mean,                           
  p2    = c(1,   1,  1,  2,   2,  2,  2,  1, 1),
  lower = c(0,   0,  0, .1,  NA, NA, NA, NA, 0),
  upper = c(NA, NA, NA, NA,  NA, NA, NA, NA, NA))

## lower and upper are taken care of by defaults.
sigma.prior <- BuildPrior(
  dists = rep("beta", 9),
  p1    = c(A=1, B.r1=1, B.r2=1, t0 = 1, mean_v.f1.true=1, mean_v.f2.true=1,
            mean_v.f1.false=1, mean_v.f2.false=1, sd_v.true = 1),
  p2    = rep(1, 9))

pp.prior <- list(mu.prior, sigma.prior)

## Visually check mu priors and sigma priors
plot(pp.prior[[1]])
plot(pp.prior[[2]])

## Initialise a small sample 
## Get the number of parameters
npar <- length(GetPNames(model)); npar
thin <- 2

## Initiate 512 new hierarchical samples and specify (randomly) only the first 
## iteration.  Others are NA or -Inf. 
## Note the number of chains is 54, a specific feature in population-based 
## Markov chain Monte Carlo
hsam0 <- ggdmc::init_newhier(512, dmi, p.prior, pp.prior, thin = thin, nchain = 54)

## This will take about 1 hr
## pm: probability of migration ; gammamult is a tuning parameter in the 
## DE-MCMC sampler
hsam0 <- ggdmc::run(hsam0, pm = .1)


```

## How to pipe DMC samples to _ggdmc_ samplers 
Diffusion-decisoin model (Ratcliff & McKoon, 2008) is one of the most popular 
cognitive models to fit choice RT data in cognitive psychology.  Here we 
show two examples, one fitting the LBA model and the other fitting DDM 
model.  The speed-up is up to 9 times quicker. 


```
###################
##   LBA model   ##
###################

## DMC could be downloaded at "osf.io/pbwx8".
setwd("~/Documents/DMCpaper/")
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")
setwd("~/Documents/ggdmc_paper/")
## load("data/dmc_pipe_LBA.rda")

model <- model.dmc(p.map = list( A = "1", B = "R", t0 = "1",
                                mean_v = c("F", "M"), sd_v = "M", st0 = "1"),
          match.map = list(M = list(s1 = 1, s2 = 2)),
          factors = list(S = c("s1", "s2"), F = c("f1", "f2")),
          constants = c(sd_v.false = 1, st0 = 0),
          responses = c("r1", "r2"),
          type = "norm")
pop.mean <- c(A=.4, B.r1=.6, B.r2=.8, t0=.3, mean_v.f1.true=1.5, 
              mean_v.f2.true=1, mean_v.f1.false=0, mean_v.f2.false=0,
              sd_v.true = .25)
pop.scale <-c(A=.1, B.r1=.1, B.r2=.1, t0=.05, mean_v.f1.true=.2, 
              mean_v.f2.true=.2, mean_v.f1.false=.2, mean_v.f2.false=.2,
              sd_v.true = .1)
pop.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,0,0,NA,NA,NA,NA,0,.1),
  upper = c(NA,NA,NA,NA,NA,NA,NA,NA,1))
raw.data <- h.simulate.dmc(model, p.prior = pop.prior, n = 250, ns = 40)
data.model <- data.model.dmc(raw.data, model)

ps <- attr(raw.data, "parameters")

p.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1    = pop.mean,
  p2    = pop.scale*5,
  lower = c(0,0,0,NA,NA,NA,NA,0,.1),
  upper = c(NA,NA,NA,NA,NA,NA,NA,NA,NA))
mu.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1    = pop.mean,
  p2    = c(1,1,1,2,2,2,2,1,1),
  lower = c(0,0,0,NA,NA,NA,NA,0,.1),
  upper = c(NA,NA,NA,NA,NA,NA,NA,NA,NA))
sigma.prior <- prior.p.dmc(
  dists = rep("beta", 9),
  p1    = c(A=1, B.r1=1, B.r2=1, t0=1, mean_v.f1.true=1, mean_v.f2.true=1, 
            mean_v.f1.false=1, mean_v.f2.false=1, sd_v.true = 1),
  p2    = rep(1,9))
pp.prior <- list(mu.prior, sigma.prior)

hsamples <- h.samples.dmc(nmc = 512, p.prior, data.model, pp.prior = pp.prior,
  thin = 64)

## Piping 
hsam0 <- ggdmc::run(hsamples, pm = .05)

## Turn off migration. Default pm = 0
hsam1 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior, samples = hsam0, 
  pp.prior = pp.prior, thin = 64))
hsam2 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior, samples = hsam1,
   pp.prior = pp.prior, thin = 64))
hsam3 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior, samples = hsam3,
  pp.prior = pp.prior, thin = 32))

## Check whether MCMC converge
plot(hsam3, hyper = TRUE)
plot(hsam3)

est <- CheckRecovery(hsam3, p.vector = ps, hyper = TRUE)
hest <- summary(hsam3, hyper = TRUE)
est <- summary(hsam3)


###################
##   DDM         ##
###################
rm(list = ls())
setwd("~/Documents/DMCpaper")
source ("dmc/dmc.R")
load_model ("DDM", "ddm.R")
setwd("~/Documents/ggdmc_paper/")
## load("data/hierarchical/dmc_pipe_DDM.rda")
model <- model.dmc(
    p.map     = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "1",
                     t0 = "1", st0 = "1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2"), F = c("f1", "f2")),
  constants = c(st0 = 0, d = 0),
  responses = c("r1", "r2"),
  type      = "rd")
  
## Population distribution
pop.mean <- c(a=2,  v.f1=4, v.f2=3, z=0.5, sz=0.3, sv=1, t0=0.3)
pop.scale <-c(a=0.5,v.f1=.5,v.f2=.5,z=0.1, sz=0.1, sv=.3,t0=0.05)
pop.prior <- prior.p.dmc(
   dists = rep("tnorm", 7),
   p1    = pop.mean,
   p2    = pop.scale,
   lower = c(0,-5, -5, 0, 0, 0, 0),
   upper = c(5, 7,  7, 1, 2, 1, 1) )

raw.data   <- h.simulate.dmc(model, p.prior = pop.prior, n = 32, ns = 8)
data.model <- data.model.dmc(raw.data, model)
ps <- attr(raw.data, "parameters")

p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,
  p2=pop.scale*5,
  lower=c(0,-5, -5, 0, 0, 0, 0),
  upper=c(5, 7,  7, 1, 2, 1, 1)
)

mu.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
  p1=pop.mean,
  p2=pop.scale*5,
  lower=c(0,-5, -5, 0, 0, 0, 0),
  upper=c(5, 7,  7, 1, 2, 1, 1)
)
sigma.prior <- prior.p.dmc(
  dists = rep("beta", length(p.prior)),
  p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),
  p2=c(1,1,1,1,1,1,1),
  upper=c(2,2,2,2,2,2,2)
)

pp.prior <- list(mu.prior, sigma.prior)
  
## Sampling  
hsam0 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior = p.prior, 
      pp.prior = pp.prior, data = data.model, thin = 2), pm = .1)

hsam1 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior = p.prior, 
      pp.prior = pp.prior, samples = hsam0, thin = 32), pm = .1)

hsam2 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior = p.prior, 
      pp.prior = pp.prior, samples = hsam1, thin = 128), pm = .3)

hsam3 <- ggdmc::run(h.samples.dmc(nmc = 512, p.prior = p.prior, 
      pp.prior = pp.prior, samples = hsam2, thin = 2))
      
save(time3, time2, time1, time0, hsam0, hsam1, hsam2, hsam3, data.model,
     pp.prior, ps, p.prior, raw.data, file = "data/dmc_pipe_DDM.rda") 
     
## MCMC converge
plot(hsam3, hyper = TRUE)
plot(hsam3)

est1 <- CheckRecovery(hsam3, p.vector = ps, hyper = TRUE)
hest <- summary(hsam3, hyper = TRUE)
est2 <- summary(hsam3)

tmp <- t(data.frame(lapply(est2, function(x){x[[1]][, 1]})))
round(ps - tmp, 2)


```


## How to conduct automatic convergence checks 
One challenge in Bayesian modeling is to make sure the posterior distribuiton
has converged.  When using the DE-MCMC sampler to fit complex models, some 
regions in a parameter space might be difficult to handle.  Here we provide one 
possible way to conduct automatic convergece checks and repeatedly run 
model fit until a proper posterior distribution has reached.

First, we convert the first stage samples (i.e., hsam0) to a generic object, 
hsam. Then, we use the _repeat_ function to iterate model fit. Meanwhile,
we use _CheckConverged_ to check whether Markov chains are flat, well-mixed,
and have enough effective samples. _flat_ is to check whether 
premature convergence happens. _premature_ convergence is a less
explored subject in MCMC computations in complex cognitive / biological models.
It has been observed in optimization works, using crossover operator 
(Ryan 1996).

```

hsam <- hsam0
counter <- 1

repeat {
  hsam <- ggdmc::run(ggdmc::init_oldhier(512, hsam, .001, thin))
  save(hsam, hsam0, dat, dmi, p.prior, pp.prior, file = "data/tmp.rda")

  converged <- matrix(NA, nrow = length(hsam), ncol = 4)
  for(i in 1:length(hsam)) converged[i,] <- ggdmc::CheckConverged(hsam[[i]])
  counter <- counter + 1
  thin <- thin * 2
  if (all(!converged) || counter > 1e2) {
    break
  }
}


```

## Prerequisites
 - R (>= 3.0.0)
 - Rcpp (>= 0.12.10), RcppArmadillo (>= 0.7.700.0.0), ggplot2 (>= 2.1.0),
   rtdists (>= 0.6-6), gridExtra (>= 2.2-1), ggmcmc (>= 0.7.3), 
   ggthemes (>= 3.0.1), stats (>= 3.2.2), loo (>= 0.1.6), coda (>= 0.16-1)
 - Windows users need Rtools (>= 3.3.0.1959) 
 - OS X users need to install Open MPI library
 - Linux/Unix users may need to install Open MPI library, if it has not 
   been installed. 
 - [Armadillo](https://CRAN.R-project.org/package=RcppArmadillo)
   requires a recent compiler; for the g++ family at least version 4.6.*
   is required. 

## Installing

```
From CRAN: install.packages("ggdmc")
From source: install.packages("ggdmc_0.2.0.0.tar.gz", repos = NULL, type="source")

```

## Citation

If you use this package, please cite the software, for example:

Lin, Y.-S., & Heathcote, A (in preparation). Distributed Genetic Monte Carlo is 
as Effective as Hamiltonian Monte Carlo in Fitting High Dimension Cognitive 
Model. Manuscript in preparation.  Retrieved from https://github.com/TasCL/ggdmc

## Contributors

The documentation, C++ codes, R helper functions and pacakging are developed by 
Yi-Shin Lin. DMC is developed by Andrew Heathcote. Please report bugs to 
[Yi-Shin Lin](mailto:yishin.lin@utas.edu.au). 

## License

GPL-2. Please see License.md/LICENSE for details.

## Acknowledgments

* Early version of density.cpp is based on Voss & Voss's (2012) density.c in 
fast-dm 30.2. 
* Early version of the truncated normal distributions is  based on Jonathan 
Olmsted's <jpolmsted@gmail.com> RcppTN 0.1-8 (https://github.com/olmjo/RcppTN) 
and Christopher Jackson's <chris.jackson@mrc-bsu.cam.ac.uk> R codes in msm package. 
* Armadillo is a collection of C++ library for performing linear
algebra <http://arma.sourceforge.net/>, authored by Conrad Sanderson. 
* Thanks to Matthew Gretton's consulation about rtdists. 
