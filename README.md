# Dynamic Models of Choice with Better Graphic Tools and Quicker Computations 

The package, evolving from dynamic model of choice (_DMC_,
Heathcote et al., 2018), is a generic tool for conducting hierarchical 
Bayesian Computations on cognitive models.  

1. Differing from DMC (Heathcote, et al., 2018), using only DE-MCMC 
(Turner, Sederberg, Brown, & Steyvers, 2013), _ggdmc_ provides a number of 
different population-based MCMC (pMCMC) samplers, although DMC may 
incorporate these pCMCM varieties too in the future. It is up to the user to 
decide which sampler works best for their models. 

2. Instead of using Gibbs or HMC, _ggdmc_ uses pMCMC (pMCMC) samplers. A 
notable Gibbs example is the Python-based HDDM (Wiecki, Sofer & Frank, 2013), 
which does not allow the user to conveniently set the variabilities of 
DDM parameters. 

3. In addition to the fast Armadillo C++ implementation, differring from _DMC_, 
the _ggdmc_ allows the user to specify a varity of pMCMC samplers, and
one such sampler is DE-MCMC. 

## Getting Started
Below is an example using the LBA Model (Brown & Heathcote, 2008). See
the tutorials in Heathcote et al., (2018) for more. 

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
## ps <- round( attr(dat, "parameters"), 2)
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
cognitive model to fit choice RT data in cognitive psychology.  Here I 
show an exampling, fitting full DDM model.  The second aim of this example
is to show how to pipe DMC samples to fast DE-MCMC sampler.  The speed-up
is usually 9 to 10 times.

```
## DMC could be downloaded at "osf.io/pbwx8".
setwd("~/Documents/DMCpaper/")
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")
setwd("~/Documents/ggdmc_paper/")
load("data/dmc_pipe_LBA.rda")

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

ps <- round( attr(raw.data, "parameters"), 2)

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
time0 <- system.time(hsam0 <- ggdmc::run(hsamples, pm = .05))

## Turn off migration. Default pm = 0
hsam1 <- h.samples.dmc(nmc = 512, p.prior, samples = hsam0, 
  pp.prior = pp.prior, thin = 64)
time2 <- system.time(hsam1 <- ggdmc::run(hsam1))

```



## Prerequisites
 - R (>= 3.0.0)
 - Rcpp (>= 0.12.10), RcppArmadillo (>= 0.7.700.0.0), ggplot2 (>= 2.1.0),
   rtdists (>= 0.6-6), gridExtra (>= 2.2-1), ggmcmc (>= 0.7.3), 
   ggthemes (>= 3.0.1), stats (>= 3.2.2), loo (>= 0.1.6), coda (>= 0.16-1)
 - Windows users need Rtools (>= 3.3.0.1959), and Microsoft Visual Studio 
   Community (>= 2015) (for Open MPI library and M_PI macro support)
 - OS X users need to install Open MPI library
 - Linux/Unix users may need to install Open MPI library, if it has not 
   been installed. 
 - [Armadillo](https://CRAN.R-project.org/package=RcppArmadillo)
   requires a recent compiler; for the g++ family at least version 4.6.*
   is required. 

Successful cases for Windows OS:
  - Microsoft Visual Studio Community 2015 (Version 14.0.25421.03 Update 3) on  
    Windows 10 64 bits.
  - Microsoft Visual Studio Community 2015 (Version 14.0.24720.1 Update 1), 
    with Rtools 3.4 on Windows 10 64 bits.
  
Unsuccseeful cases for Windows OS:
  - Microsoft Blend for Visual Studio Express 2015   

## Installing

```
From CRAN: install.packages("ggdmc")
From source: install.packages("ggdmc_0.1.9.9.tar.gz", repos = NULL, type="source")

```

## Citation

If you use this package, please cite the software, for example:

Lin, Y.-S., & Heathcote, A (in preparation). Distributed Genetic Monte Carlo is 
as Effective as Hamiltonian Monte Carlo in Fitting High Dimension Cognitive 
Model. Retrieved from https://github.com/TasCL/ggdmc

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
