# Dynamic Models of Choice with Better Graphic Tools and Quicker Computations 

The _ggdmc_ package, evolving from dynamic model of choice (_DMC_,
Heathcote et al., 2018), is a generic tool for conducting hierarchical 
Bayesian computations (BC) on cognitive process, e.g., decision diffusion 
(Ratcliff, 1978; Ratcliff & McKoon, 2008) and linear ballistic accumultion 
(Brown & Heathcote, 2008), models with the differential evolution Markov 
chain Monte Carlo (DE-MCMC) sampler 
(Turner, Sederberg, Brown, & Steyvers, 2013). In addition to using the fast 
sampler, _ggdmc_ differs from the Python-based HDDM 
(Wiecki, Sofer & Frank, 2013), by allowing the user to specify the 
variabilities of model parameters, e.g., drift rates and the start 
points in decision diffusion model  at also the hierarchical level as well as 
incorporating the features of specifying many different factorial designs 
found often in psychology. 
  
Differring from _DMC_, the _ggdmc_ package uses Rcpp 
(Eddelbuettel & Francois, 2011) and Armadillo C++ (Sanderson & Curtin, 2016) 
to programme time-critical computations and has more different genetic 
operators. 

## Getting Started
Here is a simple example extracted from Andrew Heathcote's DMC workshop 
materials. For further details, please see R help pages in this package. 

```
require(ggdmc) 

model <- ggdmc::BuildModel(p.map = list(A = "1", B = "R", t0 = "1",
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

pop.prior <- ggdmc::BuildPrior(
  dists = rep("tnorm", 9),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,   0, 0, .1, NA, NA, NA,  0,  0),
  upper = c(NA, NA, 1, NA, NA, NA, NA, NA, NA))

## Draw parameter prior distributions to visually check
ggdmc:::plot.prior(pop.prior)

## Simulate 40 participants, each condition has 250 responses.
## The true parameter vectors are drawn from parameter prior distribution
## specified in 'pop.prior' 
dat <- ggdmc:::simulate.model(model, nsim = 250, nsub = 40, p.prior = pop.prior)
dmi <- ggdmc::BindDataModel(dat, model)    ## data model instance

## A quick way to check data quality
## Upper=biased to resposne, Lower = biased away. First column = greater rate, second lesser
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
           
### FIT RANDOM EFFECTS
## Specify prior distributions at the trial/data level
p.prior <- ggdmc::BuildPrior(
  dists = rep("tnorm", 9),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,   0, 0, .1, NA, NA, NA,  0,  0),
  upper = c(NA, NA, 1, NA, NA, NA, NA, NA, NA))

## Specify prior distributions at the hyper level
mu.prior <- ggdmc::BuildPrior(
  dists = rep("tnorm",9),
  p1    = pop.mean,                           
  p2    = c(1, 1, 1, 2, 2, 2, 2, 1, 1),
  lower = c(0,   0,  0, .1,  NA,NA,NA,NA,0),
  upper = c(NA, NA, NA, NA, NA,NA,NA,NA,NA))

sigma.prior <- ggdmc::BuildPrior(
  dists = rep("beta", length(p.prior)),
  p1    = c(A=1, B.r1=1, B.r2=1, t0 = 1, mean_v.f1.true=1, mean_v.f2.true=1,
            mean_v.f1.false=1, mean_v.f2.false=1, sd_v.true = 1),
  p2    = rep(1, 9))

## Visually check mu priors and sigma priors
ggdmc:::plot.prior(mu.prior)
ggdmc:::plot.prior(sigma.prior)
pp.prior <- list(mu.prior, sigma.prior)

## Initialise a small sample 

## Get the number of parameters
## npar <- length(ggdmc::GetPNames(model)); npar
thin <- 2

## Initiate 512 new hierarchical samples and specify (randomly) only the first 
## iteration.  Others are NA or -Inf. 
## Note the number of chains is 54, a specific feature in population-based 
## Markov chain Monte Carlo
hsam0 <- ggdmc::init_newhier(512, dmi, p.prior, pp.prior, thin = thin, nchain = 54)

## This will take about 1 hr
## pm: probability of using migration sampler; gammamult a tuning parameter
## for the DE-MCMC sampler
hsam0 <- ggdmc::run_hyper_dmc(hsam0, report = 1e2, pm = 0.08, gammamult = 2.38,
                              ncore = 1)

```

### Prerequisites
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

### Installing

```
From CRAN: install.packages("ggdmc")
From source: install.packages("ggdmc_0.1.9.7.tar.gz", repos = NULL, type="source")

```

## Citation

If you use this package, please cite the software, for example:

Lin, Y.-S., & Heathcote, A (in preparation). ggdmc: An R package for 
hierarchical Bayesian evidence accumulation models, using differential
evolution Markov Chain Monte Carlo Sampler. Retrieved from
https://github.com/TasCL/ggdmc

## Contributors

The ggdmc C++ codes are developed by Yi-Shin Lin. 

If there is any bug been introduced inadvertently into DMC R codes, they are 
probably errors brought in by the first author. Please report any bugs you may 
find to the first [author](mailto:yishin.lin@utas.edu.au). 

## License

GPL-2. Please see License.md/LICENSE for details.

## Acknowledgments

* density.cpp is based on Voss & Voss's (2012) density.c in fast-dm 30.2. 
* The truncated normal distributions were originally based on Jonathan Olmsted's
<jpolmsted@gmail.com> RcppTN 0.1-8 (https://github.com/olmjo/RcppTN) and 
Christopher Jackson's <chris.jackson@mrc-bsu.cam.ac.uk> R codes in msm. 
* C++ codes in this package depend largely on Dirk Eddelbuettel, Romain 
Francois and Doug Bates's Rcpp, and RcppArmadillo.  
* Armadillo is a collection of C++ library for performing linear
algebra <http://arma.sourceforge.net/>, authored by Conrad Sanderson. 
* Thanks to Matthew Gretton's consulation about rtdists's internal. 
