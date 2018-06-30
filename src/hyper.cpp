#include <ggdmc.hpp>
using namespace Rcpp;

arma::vec UpdatePriors(arma::mat theta, std::vector<std::string> dists,
  arma::mat p1, arma::mat p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {
  // theta is one of the subject's theta; nchain x npar
  // p1 = usephi[0]; p2 = usephi[1] // nchain x npar
  unsigned int nchain = theta.n_rows;
  double tmp_lp;
  arma::vec out(nchain);
  for (size_t i = 0; i < nchain; i++) {
    tmp_lp = sumlogprior(arma::trans(theta.row(i)), dists,
      arma::trans(p1.row(i)), arma::trans(p2.row(i)), lower, upper, islog);

    if (std::isinf(tmp_lp) && tmp_lp > 0.0 ) {
      out(i) = 1e-10;
    } else {
      out(i) = tmp_lp;
    }
  }

  return out ;
}

//' Extract Start Posterior Sample
//'
//' Extract the theta's of the first MCMC iteration across chains and
//' participants. Note that the ps array in DMC is a nchain x nsubject x
//' nparameter array. Armadillo operates on slice (the third dimension), so
//' chain dimension has to be on slice.
//'
//' @param samples a MCMC sample
//' @return a nsubject x npar x nchain array
//' @examples
//' m1 <- ggdmc::BuildModel(
//'   p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
//'   match.map = list(M=list(s1="r1", s2="r2")),
//'   factors   = list(S=c("s1", "s2")),
//'   constants = c(st0=0, d=0),
//'   responses = c("r1","r2"),
//'   type      = "rd")
//'
//' ## Population distribution
//' pop.prior <- ggdmc::BuildPrior(
//'       dists = rep("tnorm", 6),
//'       p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
//'       p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
//'       lower = c(0,-5, 0, 0, 0, 0),
//'       upper = c(5, 7, 2, 2, 2, 2))
//'
//' dat <- ggdmc139:::h.simulate.dmc(m1, ns=4, n=100, p.prior=pop.prior)
//' dmi <- ggdmc139:::data.model.dmc(dat, m1)
//' ps <- attr(dat, "parameters")
//'
//' p.prior  <- ggdmc::BuildPrior(
//'         dists = rep("tnorm", 6),
//'         p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
//'         p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
//'         lower = c(0,-5, 0, 0, 0, 0),
//'         upper = c(5, 7, 2, 2, 2, 2))
//'
//' ## Make a hyper-prior list
//' mu.prior <- ggdmc::BuildPrior(
//'         dists = rep("tnorm", 6),
//'         p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
//'         p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
//'         lower = c(0,-5, 0, 0, 0, 0),
//'         upper = c(5, 7, 2, 2, 2, 2))
//'
//' sigma.prior <- ggdmc::BuildPrior(
//'           dists = rep("beta", 6),
//'           p1    = c(a=1, v=1, z=1, sz=1, sv=1, t0=1),
//'           p2    = c(1,1,1,1,1,1),
//'           upper = c(2,2,2,2,2,2))
//'
//' pp.prior <- list(mu.prior, sigma.prior)
//'
//'
//' ## Random-effect model
//' ## hs0 <- ggdmc139:::h.samples.dmc(5e2, p.prior, dmi, pp.prior, thin=1)
//' ## hs0 <- ggdmc139:::h.run.dmc(hs0, 1, 1e2, p.migrate=.05, h.p.migrate=.05)
//' theta0 <- ggdmc:::GetTheta0(hs0)
//' @export
// [[Rcpp::export]]
arma::cube GetTheta0(List samples) {
  List samples0 = samples[0];
  unsigned int nsub   = samples.size();
  unsigned int npar   = samples0["n.pars"];
  unsigned int nchain = samples0["n.chains"];
  unsigned int start  = samples0["start"];  // Note start is an R index
  arma::cube out(nsub, npar, nchain);

  for (size_t i = 0; i < nsub; i ++) {
    List subject     = samples[i];
    arma::cube theta = subject["theta"];    // nchain x npar x nmc
    for (size_t j = 0; j < nchain; j++) {
      out.slice(j).row(i) = theta.slice(start - 1).row(j);
    }
  }
  return out; // nsub x npar x nchain
}

arma::cube GetUsethetas(arma::field<arma::mat> usethetas) {
  arma::mat usetheta0 = usethetas(0);   // nchain x npar
  unsigned int nsub = usethetas.n_elem;
  unsigned int npar = usetheta0.n_cols;
  unsigned int nchain = usetheta0.n_rows;
  arma::cube out(nsub, npar, nchain);

  for (size_t i = 0; i < nsub; i++) {
    for (size_t j = 0; j < nchain; j++) {
      out.slice(j).row(i) = usethetas(i).row(j);
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
double sumloghprior(arma::vec location, arma::vec scale,
  std::vector<std::string> ldists, std::vector<std::string> sdists,
  arma::vec lp1, arma::vec sp1, arma::vec lp2, arma::vec sp2, arma::vec llower,
  arma::vec slower, arma::vec lupper, arma::vec supper, arma::uvec llog,
  arma::uvec slog) {

  return sumlogprior(location, ldists, lp1, lp2, llower, lupper, llog) +
  sumlogprior(scale, sdists, sp1, sp2, slower, supper, slog);
}

//' @export
// [[Rcpp::export]]
double sumloghlike(arma::mat thetak, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {
  double out = 0; // thetak: nsub x npar
  for(size_t i = 0; i < thetak.n_rows; i++) {
    out += sumlogprior(arma::trans(thetak.row(i)), dists, p1, p2, lower, upper,
      islog);
  }
  return out;
}


//' @export
// [[Rcpp::export]]
void StartIteration(List samples) {
  List subjecti;
  unsigned int start_C;
  arma::mat theta;
  arma::vec logprior, loglike;

  for (size_t i = 0; i < samples.size(); i++) {
    subjecti = samples[i];
    arma::cube thetai   = subjecti["theta"] ; // nchain x npar x nmc
    arma::mat logpriori = subjecti["summed_log_prior"] ; // nmc x nchain
    arma::mat loglikei  = subjecti["log_likelihoods"] ;  // nmc x nchain
    unsigned int start_R  = subjecti["start"] ; // R index
    start_C  = start_R - 1 ;        // C index
    theta    = thetai.slice(start_C) ; // nchain x npar
    logprior = vectorise(logpriori.row(start_C)) ;   // nchain x 1
    loglike  = vectorise(loglikei.row(start_C)) ;   // nchain x 1
  }
}

void Get_n(List samples, unsigned int& nallpar, unsigned int& ncondition,
            unsigned int& nresponse) {
  List subject0 = samples[0];
  List pprior   = subject0["p.prior"];
  List data     = subject0["data"];
  NumericVector modelAttr = data.attr("model");

  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  std::vector<std::string> dim1 = modelDim[0] ; // row; s1.r1, etc
  std::vector<std::string> dim3 = modelDim[2] ; // row; s1.r1, etc
  ncondition = dim1.size();
  nallpar    = allpar.n_elem;
  nresponse  = dim3.size();

}



