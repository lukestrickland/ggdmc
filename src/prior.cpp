#include <ggdmc.hpp>
using namespace Rcpp;

//' Prior Probability Density
//'
//' \code{logprior} is a C++ function. It
//' matches five string types: \code{tnorm}, \code{beta_lu}, \code{gamma_l},
//' \code{lnorm_l}, and \code{constant} to determine which density functions to
//' call (via R API). For truncated normal density, \code{logprior} calls
//' \code{dtn_scalar} (internal dtnorm) to get probability densities from the
//' truncated normal distribution. Whether taking logarithm of the probability
//' density is determined by the boolean \code{log} sent in via
//' \code{prior.p.dmc}. By default, \code{prior.p.dmc} sets \code{log} to 1.
//'
//' @param p.vector the user's supplied parameter vector or a sampler supplied
//' theta/phi vector.
//' @param p.prior a list of list usually created by prior.p.dmc to store the
//' distributional setting for prior parameters.
//' @return a named double vector with probability densities for each model
//' parameter
//' @export
//' @examples
//' ## Use Drift-diffusion model as an example
//' prior1 <- BuildPrior(
//'   dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta"),
//'   p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
//'   p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
//'   lower = c(0,-5, NA, NA, 0, NA),
//'   upper = c(2, 5, NA, NA, 2, NA))
//'
//' ggdmc::view(prior1)
//' ##      mean sd lower upper log    dist  untrans
//' ##   a     1  1     0     2   1   tnorm identity
//' ##   v     0  2    -5     5   1   tnorm identity
//' ##   z     1  1     0     1   1 beta_lu identity
//' ##   sz    1  1  -Inf   Inf   1   tnorm identity
//' ##   sv    1  1     0     2   1 beta_lu identity
//' ##   t0    1  1     0     1   1 beta_lu identity
//'
//' dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta")
//' p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1)
//' p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1)
//' lower = c(0,-5, NA, NA, 0, NA)
//' upper = c(2, 5, NA, NA, 2, NA)
//' islog = rep(1, 6)
//' pVec1 <- c(a=1.15, v=-0.10, z=0.74, sz=1.23, sv=0.11, t0=0.87)
//' pnames <- names(pVec1)
//'
//' \dontrun{
//' setwd("~/Documents/DMC_10052017")
//' source ("dmc/dmc.R")
//' source ("dmc/dmc_myfunction.R")
//' source ("dmc/models/LBA/dists.R")
//' source ("dmc/models/LBA/lba_B.R")
//' }
//'
//' ggdmc::logprior(pVec1, prior1)
//' ggdmc::dprior(pVec1, names(pVec1), dists, p1,p2, lower, upper, islog)
//' \dontrun{ log.prior.dmc(pVec1, prior1) }
//' ##         a          v          z         sz         sv         t0
//' ##-0.5484734 -1.6008386  0.0000000 -0.9453885  0.0000000  0.0000000
//' summed.log.prior(pVec1, prior1)
//'
//'\dontrun{
//' res <- microbenchmark::microbenchmark(
//'   ggdmc::logprior(pVec1, prior1),
//'   ggdmc::dprior(pVec1, pnames, dists, p1,p2, lower, upper, islog),
//'   log.prior.dmc(pVec1, prior1), times = 1e3)
//' }
//'
//' ## Use LBA model as an example
//' prior2 <- ggdmc:::prior.p.dmc(
//'   dists = c("tnorm", "tnorm", "tnorm", "tnorm", "tnorm", "tnorm"),
//'   p1    = c(A=.4, B=.6, mean_v.true=1,  mean_v.false=0,  sd_v.true=.5, t0=.3),
//'   p2    = c(A=.1, B=.1, mean_v.true=.2, mean_v.false=.2, sd_v.true=.1, t0=.05),
//'   lower = c(0,   0, NA, NA, 0, .1),
//'   upper = c(NA, NA, NA, NA, NA, 1))
//'
//' pVec2 <- c(A=0.398, B=0.614, mean_v.true=1.040,
//'   mean_v.false=-0.032, sd_v.true=0.485, t0=0.271)
//'
//' dists = c("tnorm", "tnorm", "tnorm", "tnorm", "tnorm", "tnorm")
//' p1    = c(A=.4, B=.6, mean_v.true=1,  mean_v.false=0,  sd_v.true=.5, t0=.3)
//' p2    = c(A=.1, B=.1, mean_v.true=.2, mean_v.false=.2, sd_v.true=.1, t0=.05)
//' lower = c(0,   0, NA, NA, 0, .1)
//' upper = c(NA, NA, NA, NA, NA, 1)
//' islog = rep(1, 6)
//' pnames <- names(pVec2)
//'
//'
//' ggdmc::logprior(pVec2, prior2)
//' ggdmc::dprior(pVec2, pnames, dists, p1,p2, lower, upper, islog)
//' ##    A       B  mean_v.true  mean_v.false    sd_v.true      t0
//' ## 1.38    1.37         0.67          0.68         1.37    1.91
//' summed.log.prior(pVec2, prior2)
//'
//' \dontrun{
//' log.prior.dmc(pVec2, prior2)
//' res <- microbenchmark::microbenchmark(
//'   ggdmc::logprior(pVec2, prior2),
//'   ggdmc::dprior(pVec2, pnames, dists, p1,p2, lower, upper, islog),
//'   times = 1e3)
//' }
//' ## Unit: microseconds
//' ##                 expr     min      lq      mean   median      uq      max
//' ## ggdmc::logpriordmc    10.966  12.154  13.92842   13.132  15.646   28.776
//' ## ggdmc:::log.prior.dmc  9.081  10.617  12.28825   11.595  14.109   77.875
//' ## ggdmc::dprior         10.128  11.665  13.55588   12.712  15.226   71.170
//' ## log.prior.dmc         96.173 102.774 108.06765  105.533 107.628 2080.882
//' @export
// [[Rcpp::export]]
arma::vec dprior_(arma::vec pvec, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog)
{
  unsigned int npar = dists.size();
  std::string dist1 ("tnorm");  // Available pdf' in DMC
  std::string dist2 ("beta_lu");
  std::string dist3 ("gamma_l");
  std::string dist4 ("lnorm_l");
  std::string dist5 ("constant");
  arma::vec out(npar); out.fill(NA_REAL);
  double x, l, u, tmp;

  for (size_t i = 0; i < npar; i++)
  {
    if ( dists[i].compare(dist1) == 0 ) {         // tnorm
      l = std::isnan(lower[i]) ? -INFINITY : lower[i];
      u = std::isnan(upper[i]) ? INFINITY : upper[i];
      x = pvec[i];
      tmp = dtn_scalar(x, p1[i], p2[i], l, u, islog[i]);
      out[i] = std::isnan(tmp) ? 1e-10 : tmp;
    } else if ( dists[i].compare(dist2) == 0) {  // beta_lu
      // Rcout << "beta_lu " << std::endl;
      l = std::isnan(lower[i]) ? 0 : lower[i];
      u = std::isnan(upper[i]) ? 1 : upper[i];
      x = (pvec[i] - l) / (u - l);
      tmp = islog[i] ? R::dbeta(x, p1[i], p2[i], islog[i]) - std::log(u - l) :
        R::dbeta(x, p1[i], p2[i], islog[i]) / (u - l);
      out[i] = std::isnan(tmp) ? 1e-10 : tmp;

    } else if ( dists[i].compare(dist3) == 0) {  // gamma_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      x = (std::isinf(l) || std::isinf(u)) ? pvec[i] : pvec[i] - l;
      out[i] = R::dgamma(x, p1[i], p2[i], islog[i]);
    } else if ( dists[i].compare(dist4) == 0) {  // lnorm_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      x = (std::isinf(l) || std::isinf(u)) ? pvec[i] : pvec[i] - l;
      out[i] = R::dlnorm(x, p1[i], p2[i], islog[i]);
    } else if (dists[i].compare(dist5) == 0) {  // constant
      out[i] = islog[i] ? 0 : 1;
    } else {
      Rcout << "Distribution type not yet supported" << "\n";
      out[i] = 1e-10;
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
double sumlogprior(arma::vec pvec, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {

  arma::vec den = dprior_(pvec, dists, p1, p2, lower, upper, islog);
  den.replace(arma::datum::inf, 1e-10);
  // den.replace(-arma::datum::inf, 1e-10);  // replace each -INFINITY with 1e-10
  // den.replace(arma::datum::nan, 1e-10);  // replace each nan with 1e-10
  double out = arma::accu(den);
  // if (std::isnan(out)) { out = -INFINITY; }
  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector dprior(NumericVector pvec, List prior) {
  List l1;   // a container to loop through inner list
  std::vector<std::string> pnames = pvec.names() ;
  NumericVector out = NumericVector(pvec.size());

  std::string distType1 ("tnorm");  // Available pdf' in DMC
  std::string distType2 ("beta_lu");
  std::string distType3 ("gamma_l");
  std::string distType4 ("lnorm_l");
  std::string distType5 ("constant");

  bool islog;
  double x, p1, p2, lower, upper, den;
  unsigned int npar = pvec.size();

  for (size_t i = 0; i < npar; i++) {
    l1 = prior[pnames[i]];
    std::string distName = l1.attr("dist");
    p1 = as<double>(l1[0]);  // mean; shape1; shape; meanlog
    p2 = as<double>(l1[1]);  // sd;   shape2; scale; sdlog

    // Do do.call
    if (distName.compare(distType1) == 0) {         // tnorm
      lower = as<double>(l1[2]);
      upper = as<double>(l1[3]);
      islog = as<bool>(l1[4]);
      x = pvec[i];
      den = dtn_scalar(x, p1, p2, lower, upper, islog);
    } else if (distName.compare(distType2) == 0) {  // beta_ul, shape1, shape2
      lower = as<double>(l1[2]);
      upper = as<double>(l1[3]);
      islog = as<bool>(l1[4]);
      x = (pvec[i] - lower) / (upper - lower);
      den = !islog ? R::dbeta(x, p1, p2, 0) / (upper - lower) :
        R::dbeta(x, p1, p2, 1) - std::log(upper - lower);
    } else if (distName.compare(distType3) == 0) {  // gamma_l, shape, scale
      lower = as<double>(l1[2]);
      islog = as<bool>(l1[3]);
      x = pvec[i] - lower;
      den = R::dgamma(x, p1, p2, islog);
    } else if (distName.compare(distType4) == 0) {  // lnorm_l, meanlow, sdlow
      lower = as<double>(l1[2]);
      islog = as<bool>(l1[3]);
      x = pvec[i] - lower;
      den = R::dlnorm(x, p1, p2, islog);
    } else if (distName.compare(distType5) == 0) {  // constant
      islog = Rcpp::as<bool>(l1[1]);
      den = islog ? 0 : 1;
    } else {
      Rcout << "Distribution type not yet supported" << "\n";
      den = 1e-10;
    }
    out[i] = den;
  }

  out.attr("names") = pnames;
  return out;
}

//' @export
// [[Rcpp::export]]
double sumlogpriorNV(arma::vec pvec, List prior) {
  NumericVector pvecNV  = as<NumericVector>(wrap(pvec)) ;
  std::vector<std::string> pnames = prior.names() ;
  pvecNV.names() = pnames;
  return sum(dprior(pvecNV, prior)); // sum is in Rcpp
}

//' @rdname rprior
//' @export
// [[Rcpp::export]]
NumericVector rprior_scalar(List prior)
{
  List l1;   // a list to hold each parameter's setting inside priorList
  unsigned int npar = prior.size();
  NumericVector out = Rcpp::NumericVector(npar);
  std::vector<std::string> pnames = prior.attr("names");

  std::string distType1 ("tnorm");
  std::string distType2 ("beta_lu");
  std::string distType3 ("gamma_l");
  std::string distType4 ("lnorm_l");
  std::string distType5 ("constant");

  double p1, p2, lower, upper;

  for (size_t i = 0; i < npar; i++) {
    l1 = prior[pnames[i]];
    std::string distName = l1.attr("dist");

    p1 = as<double>(l1[0]);  // parameter1: mean; shape1; shape; meanlog
    p2 = as<double>(l1[1]);  // parameter2: sd;   shape2; scale; sdlog
    lower = as<double>(l1[2]);
    upper = as<double>(l1[3]);

    if (distName.compare(distType1) == 0) {         // tnorm
      out[i] = rtn_scalar(p1, p2, lower, upper);
    } else if (distName.compare(distType2) == 0) {  // beta_ul
      out[i] = lower + R::rbeta(p1, p2) * (upper - lower);
    } else if (distName.compare(distType3) == 0) {  // gamma_l
      out[i] = lower + R::rgamma(p1, p2);
    } else if (distName.compare(distType4) == 0) {  // lnorm_l
      out[i] = lower + R::rlnorm(p1, p2);
    } else if (distName.compare(distType5) == 0) {  // constant
      out[i] = p1;
    } else {
      Rcout << "Distribution type not yet supported\n";
      out[i] = 1e-10;
    }
  }

  out.attr("names") = pnames;
  return out;
}

//' @rdname rprior
//' @export
// [[Rcpp::export]]
NumericMatrix rprior_mat(List prior, unsigned int n)
{
  unsigned int nrow = n;
  unsigned int npar = prior.size() ;
  NumericMatrix m = na_matrix(nrow, npar) ;
  NumericVector a;

  for (size_t i = 0; i < nrow; i++) {
    a = rprior_scalar(prior);
    for (size_t j = 0; j < npar; j++) m(i, j) = a[j] ;
  }

  CharacterVector pnames = prior.names();
  colnames(m) = pnames;
  return m ;
}

//' @rdname rprior
//' @export
// [[Rcpp::export]]
arma::vec rprior_vec(std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper)
{
  unsigned int npar = dists.size();
  std::string dist1 ("tnorm"); // Available pdf's Andrew has implemented
  std::string dist2 ("beta_lu");
  std::string dist3 ("gamma_l");
  std::string dist4 ("lnorm_l");
  std::string dist5 ("constant");
  arma::vec out = arma::vec(npar).fill(NA_REAL);
  double l, u;

  // [p1 p2]: [mean sd]; [shape1 shape2]; [shape scale]; [meanlog sdlog]
  for (size_t i = 0; i < npar;  i++) {
    if ( dists[i].compare(dist1) == 0 ) {         // tnorm
      l = std::isnan(lower[i]) ? -INFINITY : lower[i];
      u = std::isnan(upper[i]) ?  INFINITY : upper[i];
      out[i] = rtn_scalar(p1[i], p2[i], l, u);
    } else if ( dists[i].compare(dist2) == 0) {  // beta_ul
      l = std::isnan(lower[i]) ? 0 : lower[i];
      u = std::isnan(upper[i]) ? 1 : upper[i];
      out[i] = l + R::rbeta(p1[i], p2[i]) * (u - l);
    } else if ( dists[i].compare(dist3) == 0) {  // gamma_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      out[i] = R::rgamma(p1[i], p2[i]) + l;
    } else if ( dists[i].compare(dist4) == 0 ) {  // lnorm_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      out[i] = R::rlnorm(p1[i], p2[i]) + l;
    } else {  // constant
      out[i] = p1[i];
    }
  }
  return out;
}

arma::mat rprior(unsigned int n, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper)
{
  unsigned int npar = dists.size();
  arma::mat m = arma::mat(npar, n).zeros();
  for (size_t i = 0; i < n; i++) m.col(i) = rprior_vec(dists, p1, p2, lower, upper);
  return m;
}

