#include <ggdmc.hpp>
#include <chrono>
#include <random>
using namespace Rcpp;

//' Generate a Gamma Vector
//'
//' This is part of DE-MCMC algorithm. This function generates a gamma vector
//' for element-wise computation in Armadillo C++. This function is based on
//' p242 ter Braak (2006) who cited Roberts and Rosenthal (2001)
//'
//' @param npar number of parameters.
//' @param gammamult a tuning parameter stands for for gamma mutation. Default
//' value is 2.38.
//' @return a vector
//' @examples
//' pVec <- c(A = 1.51, b = 2.7, muv1 = 3.32, muv2 = 2.24, t_ND = 0.08,
//'           muw1 = 1.51, muw2 = 3.69, t_delay = 0.31, sv = 1, swt = 0.5)
//' gamma <- GetGamma(length(pVec), 2.38)
//' @export
// [[Rcpp::export]]
arma::vec GetGamma(unsigned int npar, double gammamult = 2.38,
  bool hyper = false) {
  arma::vec out(npar);
  double divisor = hyper ? std::sqrt(4.0*(double)npar) : std::sqrt(2.0*(double)npar);
  for (size_t i = 0; i < npar; i++) {
    out(i) = std::isnan(gammamult) ? R::runif(0.5, 1.0) : gammamult/divisor;
  }
  return out;
}

//' Draw n other chains and shuffle them
//'
//' This is part of DE-MCMC algorithm. \code{PickChains} draws \code{n}
//' chains out of \code{length(chains)} chains, excluding the kth chain.
//' \code{GetSubchains} is used in \code{migration} operator. It draws a subset
//' of chains in \code{nchain} chains.
//'
//' @param k the kth processed chain. Must be an integer within the range of 0
//' to \code{nchain - 1}. No check for errorly using R index.
//' @param n numbers of chain to draw.
//' @param chains an integer vector, indicating chain index, e.g., 0:23
//' @param nchain number of chains. Must be an integer.
//' @return a column vector
//' @keywords PickChains, getsubchains
//' @export
//' @examples
//' chains <- 0:23
//'
//' ## Presuming current processing chain is the 1st chain (C index = 0)
//' ## pick 2 chains out of 24 chains, excluding current chain.
//' PickChains(0, 2, chains)
//'
//' ## Example outputs
//' ##      [,1]
//' ## [1,]   17
//' ## [2,]   12
//' ##      [,1]
//' ## [1,]    2
//' ## [2,]    5
//' ##      [,1]
//' ## [1,]    5
//' ## [2,]    3
//' ##      [,1]
//' ## [1,]   10
//' ## [2,]    8
//' ##      [,1]
//' ## [1,]   15
//' ## [2,]    8
//'
//' ## get a random number of subchains
//' GetSubchains(nchain)
//' ##       [,1]
//' ##  [1,]    0
//' ##  [2,]    3
//' ##  [3,]    5
//' ##  [4,]    9
//' ##  [5,]   10
//' ##  [6,]   12
//' ##  [7,]   14
//' ##  [8,]   15
//' ##  [9,]   18
//' ## [10,]   20
//' ## [11,]   21
//' ## [12,]   22
//'
//' @rdname PickChains
//' @export
// [[Rcpp::export]]
arma::uvec PickChains(unsigned int k, unsigned int n, arma::uvec chains) {
  chains.shed_row(k);
  arma::uvec rchains = arma::shuffle(chains);
  return rchains.rows(0, n - 1);
}

//' @rdname PickChains
//' @export
// [[Rcpp::export]]
arma::uvec GetSubchains(unsigned int nchain, bool debug = false) {
  // Migration algroithm - two-step shuffling
  // Step 1. Select a number l (integer) uniformly between 1 and k to be the
  // number of subpopulations for migration. Note in Turner's DE-MCMC, he
  // treats each chain as a subgroup.
  unsigned int nsubchain;
  if (debug) {
    nsubchain = 1;
  } else {
    nsubchain = (unsigned int)std::ceil((double)nchain * R::runif(0.0, 1.0));
  }

  // Step 2. Generate the chain index (0 to nchain - 1) and shuffle them
  arma::uvec chains  = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  arma::uvec rchains = arma::shuffle(chains);

  // From the shuffled chains, take the first nsubchain out of them.
  arma::uvec subchains = rchains.rows(0, nsubchain - 1);
  return arma::sort(subchains);   // Return shuffled, sorted,  subchains
}

//' @rdname PickChains
//' @export
// [[Rcpp::export]]
arma::uvec SelectEmigrants(unsigned int ngroup, unsigned int k) {

  arma::uvec groups, groups_, rgroups, emigrants;
  // randomly determine how many groups to emigrate
  unsigned int l = (unsigned int)std::ceil((double)(ngroup-1) *
    R::runif(0.0, 1.0));

  // remembering where the current group is in the original order
  groups = arma::linspace<arma::uvec>(0, ngroup - 1, ngroup);
  groups_ = groups;

  groups.shed_row(k); // remove current group
  rgroups = arma::shuffle(groups); // shuffle the rest of group members
  emigrants = rgroups.rows(0, l - 1); // pick the first to the l groups

  // glue the current group back to the selected emigrating groups
  arma::uvec out = arma::join_cols(groups_.row(k), emigrants);
  return arma::sort(out); // sorting them. E.g. 0, 2, 3, 5 within 0, 1,..., 6
}


void MutateDGMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int ngroup, double rp,
  unsigned int npda, double bw, unsigned int ncore,
  unsigned int gpuid) {

  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::mat theta  = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int start, end, k;
  unsigned int m = nchain / ngroup;
  arma::vec tmp;
  arma::uvec subchains;

  for (size_t i = 0; i < ngroup; i++) {

    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t ii = 0; ii < m; ii++) {
      k = subchains(ii);
      tmp = theta.col(k) + R::runif(-rp, rp);
      tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore,
        gpuid, false);
      tmp_logpos = tmp_lp + tmp_ll;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

      cur_logpos = usell(k) + uselp(k);
      if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k) = tmp;
        uselp(k) = tmp_lp;
        usell(k) = tmp_ll;
      }
    }
  }
  usetheta = arma::trans(theta);
}


void MutateDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube thetas,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, unsigned int ngroup,
  double rp) {

  unsigned int k, start, end;
  double tmp_logpos, cur_logpos, tmp_hll, tmp_hlp;
  arma::mat useloc = arma::trans(usephi[0]); //npar x nchain
  arma::mat usesca = arma::trans(usephi[1]);
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  unsigned int m = nchain / ngroup;
  arma::uvec subchains;
  arma::vec tmp_loc, tmp_sca;

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t j = 0; j < m; j++) {
      k = subchains(i);
      cur_logpos = usehlp(k) + usehll(k);

      tmp_loc = useloc.col(k) + R::runif(-rp, rp);
      tmp_sca = usesca.col(k) + R::runif(-rp, rp);
      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
                             sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(thetas.slice(k), pdists, tmp_loc, tmp_sca,
                            plower, pupper, plog);
      tmp_logpos = tmp_hlp + tmp_hll;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        useloc.col(k) = tmp_loc; // npar x nchain
        usesca.col(k) = tmp_sca;
        usehlp(k)  = tmp_hlp;
        usehll(k)  = tmp_hll;
      }
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void MutateDGMCDataChains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists,
  arma::mat usephi0, arma::mat usephi1,
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int ngroup, unsigned int npda = 16384, double bw = .01,
  unsigned int ncore = 1, unsigned int gpuid = 0, double rp = .001) {

  unsigned int start, end, k;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::mat p1_ = arma::trans(usephi0); // npar x nchain
  arma::mat p2_ = arma::trans(usephi1); // npar x nchain
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int m = nchain / ngroup;
  arma::vec tmp;
  arma::uvec subchains;

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t j = 0; j < m; j++) {
      k = subchains(j);
      cur_logpos = usell(k) + uselp(k);
      tmp = theta.col(k) + R::runif(-rp, rp);
      tmp_lp = sumlogprior(tmp, dists, p1_.col(k), p2_.col(k), lower, upper,
                           islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
                          dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
                          npda, bw, ncore, gpuid, false);
      tmp_logpos = tmp_lp + tmp_ll;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k) = tmp;
        uselp(k) = tmp_lp;
        usell(k) = tmp_ll;
      }
    }
  }
  usetheta = arma::trans(theta);
}


void CrossoverDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, unsigned int ngroup,
  double rp, double gammaMult) {

  unsigned int k0, k1, k2, start, end, m;
  double tmp_hll, tmp_hlp, tmp_logpos, cur_logpos;
  arma::uvec dchains, subchains; // direction chains
  arma::vec hgamma, tmp_loc, tmp_sca;
  arma::mat useloc = arma::trans(usephi(0)); //npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  hgamma = GetGamma(npar, gammaMult, true);
  m = nchain / ngroup;

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);
    for (size_t ii = 0; ii < m; ii++) {
      dchains = PickChains(ii, 2, subchains);
      k0 = subchains(ii);
      k1 = dchains(0);
      k2 = dchains(1);
      cur_logpos = usehlp(k0) + usehll(k0);
      tmp_loc = useloc.col(k0) + (hgamma % (useloc.col(k1) - useloc.col(k2))) +
        R::runif(-rp, rp);
      tmp_sca = usesca.col(k0) + (hgamma % (usesca.col(k1) - usesca.col(k2))) +
        R::runif(-rp, rp);
      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
        sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(theta.slice(k0), pdists, tmp_loc, tmp_sca,
        plower, pupper, plog);


      tmp_logpos = tmp_hlp + tmp_hll;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        useloc.col(k0) = tmp_loc; // npar x nchain
        usesca.col(k0) = tmp_sca;
        usehlp(k0)  = tmp_hlp;
        usehll(k0)  = tmp_hll;
      }
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void CrossoverDGMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists,
  arma::mat usephi0, arma::mat usephi1,
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int ngroup,
  unsigned int npda, double bw, unsigned int ncore, unsigned int gpuid,
  double rp, double gammamult) {

  unsigned int k0, k1, k2, start, end, m;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::uvec dchains, subchains;
  arma::vec gamma, tmp;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain matrix
  arma::mat phi0  = arma::trans(usephi0); // npar x nchain
  arma::mat phi1  = arma::trans(usephi1); // npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma  = GetGamma(npar, gammamult);
  m = nchain / ngroup;

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);
    for (size_t j = 0; j < m; j++) {
      dchains = PickChains(j, 2, subchains); // Pick 2 within 0, 1, 2, 3, e.g.
      k0 = subchains(j);
      k1 = dchains(0);
      k2 = dchains(1);
      cur_logpos = uselp(k0) + usell(k0);
      tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) +
        R::runif(-rp, rp);
      tmp_lp = sumlogprior(tmp, dists, phi0.col(k0), phi1.col(k0), lower, upper,
                           islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
                          dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
                          npda, bw, ncore, gpuid, false);
      tmp_logpos = tmp_ll + tmp_lp;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
      }
    }
  }
  usetheta = arma::trans(theta);
}


void CrossoverDMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
  double gammaMult) {

  unsigned int k0, k1, k2;
  double tmp_hll, tmp_hlp, tmp_logpos, cur_logpos;
  arma::uvec chains, subchains;
  arma::vec hgamma, tmp_loc, tmp_sca, noise;
  arma::mat useloc = arma::trans(usephi(0)); // npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  hgamma = GetGamma(npar, gammaMult, true); // true = extras *2 as p1 and p2
  // chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));

  for (size_t i = 0; i < nchain; i++) {
    subchains = PickChains(chains(i), 2, chains); // (b-a) * R::runif(1) + a;
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    // update mu and sigma
    tmp_loc = useloc.col(k0) + (hgamma % (useloc.col(k1) - useloc.col(k2))) +
      noise;
    tmp_sca = usesca.col(k0) + (hgamma % (usesca.col(k1) - usesca.col(k2))) +
      noise;

    /* !!! This is critical !!!*/
    // Update usehll for new theta; nsub x npar x nchain
    usehll(k0) = sumloghlike(theta.slice(k0), pdists, useloc.col(k0),
       usesca.col(k0), plower, pupper, plog);

    cur_logpos = usehlp(k0) + usehll(k0); // Calcualte use.post

    // theta: nsub x npar x nchain == ps: nchain x nsub x npar
    tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
      sp2, llower, slower, lupper, supper, llog, slog);
    tmp_hll = sumloghlike(theta.slice(k0), pdists, tmp_loc, tmp_sca, plower,
      pupper, plog);

    // if (std::isinf(tmp_hlp) && tmp_hlp > 0.0 ) tmp_hlp = 1e-10;
    // if (std::isinf(tmp_hll) && tmp_hll > 0.0 ) tmp_hll = 1e-10;

    tmp_logpos = tmp_hlp + tmp_hll;

    if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
      useloc.col(k0) = tmp_loc; // npar x nchain
      usesca.col(k0) = tmp_sca;
      usehlp(k0) = tmp_hlp;
      usehll(k0) = tmp_hll;
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void CrossoverDMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, std::vector<std::string> dists,
  arma::mat usephi0, arma::mat usephi1, arma::vec lower,
  arma::vec upper, arma::uvec islog, arma::vec allpar,
  std::vector<std::string> parnames, arma::ucube model, std::string type,
  std::vector<std::string> dim1, std::vector<std::string> dim2,
  std::vector<std::string> dim3, arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT, arma::uvec matchcell, arma::uvec isr1,
  unsigned int npda = 16384, double bw = .001, unsigned int ncore = 1,
  unsigned int gpuid = 0, double rp = .001, double gammamult = 2.38,
  bool debug = false) {

  unsigned int k0, k1, k2;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::uvec chains, subchains;
  arma::vec gamma, tmp, noise;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain
  arma::mat phi0  = arma::trans(usephi0);  // npar x nchain
  arma::mat phi1  = arma::trans(usephi1);  // npar x nchain

  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma  = GetGamma(npar, gammamult);
  chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  // chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));

  for (size_t i = 0; i < nchain; i++) {
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    // uselp(k0) = sumlogprior(arma::trans(usetheta.row(k0)), dists, phi0.col(k0),
    //   phi1.col(k0), lower, upper,
    //   islog);

    cur_logpos = uselp(k0) + usell(k0);
    tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) +
      R::runif(-rp, rp);
    tmp_lp = sumlogprior(tmp, dists, phi0.col(k0), phi1.col(k0), lower, upper,
      islog);
    tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1, dim2,
      dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore, gpuid,
      debug);
    // if (std::isinf(tmp_lp) && tmp_lp > 0.0 ) tmp_lp = 1e-10;
    tmp_logpos = tmp_ll + tmp_lp;

    if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

    if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
       theta.col(k0) = tmp;
       uselp(k0)     = tmp_lp;
       usell(k0)     = tmp_ll;
     }
  }
  usetheta = arma::trans(theta);
}

void CrossoverDGMCChains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,  // nchain x 1
  std::vector<std::string> pnames, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog, arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int ngroup,
  unsigned int npda, double bw, unsigned int ncore, unsigned int gpuid,
  double rp, double gammamult) {

  // Rcout << "CrossoverDGMCChains" << std::endl;
  unsigned int k0, k1, k2, start, end, m;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::uvec dchains, subchains; // direction chains
  arma::vec gamma, tmp;

  arma::mat theta  = arma::trans(usetheta); // npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma = GetGamma(npar, gammamult); // .5 * arma::randu(npar) + .5;
  m = nchain / ngroup;

  for (size_t i = 0; i < ngroup; i++) {

    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t ii = 0; ii < m; ii++) {
      dchains = PickChains(ii, 2, subchains); //  dchains are on 0, 1, ... m index
      k0 = subchains(ii);
      k1 = dchains(0);
      k2 = dchains(1);
      cur_logpos = usell(k0) + uselp(k0);
      tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) +
        R::runif(-rp, rp);
      tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore,
        gpuid, false);
      tmp_logpos = tmp_lp + tmp_ll;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
      }
    }

  }
  usetheta = arma::trans(theta);
}


void CrossoverDMCChains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,  // nchain x 1
  std::vector<std::string> pnames,
  std::vector<std::string> dists, arma::vec p1, arma::vec p2,
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  double rp, double gammamult, bool force, unsigned int npda,
  double bw, unsigned int ncore, unsigned int gpuid, bool debug)
{
  // Ter Braak's crossover (2006), expect I randomly select a current chain
  unsigned int k0, k1, k2;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::uvec chains, subchains;
  arma::vec gamma, tmp, noise;
  arma::mat theta = arma::trans(usetheta);   // theta: npar x nchain;
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma = GetGamma(npar, gammamult); // ter Braak's step size
  chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  //chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));

  for (size_t i = 0; i < nchain; i++) {
    // Among the rest of the subpopulation xi \ {x^i_j}, sample two
    // chain/chromosomes, (xi_s and xi_t)
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i); k1 = subchains(0); k2 = subchains(1);
    // noise = 2.0*rp*arma::randu<arma::vec>(npar) + rp;

    // Hu & Tsai and Liang & Wong use norm, but ter Braak does not
    // Uniform(-rp, rp); (b-a) * R::runif(1) + a;
    // tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;
    tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + R::runif(-rp, rp);

    // PDA re-calculation
    if (force) {
      uselp(k0) = sumlogprior(theta.col(k0), dists, p1, p2, lower, upper, islog);
      usell(k0) = sumloglike (theta.col(k0), pnames, allpar, parnames,
        model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
        npda, bw, ncore, gpuid, debug);
    }
    cur_logpos = usell(k0) + uselp(k0);

    tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
    tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1, dim2,
      dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore, gpuid,
      debug);
    tmp_logpos = tmp_lp + tmp_ll;
    if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
      theta.col(k0) = tmp;
      uselp(k0) = tmp_lp;
      usell(k0) = tmp_ll;
    }
  }
  usetheta = arma::trans(theta);
}





// void MigrateDMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
//   std::vector<std::string> pnames,
//   std::vector<std::string> dists, arma::vec p1, arma::vec p2, arma::vec lower,
//   arma::vec upper, arma::uvec islog, arma::vec allpar,
//   std::vector<std::string> parnames, arma::ucube model, std::string type,
//   std::vector<std::string> dim1,
//   std::vector<std::string> dim2,
//   std::vector<std::string> dim3,
//   arma::umat n1idx, arma::uvec ise,
//   arma::umat cellidx, arma::vec RT,
//   arma::uvec matchcell, arma::uvec isr1,
//   double rp, double gammamult, bool force, unsigned int nsim,
//   double bw, unsigned int ncore, unsigned int gpuid, bool debug) {
//   double tmp_logpos, cur_logpos;
//   arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
//   unsigned int npar   = theta.n_rows;
//   unsigned int nchain = theta.n_cols;
//   arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8; could be just 1 chain
//   unsigned int nsubchain = subchains.n_elem;
//
//   arma::mat tmp(npar, nsubchain);
//   arma::vec cur_lp(nsubchain), cur_ll(nsubchain), noise;
//   arma::vec tmp_lp(nsubchain), tmp_ll(nsubchain);
//
//   for(size_t i = 0; i < nsubchain; i++) {
//     noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
//     tmp.col(i) = theta.col(subchains(i)) + noise; // proposal
//     cur_lp(i) = uselp(subchains(i));
//     cur_ll(i) = usell(subchains(i));
//     tmp_lp(i) = sumlogprior(tmp.col(i), dists, p1, p2, lower, upper, islog);
//     tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames,
//       model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT,
//       matchcell, isr1, nsim, bw, ncore, gpuid, debug);
//   }
//
//   // Conduct topological migration within subchains;
//   // starting from the last subchain;
//   // each individual chain/chromosome in subchains is treated as a subgroup
//   tmp_logpos = tmp_ll(nsubchain - 1) + tmp_lp(nsubchain - 1);
//   cur_logpos = cur_ll(0) + cur_lp(0);   // migrate to the first subchain
//   if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
//   if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
//     theta.col(subchains(0)) = tmp.col(nsubchain - 1);
//     uselp(subchains(0)) = tmp_lp(nsubchain - 1);
//     usell(subchains(0)) = tmp_ll(nsubchain - 1);
//   }
//
//   // Continue migration, if nsubchain > 1
//   if (nsubchain != 1) {
//     for(size_t k = 0; k < (nsubchain - 2); k++) {
//       tmp_logpos = tmp_ll(k) + tmp_lp(k);
//       cur_logpos = cur_ll(k + 1) + cur_lp(k + 1);
//       if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
//       if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
//         theta.col(subchains(k + 1))   = tmp.col(k);
//         uselp(subchains(k + 1)) = tmp_lp(k);
//         usell(subchains(k + 1))  = tmp_ll(k);
//       }
//     }
//   }
//   usetheta = arma::trans(theta);
// }

void MigrateDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll,  arma::cube& thetas, // nsub x npar x nchain
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog,
  unsigned int ngroup, double rp) {
  // Rcout << "MigrateDGMC Hyperchains" << std::endl;
  unsigned int start, end, next, k, l;
  double tmp_hlp, tmp_hll, tmp_logpos, cur_hlp, cur_hll, cur_logpos;

  arma::mat useloc = arma::trans(usephi(0)); // nchain x npar to npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  unsigned int m = nchain / ngroup;    // number of chains in each group
  arma::uvec subchains, subgroups;
  arma::vec tmp_loc(npar), tmp_sca(npar);


  for (size_t i = 0; i < ngroup; i++) {
    subgroups = SelectEmigrants(ngroup, i);
    l = subgroups.n_elem;

    for (size_t j = 0; j < l; j++) {
      start = subgroups(j) * m;              // 0, 10, 15, 25
      end   = ((1 + subgroups(j)) * m - 1);  // 4, 14, 19, 29
      subchains = arma::linspace<arma::uvec>(start, end, m);

      for (size_t jj = 0; jj < m; jj++) {
        next = (j == (l-1)) ? subgroups(0)*m + jj : subgroups(j+1) * m + jj;
        k = subchains(jj);

        for (size_t jjj = 0; jjj < npar; jjj++) {
          tmp_loc(jjj) = useloc(jjj, k) + R::rnorm(useloc(jjj, k), rp);
          tmp_sca(jjj) = usesca(jjj, k) + R::rnorm(usesca(jjj, k), rp);
        }

        tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
          sp2, llower, slower, lupper, supper, llog, slog);
        tmp_hll = sumloghlike(thetas.slice(k), pdists, tmp_loc, tmp_sca, plower,
          pupper, plog);
        tmp_logpos = tmp_hlp + tmp_hll;
        cur_logpos = usehlp(next) + usehll(next);

        if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
          useloc.col(next) = tmp_loc;
          usesca.col(next) = tmp_sca;
          usehlp(next)    = tmp_hlp;
          usehll(next)    = tmp_hll;
        }
      }
    }
  }

  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void MigrateDGMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists,
  arma::mat usephi0, arma::mat usephi1, // nchain x npar
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int ngroup, double rp, unsigned int npda = 16384,
  double bw = .01, unsigned int ncore = 1, unsigned int gpuid = 0,
  bool debug = false)
{
  // Rcout << "MigrateDGMC Data chains" << std::endl;
  unsigned int start, end, next, k, l;
  double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos;

  arma::mat phi0 = arma::trans(usephi0);
  arma::mat phi1 = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int m      = nchain / ngroup;
  arma::uvec subchains, subgroups;
  arma::vec tmp(npar);

  for (size_t i = 0; i < ngroup; i++) {
    subgroups = SelectEmigrants(ngroup, i);
    l = subgroups.n_elem;

    for (size_t j = 0; j < l; j++) {
      start = subgroups(j) * m;              // 0, 10, 15, 25
      end   = ((1 + subgroups(j)) * m - 1);  // 4, 14, 19, 29
      subchains = arma::linspace<arma::uvec>(start, end, m);
      for (size_t jj = 0; jj < m; jj++) {
        next  = (j == (l-1)) ? subgroups(0)*m + jj : subgroups(j + 1)*m + jj;
        k = subchains(jj);

        for (size_t jjj = 0; jjj < npar; jjj++) {
          tmp(jjj) = theta(jjj, k) + R::rnorm(theta(jjj, k), rp);
        }

        tmp_lp = sumlogprior(tmp, dists, phi0.col(k), phi1.col(k), lower, upper,
          islog);
        tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore,
          gpuid, debug);
        tmp_logpos = tmp_lp + tmp_ll;
        cur_logpos = uselp(next) + usell(next);
        if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
          theta.col(next) = tmp;
          uselp(next) = tmp_lp;
          usell(next) = tmp_ll;
        }
      }
    }
  }
  usetheta = arma::trans(theta);
}

void MigrateDMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp) {

  arma::mat useloc = arma::trans(usephi(0));
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  arma::uvec subchains = GetSubchains(nchain);
  unsigned int nsubchain = subchains.n_elem;
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain), cur_loc, cur_sca,
     tmp_loc(npar), tmp_sca(npar), noise;

  unsigned int next_chain, k;
  double tmp_hlp, tmp_hll, tmp_logpos, cur_logpos;

  for(size_t i = 0; i < nsubchain; i++) { // 0, 1, 2, 5, ...
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);

    k = subchains(i);

    // cur_loc = useloc.col(next_chain);
    // cur_sca = usesca.col(next_chain);
    tmp_loc = useloc.col(k) + noise;
    tmp_sca = usesca.col(k) + noise;

    tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists,
      lp1, sp1, lp2, sp2, llower, slower, lupper, supper, llog, slog);
    tmp_hll = sumloghlike(theta.slice(k), pdists, tmp_loc, tmp_sca, plower,
      pupper, plog);

    // if (std::isinf(tmp_hlp) && tmp_hlp > 0.0 ) tmp_hlp = 1e-10;
    // if (std::isinf(tmp_hll) && tmp_hll > 0.0 ) tmp_hll = 1e-10;

    tmp_logpos = tmp_hlp + tmp_hll;
    if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

    // Update usehll for new theta; nsub x npar x nchain
    usehll(next_chain) = sumloghlike(theta.slice(next_chain), pdists,
      useloc.col(next_chain), usesca.col(next_chain), plower, pupper, plog);
    cur_logpos = usehlp(next_chain) + usehll(next_chain);


    if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
      useloc.col(next_chain) = tmp_loc;
      usesca.col(next_chain) = tmp_sca;
      usehlp(next_chain) = tmp_hlp;
      usehll(next_chain) = tmp_hll;
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void MigrateDMCHyperchains_old(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll,  arma::cube theta,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp) {

  // Rcout << "MigrateDMC Hyperchains_old" << std::endl;
  arma::mat useloc = arma::trans(usephi(0)); // useloc: npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  arma::uvec subchains = GetSubchains(nchain);
  unsigned int nsubchain = subchains.n_elem;
  unsigned int k;

  double tmp_logpos, cur_logpos;
  arma::mat tmp_loc(npar, nsubchain), tmp_sca(npar, nsubchain);
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain), noise;
  arma::vec tmp_hlp(nsubchain), tmp_hll(nsubchain);

  for(size_t i = 0; i < nsubchain; i++) { // 0, 1, 2, 5, ...
    k = subchains(i);
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    tmp_loc.col(i) = useloc.col(k) + noise; // proposal
    tmp_sca.col(i) = usesca.col(k) + noise; // proposal

    /* !!! This is critical !!!*/
    // Update usehll for new theta; nsub x npar x nchain
    usehll(k) = sumloghlike(theta.slice(k), pdists, useloc.col(k),
      usesca.col(k), plower, pupper, plog);
    cur_hlp(i) = usehlp(k);
    cur_hll(i) = usehll(k);

    tmp_hlp(i) = sumloghprior(tmp_loc.col(i), tmp_sca.col(i), ldists,
      sdists, lp1, sp1, lp2, sp2, llower, slower, lupper, supper, llog, slog);

    // nsub x npar x nchain
    tmp_hll(i) = sumloghlike(theta.slice(k), pdists, tmp_loc.col(i),
      tmp_sca.col(i), plower, pupper, plog);
  }

  tmp_logpos = tmp_hll(nsubchain - 1) + tmp_hlp(nsubchain - 1);
  cur_logpos = cur_hll(0) + cur_hlp(0);   // migrate to the first subchain
  if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
    useloc.col(subchains(0)) = tmp_loc.col(nsubchain - 1);
    usesca.col(subchains(0)) = tmp_sca.col(nsubchain - 1);
    usehlp(subchains(0)) = tmp_hlp(nsubchain - 1);
    usehll(subchains(0)) = tmp_hll(nsubchain - 1);
  }

  // Continue migration, if nsubchain > 1
  if (nsubchain != 1) {
    for(size_t k = 0; k < (nsubchain - 2); k++) {
      tmp_logpos = tmp_hll(k) + tmp_hlp(k);
      cur_logpos = cur_hll(k + 1) + cur_hlp(k + 1);
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        useloc.col(subchains(k + 1))   = tmp_loc.col(k);
        usesca.col(subchains(k + 1))   = tmp_sca.col(k);
        usehlp(subchains(k + 1)) = tmp_hlp(k);
        usehll(subchains(k + 1)) = tmp_hll(k);
      }
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}


void MigrateDMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists,
  arma::mat usephi0, arma::mat usephi1, // nchain x npar
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int nsim = 16384, double bw = .01, unsigned int ncore = 1,
  unsigned int gpuid = 0, bool debug = false, double rp = .001,
  double gammamult = 2.38) {
  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar    = theta.n_rows;
  unsigned int nchain  = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain); // eg. 0, 1, 3, 4, 8
  unsigned int nsubchain = subchains.n_elem;

  unsigned int next_chain, k;
  arma::vec theta_cur, theta_star(npar), noise;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;

  for (size_t i = 0; i < nsubchain; i++) {
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);
    k = subchains(i);

    theta_cur  = theta.col(next_chain);
    theta_star = theta.col(k) + noise;

    tmp_lp = sumlogprior(theta_star, dists, phi0.col(k), phi1.col(k), lower,
      upper, islog);
    tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
      dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw, ncore,
      gpuid, debug);
    if (std::isinf(tmp_lp) && tmp_lp > 0.0) tmp_lp = 1e-10;
    tmp_logpos = tmp_lp + tmp_ll;

    if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

    cur_logpos = uselp(next_chain) + usell(next_chain);

    if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
      theta.col(next_chain) = theta_star;
      uselp(next_chain) = tmp_lp;
      usell(next_chain) = tmp_ll;
    }
  }
  usetheta = arma::trans(theta);
}

void MigrateDMCDatachains_old(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists,
  arma::mat usephi0, arma::mat usephi1, // nchain x npar
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int nsim = 16384, double bw = .01, unsigned int ncore = 1,
  unsigned int gpuid = 0, double rp = .001,
  double gammamult = 2.38) {
  // Rcout << "MigrateDMC Datachains_old" << std::endl;

  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar    = theta.n_rows;
  unsigned int nchain  = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain); // eg. 0, 1, 3, 4, 8
  unsigned int nsubchain = subchains.n_elem;


  // unsigned int next_chain, k;
  arma::vec cur_lp(nsubchain), cur_ll(nsubchain), tmp_lp(nsubchain),
            tmp_ll(nsubchain), noise;
  double tmp_logpos, cur_logpos, tmp_lp_;
  arma::mat tmp(npar, nsubchain);

  for (size_t i = 0; i < nsubchain; i++) {
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    tmp.col(i) = theta.col(subchains(i)) + noise; // proposal
    cur_lp(i) = uselp(subchains(i));
    cur_ll(i) = usell(subchains(i));

    tmp_lp_ = sumlogprior(tmp.col(i), dists, phi0.col(subchains(i)),
      phi1.col(subchains(i)), lower, upper, islog);
    // if (std::isinf(tmp_lp_) && tmp_lp_ > 0.0) tmp_lp_ = 1e-10;

    tmp_lp(i) = tmp_lp_;
    tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames,
      model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT,
      matchcell, isr1, nsim, bw, ncore, gpuid, false);
  }

  tmp_logpos = tmp_ll(nsubchain - 1) + tmp_lp(nsubchain - 1);
  cur_logpos = cur_ll(0) + cur_lp(0);   // migrate to the first subchain
  if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
    theta.col(subchains(0)) = tmp.col(nsubchain - 1);
    uselp(subchains(0)) = tmp_lp(nsubchain - 1);
    usell(subchains(0)) = tmp_ll(nsubchain - 1);
  }

  if (nsubchain != 1) {
    for(size_t k = 0; k < (nsubchain - 2); k++) {
      tmp_logpos = tmp_ll(k) + tmp_lp(k);
      cur_logpos = cur_ll(k + 1) + cur_lp(k + 1);
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(subchains(k + 1)) = tmp.col(k);
        uselp(subchains(k + 1)) = tmp_lp(k);
        usell(subchains(k + 1)) = tmp_ll(k);
      }
    }
  }
  usetheta = arma::trans(theta);
}

void MigrateDGMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int ngroup, double rp, unsigned int npda, double bw,
  unsigned int ncore, unsigned int gpuid) {

  unsigned int start, end, next, k, l;
  double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos;
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int m = nchain / ngroup; // m = 30 / 6 = 5
  arma::uvec subchains, subgroups;
  arma::vec tmp(npar);

  for (size_t i = 0; i < ngroup; i++) {
    subgroups = SelectEmigrants(ngroup, i); // 6 groups: 0, 2, 3, 5
    l = subgroups.n_elem; // l = 4

    for (size_t j = 0; j < l; j++) { // j = 0, 1, 2, 3
      start = subgroups(j) * m;              // 0, 10, 15, 25
      end   = ((1 + subgroups(j)) * m - 1);  // 4, 14, 19, 29
      subchains = arma::linspace<arma::uvec>(start, end, m);

      for (size_t jj = 0; jj < m; jj++) { // jj = 0, 1, 2, 3, 4
        // corresponding index in the next group
        next = (j == (l - 1)) ? subgroups(0)*m + jj : subgroups(j+1)*m + jj;

        k = subchains(jj);
        for (size_t jjj = 0; jjj < npar; jjj++) {
          tmp(jjj) = theta(jjj, k) + R::rnorm(theta(jjj, k), rp);
        }

        tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
        tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT,matchcell, isr1, npda, bw, ncore,
          gpuid, false);
        tmp_logpos = tmp_lp + tmp_ll;
        if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        cur_logpos = uselp(next) + usell(next);
        if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
          theta.col(next) = tmp;
          uselp(next) = tmp_lp;
          usell(next) = tmp_ll;
        }
      }

    }

  }
  usetheta = arma::trans(theta);
}

// void MigrateDMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
//   std::vector<std::string> pnames,
//   std::vector<std::string> dists, arma::vec p1, arma::vec p2, arma::vec lower,
//   arma::vec upper, arma::uvec islog, arma::vec allpar,
//   std::vector<std::string> parnames, arma::ucube model, std::string type,
//   std::vector<std::string> dim1,
//   std::vector<std::string> dim2,
//   std::vector<std::string> dim3,
//   arma::umat n1idx, arma::uvec ise,
//   arma::umat cellidx, arma::vec RT,
//   arma::uvec matchcell, arma::uvec isr1,
//   double rp, double gammamult, bool force, unsigned int nsim,
//   double bw, unsigned int ncore, unsigned int gpuid) {
//   double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos;
//   arma::vec theta_cur;
//
//   arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
//   unsigned int npar   = theta.n_rows;
//   unsigned int nchain = theta.n_cols;
//   arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8;
//   unsigned int nsubchain = subchains.n_elem;   // could be just 1 chain
//   unsigned int next_chain, k;
//   arma::vec theta_star(npar);
//
//   for(size_t i = 0; i < nsubchain; i++) {
//     next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);
//
//     k = subchains(i);
//     theta_cur = theta.col(next_chain);
//     for (size_t j = 0; j < npar; j++) {
//       theta_star(j) = theta(j, k) + R::rnorm(theta(j, k), rp);
//     }
//
//     tmp_lp = sumlogprior(theta_star, dists, p1, p2, lower, upper, islog);
//     tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
//       dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw, ncore,
//       gpuid, false);
//     tmp_logpos = tmp_lp + tmp_ll;
//     if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
//     cur_logpos = uselp(next_chain) + usell(next_chain);
//
//     if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
//       theta.col(next_chain) = theta_star;
//       uselp(next_chain) = tmp_lp;
//       usell(next_chain) = tmp_ll;
//     }
//   }
//   usetheta = arma::trans(theta);
// }

void MigrateDMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists, arma::vec p1, arma::vec p2, arma::vec lower,
  arma::vec upper, arma::uvec islog, arma::vec allpar,
  std::vector<std::string> parnames, arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  double rp, double gammamult, bool force, unsigned int nsim,
  double bw, unsigned int ncore, unsigned int gpuid) {

  double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos;
  arma::vec theta_cur;

  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  // arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8;
  // unsigned int nsubchain = subchains.n_elem;   // could be just 1 chain
  unsigned int next_chain, k, l;
  arma::vec theta_star(npar);
  arma::uvec subgroups;

  for(size_t i = 0; i < nchain; i++) {
    subgroups = SelectEmigrants(nchain, i); // ngroup == nchain
    l = subgroups.n_elem; // l = 4

    for (size_t j = 0; j < l; j++) { // j = 0, 1, 2, 3
        next_chain = (j == (l - 1)) ? subgroups(0) : subgroups(j+1);
        k = subgroups(j);
        for (size_t jj = 0; jj < npar; jj++) {
          theta_star(jj) = theta(jj, k) + R::rnorm(theta(jj, k), rp);
        }

        tmp_lp = sumlogprior(theta_star, dists, p1, p2, lower, upper, islog);
        tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT,matchcell, isr1, nsim, bw, ncore,
          gpuid, false);
        tmp_logpos = tmp_lp + tmp_ll;
        if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        cur_logpos = uselp(next_chain) + usell(next_chain);
        if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
          theta.col(next_chain) = theta_star;
          uselp(next_chain) = tmp_lp;
          usell(next_chain) = tmp_ll;
        }
    }
  }
  usetheta = arma::trans(theta);
}


void MigrateDMCChains_old(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists, arma::vec p1, arma::vec p2, arma::vec lower,
  arma::vec upper, arma::uvec islog, arma::vec allpar,
  std::vector<std::string> parnames, arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  double rp, double gammamult, bool force, unsigned int nsim,
  double bw, unsigned int ncore, unsigned int gpuid) {

  // Rcout << "Migrate old" << std::endl;
  double tmp_logpos, cur_logpos;
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  // arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8; could be just 1 chain
  arma::uvec subchains = GetSubchains(nchain, true);
  unsigned int nsubchain = subchains.n_elem;

  arma::mat tmp(npar, nsubchain);
  arma::vec cur_lp(nsubchain), cur_ll(nsubchain), noise;
  arma::vec tmp_lp(nsubchain), tmp_ll(nsubchain);

  for(size_t i = 0; i < nsubchain; i++) {
    // tmp.col(i) = theta.col(subchains(i)) + arma::randn<arma::vec>(npar) * rp;
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    tmp.col(i) = theta.col(subchains(i)) + noise; // proposal
    // for(size_t ii = 0; ii < nsubchain; ii++) {
    //   tmp(ii, i) = theta(ii, subchains(i)) + R::rnorm(theta(ii, subchains(i)), rp);
    // }

    cur_lp(i) = uselp(subchains(i));
    cur_ll(i) = usell(subchains(i));
    tmp_lp(i) = sumlogprior(tmp.col(i), dists, p1, p2, lower, upper, islog);
    tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames,
      model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT,
      matchcell, isr1, nsim, bw, ncore, gpuid, false);
  }

  // Conduct topological migration within subchains;
  // starting from the last subchain;
  // each individual chain/chromosome in subchains is treated as a subgroup
  tmp_logpos = tmp_ll(nsubchain - 1) + tmp_lp(nsubchain - 1);
  cur_logpos = cur_ll(0) + cur_lp(0);   // migrate to the first subchain
  if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
    theta.col(subchains(0)) = tmp.col(nsubchain - 1);
    uselp(subchains(0)) = tmp_lp(nsubchain - 1);
    usell(subchains(0)) = tmp_ll(nsubchain - 1);
  }

  // Continue migration, if nsubchain > 1
  if (nsubchain != 1) {
    for(size_t k = 0; k < (nsubchain - 2); k++) {
      tmp_logpos = tmp_ll(k) + tmp_lp(k);
      cur_logpos = cur_ll(k + 1) + cur_lp(k + 1);
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(subchains(k + 1))   = tmp.col(k);
        uselp(subchains(k + 1)) = tmp_lp(k);
        usell(subchains(k + 1))  = tmp_ll(k);
      }
    }
  }
  usetheta = arma::trans(theta);
}

//' @export
// [[Rcpp::export]]
List run_dgmc(List samples, arma::uvec force, unsigned int report, double pm,
  double qm, double gammamult, unsigned int ncore, unsigned int ngroup) {

  List samples_in(clone(samples)); // so R' original samples stays
  CheckPnames(samples_in);

  List data          = samples_in["data"];
  double rp          = samples_in["rp"];
  arma::ucube model  = data.attr("model");
  unsigned int npda  = data.attr("n.pda");
  double bw          = data.attr("bw");
  bool debug         = data.attr("debug");
  unsigned int gpuid = data.attr("gpuid");
  arma::vec RT       = data["RT"];
  unsigned int nmc   = samples_in["nmc"];
  unsigned int start_R = samples_in["start"];
  unsigned int nthin = samples_in["thin"];
  unsigned int store_i = start_R - 1;
  unsigned int nsamp = 1 + (nmc - start_R) * nthin;
  arma::cube theta = samples_in["theta"] ; // nchain x npar x nmc
  arma::mat lp     = samples_in["summed_log_prior"] ; // nmc x nchain
  arma::mat ll     = samples_in["log_likelihoods"] ;  // nmc x nchain
  arma::mat usetheta = theta.slice(store_i) ; // nchain x npar
  arma::vec uselp = arma::trans(lp.row(store_i)) ; // nchain x 1
  arma::vec usell = arma::trans(ll.row(store_i)) ;  // nchain x 1

  // extract data-model arguments
  NumericVector modelAttr = data.attr("model");
  std::string type = modelAttr.attr("type");
  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  arma::umat n1idx  = modelAttr.attr("n1.order");
  arma::uvec mc    = modelAttr.attr("match.cell");
  arma::uvec ise   = data.attr("cell.empty") ;
  arma::umat cellidx = cellIdx2Mat(data);
  std::vector<std::string> parnames = modelAttr.attr("par.names");
  std::vector<std::string> dim1 = modelDim[0] ; // row
  std::vector<std::string> dim2 = modelDim[1] ; // col; para
  std::vector<std::string> dim3 = modelDim[2] ; // slice; r1, r2
  arma::uvec isr1 = GetIsR1(modelAttr, type);
  NumericVector pvec = modelAttr.attr("p.vector");
  unsigned int npar = pvec.size();
  std::vector<std::string> pnames = pvec.attr("names");

  std::vector<std::string> pdists(npar);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar);
  arma::uvec plog(npar);
  List pprior = samples_in["p.prior"];
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);

  arma::uvec subgroups, subchains;
  unsigned int nchain = usetheta.n_rows;
  unsigned int start, end, m;
  m = nchain / ngroup;

  for (size_t i = 1; i < nsamp; i++) {

      if (R::runif(0.0, 1.0) < pm) {
        MigrateDGMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
          plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
          n1idx, ise, cellidx, RT, mc, isr1, ngroup, rp, npda, bw,
          ncore, gpuid);
      } else if (R::runif(0.0, 1.0) <= qm) {
        MutateDGMCChains(usetheta, uselp, usell, pnames, pdists, pp1,
          pp2, plower, pupper, plog, allpar, parnames, model, type, dim1, dim2,
          dim3, n1idx, ise, cellidx, RT, mc, isr1, ngroup, rp, npda,
          bw, ncore, gpuid);
      } else {
        CrossoverDGMCChains(usetheta, uselp, usell, pnames, pdists,
          pp1, pp2, plower, pupper, plog, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT, mc, isr1, ngroup,
          npda, bw, ncore, gpuid, rp, gammamult);
      }

    if (i % nthin == 0) {
         store_i++;
         if ((store_i + 1) % report == 0) Rcout << store_i + 1 << " ";
         lp.row(store_i) = uselp.t();  // nmc x nchain
         ll.row(store_i) = usell.t();   // nmc x nchain
         theta.slice(store_i) = usetheta ;
    }
  }
  /* ------------------Output ------------------------------------- */
  samples_in["summed_log_prior"] = lp; // nmc x nchain
  samples_in["log_likelihoods"]  = ll;  // nmc x nchain
  samples_in["theta"]            = theta;    // nchain x npar x nmc
  Rcout << std::endl;
  return samples_in;
}

//' @export
// [[Rcpp::export]]
List run_dmc(List samples, arma::uvec force, unsigned int report, double pm,
  double gammamult, unsigned int ncore, bool debug = false) {

  List samples_in(clone(samples)); // so R' original samples stays
  CheckPnames(samples_in);

  List data          = samples_in["data"];
  double rp          = samples_in["rp"];
  arma::ucube model  = data.attr("model");
  unsigned int npda  = data.attr("n.pda");
  double bw          = data.attr("bw");
  // bool debug         = data.attr("debug");
  unsigned int gpuid = data.attr("gpuid");
  arma::vec RT       = data["RT"];
  unsigned int nmc   = samples_in["nmc"];
  unsigned int start_R = samples_in["start"];
  unsigned int thin = samples_in["thin"];
  unsigned int store_i = start_R - 1;
  unsigned int nsamp = 1 + (nmc - start_R) * thin;

  arma::cube theta = samples_in["theta"];    // nchain x npar x nmc
  arma::mat lp = samples_in["summed_log_prior"];  // nmc x nchain
  arma::mat ll = samples_in["log_likelihoods"];
  arma::mat usetheta = theta.slice(store_i);   // nchain x npar
  arma::vec uselp = arma::trans(lp.row(store_i)); // nchains x 1
  arma::vec usell = arma::trans(ll.row(store_i));

  // extract data-model options
  NumericVector modelAttr = data.attr("model");
  std::string type = modelAttr.attr("type");
  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  arma::umat n1idx = modelAttr.attr("n1.order");
  arma::uvec mc    = modelAttr.attr("match.cell");
  arma::uvec ise   = data.attr("cell.empty") ;
  arma::umat cellidx = cellIdx2Mat(data);
  std::vector<std::string> parnames = modelAttr.attr("par.names");
  std::vector<std::string> dim1 = modelDim[0] ; // row
  std::vector<std::string> dim2 = modelDim[1] ; // col; para
  std::vector<std::string> dim3 = modelDim[2] ; // slice; r1, r2
  arma::uvec isr1 = GetIsR1(modelAttr, type);

  NumericVector pvec = modelAttr.attr("p.vector");
  std::vector<std::string> pnames = pvec.attr("names");
  unsigned int npar = pnames.size();

  List pprior = samples_in["p.prior"];
  std::vector<std::string> pdists(npar);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar);
  arma::uvec plog(npar);
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);

  for (size_t i = 1; i < nsamp; i++) {// From 1 not 0
    if (R::runif(0.0, 1.0) < pm) {

      if (debug) {
        MigrateDMCChains_old(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
          plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
          n1idx, ise, cellidx, RT, mc, isr1, rp, gammamult, force(i), npda, bw,
          ncore, gpuid);
      } else {
        MigrateDMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2, plower,
          pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3, n1idx,
          ise, cellidx, RT, mc, isr1, rp, gammamult, force(i), npda, bw, ncore,
          gpuid);
      }

    } else {
      CrossoverDMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
        plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
        n1idx, ise, cellidx, RT, mc, isr1, rp, gammamult, force(i), npda, bw,
        ncore, gpuid, debug);
    }

    if (i % thin == 0) {
      store_i++;
      if ((store_i + 1) % report == 0) Rcout << store_i + 1 << " ";
      lp.row(store_i) = uselp.t();  // nmc x nchain
      ll.row(store_i) = usell.t();   // nmc x nchain
      theta.slice(store_i)  = usetheta;
    }
  }

  /* ------------------Output ------------------------------------- */
  samples_in["summed_log_prior"] = lp; // nmc x nchain
  samples_in["log_likelihoods"]  = ll;  // nmc x nchain
  samples_in["theta"]            = theta;    // nchain x npar x nmc
  Rcout << std::endl;
  return samples_in;
}

//' @export
// [[Rcpp::export]]
List run_hyper_dmc(List samples, unsigned int report, double pm, double hpm,
  double gammamult, unsigned int ncore, bool debug) {

  List samples_in(clone(samples));
  CheckHyperPnames(samples_in);

  List hyper = samples_in.attr("hyper");
  unsigned int npar = hyper["n.pars"];
  unsigned int nchain = hyper["n.chains"];
  unsigned int nmc  = hyper["nmc"];
  unsigned int thin = hyper["thin"];
  unsigned int start= hyper["start"]; // start_R == 1;
  double rp = hyper["rp"];            // rp is defined in initialise
  unsigned int nsub = samples.size();
  unsigned int start_C = start - 1;   // start_C == 0;
  unsigned int store_i = start_C;    // store_i == 0;
  unsigned int nsamp = 1 + (nmc - start) * thin;

  /* data_hyper/thetas (nsub x npar x nchain) == cps (nchain x nsub x npar) */
  std::vector<std::string> pnames = hyper("p.names");
  List phi      = hyper["phi"];
  arma::mat hlp = hyper["h_summed_log_prior"]; // nmc x nchain
  arma::mat hll = hyper["h_log_likelihoods"];
  arma::cube location = phi[0]; // nchain x npar x nmc
  arma::cube scale    = phi[1];
  arma::cube theta0   = GetTheta0(samples_in); // nsub x npar x nchain
  // Hyper-level data: theta0 == ps (nchain x nsub x npar)
  // extract first nmc thetas from each participants

  arma::field<arma::vec> blocks(npar);  // reserved blocks variable
  arma::field<arma::mat> usephi(2);
  usephi(0) = location.slice(start_C); // nchain x npar
  usephi(1) = scale.slice(start_C);
  arma::vec usehlp = arma::trans(hlp.row(start_C)); // nchain
  arma::vec usehll = arma::trans(hll.row(start_C));

  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];
  List ppprior  = hyper["pp.prior"];     /* Extract pprior & ppprior */
  List lprior   = ppprior[0];
  List sprior   = ppprior[1];
  std::vector<std::string> pdists(npar), ldists(npar), sdists(npar),
     types(nsub);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar), lp1(npar),
     lp2(npar), llower(npar), lupper(npar), sp1(npar), sp2(npar), slower(npar),
     supper(npar), bws(nsub);
  arma::uvec llog(npar), slog(npar), plog(npar), substore(nsub), npdas(nsub),
     gpuids(nsub);
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);
  GetPrior(lprior, ldists, lp1, lp2, llower, lupper, llog);
  GetPrior(sprior, sdists, sp1, sp2, slower, supper, slog);

  // Extract subject level data before entering the loop.
  arma::field<arma::cube> subtheta(nsub); // nchains x npar x nmc; nmc x nchain
  arma::field<arma::mat> usetheta(nsub), lp(nsub), ll(nsub);
  arma::field<arma::vec> uselp(nsub), usell(nsub), allpars(nsub), RTs(nsub);
  arma::field<arma::umat> n1idx(nsub), cellidx(nsub);
  arma::field<arma::uvec> matchcells(nsub), emptycells(nsub), isr1(nsub);
  arma::field<std::vector<std::string>> parnames(nsub), dim1s(nsub), dim2s(nsub),
               dim3s(nsub);
  arma::field<arma::ucube> models(nsub);

  TransformSubjects(samples_in, subtheta, usetheta, lp, uselp, ll, usell,
    substore, types, allpars, n1idx, matchcells, emptycells, cellidx, parnames,
    dim1s, dim2s, dim3s, isr1, models, npdas, bws, gpuids, RTs);

   for (size_t i = 1; i < nsamp; i++) {

       if (R::runif(0, 1) < hpm) {

         if (debug) {
           MigrateDMCHyperchains_old(usephi, usehlp, usehll, theta0, pdists,
             plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
             sdists, sp1, sp2, slower, supper, slog, rp);
         } else {
           MigrateDMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
             plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
             sdists, sp1, sp2, slower, supper, slog, rp);
         }

       } else {
          CrossoverDMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
            plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
            sdists, sp1, sp2, slower, supper, slog, rp, gammamult);
       }

      /* usephi(0) and usephi(1): nchain x npar */
      for (size_t j = 0; j < nsub; j++) { // usethetas(j): nchain x npar

          // Because usephi may change in the CrossoverHyperchains
          // usephi: nchain x npar
          // uselp(j) = UpdatePriors(usetheta(j), pdists, usephi(0), usephi(1),
          //    plower, pupper, plog);

          if (R::runif(0, 1) < pm) {

            if (debug) {
              MigrateDMCDatachains_old(usetheta(j), uselp(j), usell(j), pnames,
                pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
                parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
                n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j), isr1(j));
            } else {
              MigrateDMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
                pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
                parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
                n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j), isr1(j));
            }
          } else {
            CrossoverDMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
              isr1(j));
          }

          if ( i % thin == 0 ) {
            substore(j)++;
            // nmc x nchain; nchain x 1
            lp(j).row(substore(j)) = uselp(j).t();
            ll(j).row(substore(j)) = usell(j).t();  // usetheta: nchain x npar
            subtheta(j).slice(substore(j)) = usetheta(j);
          }
          // theta0s: nsub x npar x nchain == ps: nchain x nsub x npar
          // GetUseTheta0s; theta0s: nsub x npar x nchain
          for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
            theta0.slice(k).row(j) = usetheta(j).row(k);
          }
      }

      // theta0 = GetUsethetas(usethetas);
      if (i % thin == 0) {
        store_i++;
        if ((store_i+1) % report == 0) Rcout << store_i + 1 << " ";
        hlp.row(store_i)  = usehlp.t(); // nmc x nchain = nchain x 1
        hll.row(store_i)  = usehll.t();
        location.slice(store_i) = usephi(0); // nchain x npar x nmc = nchain x npar
        scale.slice(store_i)    = usephi(1);
      }
   }

   /* Reconstruct DMC hypersamples */
   for (size_t ii = 0; ii < nsub; ii++) {
     List subject = samples_in[ii];
     subject["summed_log_prior"] = lp(ii);
     subject["log_likelihoods"]  = ll(ii);
     subject["theta"]  = subtheta(ii);
     samples_in[ii] = subject;
   }

   Rcpp::List newphi      = Rcpp::List::create(
     Rcpp::Named("location")   = location,
     Rcpp::Named("scale")      = scale);

   hyper["h_log_likelihoods"]  = hll;
   hyper["h_summed_log_prior"] = hlp;
   hyper["phi"]                = newphi;

   samples_in.attr("hyper") = hyper;
   return samples_in;
}

//' @export
// [[Rcpp::export]]
List run_hyper_dgmc(List samples, unsigned int report, double pm, double qm,
  double gammamult, unsigned int ngroup, unsigned int ncore) {

  List samples_in(clone(samples));
  CheckHyperPnames(samples_in);

  List hyper = samples_in.attr("hyper");
  unsigned int npar  = hyper["n.pars"];
  unsigned int nmc   = hyper["nmc"];
  unsigned int thin  = hyper["thin"];
  unsigned int startR= hyper["start"]; // start_R == 1;
  double rp = hyper["rp"];            // rp is defined in initialise
  unsigned int nsub = samples.size();
  unsigned int startC = startR - 1;   // start_C == 0;
  unsigned int store_i = startC;      // store_i == 0;
  unsigned int nsamp = 1 + (nmc - startR) * thin;

  /* Hyperparameters, hyper logprior and hyper loglikelihood */
  std::vector<std::string> pnames = hyper("p.names");
  List phi      = hyper["phi"];
  arma::mat hlp = hyper["h_summed_log_prior"]; // nmc x nchain
  arma::mat hll = hyper["h_log_likelihoods"];
  arma::cube location = phi[0]; // nchain x npar x nmc
  arma::cube scale    = phi[1];
  // extract first nmc thetas from each participants; nsub x npar x nchain
  arma::cube theta0  = GetTheta0(samples_in);

  arma::field<arma::vec> blocks(npar); // reserved blocks variable
  arma::field<arma::mat> usephi(2);
  usephi(0) = location.slice(startC); // nchain x npar
  usephi(1) = scale.slice(startC);
  arma::vec usehlp = arma::trans(hlp.row(startC)); // nchain
  arma::vec usehll = arma::trans(hll.row(startC));

  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];
  List ppprior  = hyper["pp.prior"];
  List lprior   = ppprior[0];
  List sprior   = ppprior[1];
  std::vector<std::string> pdists(npar), ldists(npar), sdists(npar), types(nsub);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar), lp1(npar),
     lp2(npar), llower(npar), lupper(npar), sp1(npar), sp2(npar), slower(npar),
     supper(npar), bws(nsub);
  arma::uvec llog(npar), slog(npar), plog(npar), substore(nsub), npdas(nsub),
     gpuids(nsub);
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);
  GetPrior(lprior, ldists, lp1, lp2, llower, lupper, llog);
  GetPrior(sprior, sdists, sp1, sp2, slower, supper, slog);

  // Extract subject level data before entering for loop.
  arma::field<arma::cube> subtheta(nsub); // nchains x npar x nmc
  arma::field<arma::mat> usetheta(nsub), lp(nsub), ll(nsub); // nmc x nchain
  arma::field<arma::vec> uselp(nsub), usell(nsub), allpars(nsub), RTs(nsub);
  arma::field<arma::umat> n1idx(nsub), cellidx(nsub);
  arma::field<arma::uvec> matchcells(nsub), emptycells(nsub), isr1(nsub);
  arma::field<std::vector<std::string>> parnames(nsub), dim1s(nsub),
     dim2s(nsub), dim3s(nsub);
  arma::field<arma::ucube> models(nsub);

  TransformSubjects(samples_in, subtheta, usetheta, lp, uselp, ll, usell,
    substore, types, allpars, n1idx, matchcells, emptycells, cellidx, parnames,
    dim1s, dim2s, dim3s, isr1, models, npdas, bws, gpuids, RTs);

  unsigned int nchain = theta0.n_slices;
  for (size_t i = 1; i < nsamp; i++) {
      if (R::runif(0, 1) < pm) {
           MigrateDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
              plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
              sdists, sp1, sp2, slower, supper, slog, ngroup, rp);
      } else if (R::runif(0, 1) <= qm) {
           MutateDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists, plower,
              pupper, plog, ldists, lp1, lp2, llower, lupper, llog, sdists, sp1,
              sp2, slower, supper, slog, ngroup, rp);
      } else {
           CrossoverDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
              plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
              sdists, sp1, sp2, slower, supper, slog, ngroup, rp, gammamult);
      }

      /* ----------------------------------------------------------
       *  Update data level and hyper data
       * ----------------------------------------------------------*/
      for (size_t j = 0; j < nsub; j++) {
        uselp(j) = UpdatePriors(usetheta(j), pdists, usephi(0), usephi(1),
          plower, pupper, plog); // usephi: nchain x npar

        if (R::runif(0, 1) < pm) {
            MigrateDGMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
              isr1(j), ngroup, rp);
        } else if (R::runif(0, 1) <= qm) {
            MutateDGMCDataChains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j), isr1(j),
              ngroup, npdas(j), bws(j), ncore, 0, rp);
        } else {
            CrossoverDGMCDatachains(usetheta(j), uselp(j), usell(j),
              pnames, pdists, usephi(0), usephi(1), plower, pupper, plog,
              allpars(j), parnames(j), models(j), types[j], dim1s(j), dim2s(j),
              dim3s(j), n1idx(j), emptycells(j), cellidx(j), RTs(j),
              matchcells(j), isr1(j), ngroup, npdas(j), bws(j), ncore, 0, rp,
              gammamult);
        }

        // Update theta0: nsub x npar x nchain
        for (size_t jj = 0; jj < nchain; jj++) { // usetheta is nchain x npar
          theta0.slice(jj).row(j) = usetheta(j).row(jj);
        }
      } // end of subject loop
    // }   // end of group loop

    for (size_t j = 0; j < nsub; j++) {
      if ( i % thin == 0 ) {
        substore(j)++;
        // Rcout << substore(j) << " ";

        lp(j).row(substore(j)) = uselp(j).t();
        ll(j).row(substore(j)) = usell(j).t();
        subtheta(j).slice(substore(j)) = usetheta(j);
      }
    }


    if (i % thin == 0) {
        store_i++;
        if ((store_i+1) % report == 0) Rcout << store_i + 1 << " ";
        hlp.row(store_i)  = usehlp.t(); // nmc x nchain = nchain x 1
        hll.row(store_i)  = usehll.t();
        location.slice(store_i) = usephi(0); // nchain x npar x nmc = nchain x npar
        scale.slice(store_i)    = usephi(1);
    }
  } // end of MCMC iteration

  // Rcout << "Finish MCMC " << std::endl;

  /* ----------------------------------------------------------
   *  Reconstruct DMC hypersamples
   * ----------------------------------------------------------*/
  for (size_t j = 0; j < nsub; j++) {
    List subject = samples_in[j];
    subject["summed_log_prior"] = lp(j);
    subject["log_likelihoods"]  = ll(j);
    subject["theta"]            = subtheta(j);
    samples_in[j]               = subject;
  }

  Rcpp::List newphi      = Rcpp::List::create(
    Rcpp::Named("location")   = location,
    Rcpp::Named("scale")      = scale);
  hyper["h_log_likelihoods"]  = hll;
  hyper["h_summed_log_prior"] = hlp;
  hyper["phi"]                = newphi;

  samples_in.attr("hyper") = hyper;
  return samples_in;
}


void NewMigrateDMCChains_bk(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  std::vector<std::string> dists, arma::vec p1, arma::vec p2, arma::vec lower,
  arma::vec upper, arma::uvec islog, arma::vec allpar,
  std::vector<std::string> parnames, arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  double rp, double gammamult, bool force, unsigned int nsim,
  double bw, unsigned int ncore, unsigned int gpuid, bool debug)
{
  double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos;
  arma::vec theta_cur;

  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8;
  unsigned int nsubchain = subchains.n_elem;   // could be just 1 chain
  unsigned int next_chain, k;
  // Rcout << "New Migration\n";
  arma::vec noise(npar), theta_star(npar);
  // for (size_t j = 0; j < npar; j++) {
  //   noise(j) = R::rnorm();
  // }

  for(size_t i = 0; i < nsubchain; i++) {
    if ((i+1) == nsubchain) {
      next_chain = subchains(0);
    } else {
      next_chain = subchains(i+1);
    }

    k = subchains(i);
    theta_cur = theta.col(next_chain);
    for (size_t j = 0; j < npar; j++) {
      double tmp = theta(j, k);
      theta_star(j) = tmp + R::rnorm(tmp, 1);
    }


    tmp_lp = sumlogprior(theta_star, dists, p1, p2, lower, upper, islog);
    tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
      dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw, ncore,
      gpuid, debug);
    tmp_logpos = tmp_lp + tmp_ll;
    if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

    // Current posterior
    cur_lp = sumlogprior(theta_cur, dists, p1, p2, lower, upper, islog);
    cur_ll = sumloglike(theta_cur, pnames, allpar, parnames, model, type, dim1,
      dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw, ncore,
      gpuid, debug);
    cur_logpos = cur_lp + cur_ll;
    if (std::isnan(cur_logpos)) cur_logpos = -INFINITY;

    // uselp(next_chain) + usell(next_chain);

    if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
      theta.col(next_chain) = theta_star;
      uselp(next_chain) = tmp_lp;
      usell(next_chain) = tmp_ll;
    }
  }


  usetheta = arma::trans(theta);
}

