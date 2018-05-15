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
//' \code{GetSubchains} is used in \code{migration} sampler. It draws a subset of
//' chains in \code{nchain} chains.
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
  return rchains.rows(0, n-1);
}

//' @rdname PickChains
//' @export
// [[Rcpp::export]]
arma::uvec GetSubchains(unsigned int nchain) {
  // Migration algroithm - two-step shuffling
  // Step 1. Select a number l (integer) uniformly between 1 and k to be the
  // number of subpopulations for migration. Note in Turner's DE-MCMC, he
  // treats each chain as a subgroup.
  unsigned int nsubchain = (unsigned int)std::ceil((double)nchain * R::runif(0.0, 1.0));

  // Step 2. Generate the chain index (0 to nchain - 1) and shuffle them
  arma::uvec chains  = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  arma::uvec rchains = arma::shuffle(chains);

  // From the shuffled chains, take the first nsubchain out of them.
  arma::uvec subchains = rchains.rows(0, nsubchain - 1);
  return arma::sort(subchains);   // Return shuffled, sorted,  subchains
}

//' @export
// [[Rcpp::export]]
arma::uvec SelectEmigrants(unsigned int ngroup, unsigned int k) {
  arma::uvec groups, groups_, rgroups, emigrants;
  // randomly determine how many groups to emigrate
  unsigned int l = (unsigned int)std::ceil((double)(ngroup-1) * R::runif(0.0, 1.0));
  groups = arma::linspace<arma::uvec>(0, ngroup - 1, ngroup);
  // remembering where the current group is in the original grouping order
  groups_ = arma::linspace<arma::uvec>(0, ngroup - 1, ngroup);
  groups.shed_row(k); // remove current group
  rgroups = arma::shuffle(groups); // shuffle the rest of groups
  emigrants = rgroups.rows(0, l - 1); // pick the first to the l groups
  // glue the current group back to the selected emigrating groups
  arma::uvec out = arma::join_cols(groups_.row(k), emigrants);
  return arma::sort(out); // sorting them. E.g. 0, 2, 3, 5 within 0, 1,..., 6
}


void MigrateDGMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, arma::uvec subgruops,
  unsigned int ngroup, double rp, bool force, unsigned int npda, double bw,
  unsigned int ncore, unsigned int gpuid, bool debug) {

  unsigned int npar, nchain, m, start, end, next, k0, nsubgroup;
  double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos;
  arma::uvec subchains, subgroups;
  arma::vec tmp;
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  npar   = theta.n_rows;
  nchain = theta.n_cols;
  m = nchain / ngroup;  // how many members in a subgroup
  nsubgroup = subgroups.n_elem;

  for (size_t i = 0; i < nsubgroup; i++) {
    // subgroups = SelectEmigrants(ngroup, i); // emigrating groups + current group
    // for (size_t j = 0; j < subgroups.n_elem; j++) { // loop only a subset of groups
      start = i*m;
      end   = ((1 + i)*m - 1);
      subchains = arma::linspace<arma::uvec>(start, end, m); // overall index
      for (size_t j = 0; j < m; j++) { // loop group members
        next = (i == subgroups.n_elem-1) ? subgroups(0)*m + j : subgroups(i+1)*m + j;
        k0 = subchains(j); // so does k0
        tmp = theta.col(k0) + R::runif(-rp, rp);
        tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
        tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT,matchcell, isr1, npda, bw, ncore,
          gpuid, debug);
        tmp_logpos = tmp_lp + tmp_ll;
        if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        // Reference log-posterior; Same location, but in the next group
        // theta_cur = theta.col(next);
        cur_logpos = uselp(next) + usell(next);
        if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
          theta.col(next) = tmp;
          uselp(next) = tmp_lp;
          usell(next) = tmp_ll;
        }
      }
    //}
  }
  usetheta = arma::trans(theta);
}

void MutateDGMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int gropuidx,
  unsigned int ngroup, double rp, bool force, unsigned int npda, double bw,
  unsigned int ncore,
  unsigned int gpuid, bool debug) {

  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::mat theta  = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int start, end, k0;
  unsigned int m = nchain / ngroup;
  arma::vec tmp;
  arma::uvec subchains;

  // for (size_t i = 0; i < ngroup; i++) {
    start = gropuidx*m;
    end   = ((1 + gropuidx)*m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m); // overall index
    for (size_t i = 0; i < m; i++) {
      k0 = subchains(i);
      tmp = theta.col(k0) + R::runif(-rp, rp);
      tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore,
        gpuid, debug);
      tmp_logpos = tmp_lp + tmp_ll;
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

      cur_logpos = usell(k0) + uselp(k0);
      if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
      }
    }
  // }
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
  arma::uvec matchcell, arma::uvec isr1, unsigned int groupidx,
  unsigned int ngroup, double rp, bool force, unsigned int npda, double bw,
  unsigned int ncore, unsigned int gpuid, double gammamult, bool debug) {

  unsigned int k0, k1, k2, start, end, m;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::uvec subchains, idx, dchains; // direction chains
  arma::vec gamma, tmp;
  arma::mat theta  = arma::trans(usetheta); // npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma = GetGamma(npar, gammamult); // .5 * arma::randu(npar) + .5;
  // chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  m = nchain / ngroup;
  subchains = arma::linspace<arma::uvec>(0, m - 1, m); // sub index

  // for (size_t i = 0; i < ngroup; i++) {
    start = groupidx*m;
    end   = ((1 + groupidx)*m - 1);
    // a subset of overall index
    idx = arma::linspace<arma::uvec>(start, end, m);
    for (size_t i = 0; i < m; i++) { // j is sub-index
      dchains = PickChains(i, 2, subchains); //  dchains are on 0, 1, ... m index
      // k0, k1, & k2 return to overall index via idx;
      k0 = idx(i);
      k1 = idx(dchains(0));
      k2 = idx(dchains(1));
      tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + R::runif(-rp, rp);

      if (force) {     // PDA re-calculation
        uselp(k0) = sumlogprior(theta.col(k0), dists, p1, p2, lower, upper, islog);
        usell(k0) = sumloglike (theta.col(k0), pnames, allpar, parnames, model,
          type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
          npda, bw, ncore, gpuid, debug);
      }

      cur_logpos = usell(k0) + uselp(k0);

      tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore,
        gpuid, debug);
      tmp_logpos = tmp_lp + tmp_ll;

      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0.0, 1.0) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
      }
    }
  // }
  usetheta = arma::trans(theta);
}


void MutateDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube thetas,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog,
  unsigned int ngroup, double rp) {

  double tmp_logpos, cur_logpos, tmp_hll, tmp_hlp;
  arma::mat useloc = arma::trans(usephi[0]); //npar x nchain
  arma::mat usesca = arma::trans(usephi[1]);
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  unsigned int start, end, k0;
  unsigned int m = nchain / ngroup;
  // arma::uvec chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  // arma::uvec subchains = arma::linspace<arma::uvec>(0, m - 1, m); // sub index
  arma::vec noise, tmp_loc, tmp_sca;
  arma::uvec subchains;

  for (size_t i = 0; i < ngroup; i++) {
    start = i*m;
    end   = ( (1 + i)*m - 1);
    subchains   = arma::linspace<arma::uvec>(start, end, m); // overall index
    for (size_t j = 0; j < m; j++) {
      noise = 2.0*rp*arma::randu<arma::vec>(npar) + rp;
      k0 = subchains(j);
      tmp_loc = useloc.col(k0) + noise;
      tmp_sca = usesca.col(k0) + noise;
      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
        sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(thetas.slice(k0), pdists, tmp_loc, tmp_sca,
        plower, pupper, plog);
      tmp_logpos = tmp_hlp + tmp_hll;

      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      cur_logpos = usehlp(k0) + usehll(k0);
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

  unsigned int npar, nchain, m, start, end, next, k0;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos;
  arma::mat p1_ = arma::trans(usephi0); // npar x nchain
  arma::mat p2_ = arma::trans(usephi1); // npar x nchain
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  npar   = theta.n_rows;
  nchain = theta.n_cols;
  m = nchain / ngroup;
  arma::vec noise, tmp;
  arma::uvec subchains;

  for (size_t i = 0; i < ngroup; i++) {
    start = i*m;
    end   = ( (1 + i)*m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m); // overall index
    for (size_t j = 0; j < m; j++) {
      k0 = subchains(j);
      noise = (2.0*rp*arma::randu<arma::vec>(npar) + rp);
      tmp = theta.col(k0) + noise;
      tmp_lp = sumlogprior(tmp, dists, p1_.col(k0), p2_.col(k0), lower, upper,
        islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore,
        gpuid, false);
      tmp_logpos = tmp_lp + tmp_ll;

      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      cur_logpos = usell(k0) + uselp(k0);
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
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
  arma::uvec chains, subchains, dchains; // direction chains
  arma::vec hgamma, tmp_loc, tmp_sca, noise;
  arma::mat useloc = arma::trans(usephi[0]); //npar x nchain
  arma::mat usesca = arma::trans(usephi[1]);
  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  hgamma = GetGamma(npar, gammaMult, true);
  // arma::uvec chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  m = nchain / ngroup;
  subchains = arma::linspace<arma::uvec>(0, m - 1, m); // sub index

  for (size_t i = 0; i < ngroup; i++) { // loop through all groups
    start = i*m;
    end   = ((1+i)*m - 1);
    chains   = arma::linspace<arma::uvec>(start, end, m); // overall index
    for (size_t j = 0; j < m; j++) {
      // noise = 2.0*rp*arma::randu<arma::vec>(npar) + rp;
      dchains = PickChains(j, 2, subchains); // sub index
      k0 = chains(j); k1 = chains(dchains(0)); k2 = chains(dchains(1));

      tmp_loc = useloc.col(k0) + (hgamma % (useloc.col(k1) - useloc.col(k2))) +
        R::runif(-rp, rp);
      tmp_sca = usesca.col(k0) + (hgamma % (usesca.col(k1) - usesca.col(k2))) +
        R::runif(-rp, rp);
      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
        sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(theta.slice(k0), pdists, tmp_loc, tmp_sca,
        plower, pupper, plog);
      tmp_logpos = tmp_hlp + tmp_hll;
      cur_logpos = usehlp(k0) + usehll(k0);

      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0,1) < std::exp(tmp_logpos - cur_logpos)) {
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
  arma::uvec chains, subchains, dchains;
  arma::vec gamma, tmp, noise;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain matrix
  arma::mat phi0  = arma::trans(usephi0); // npar x nchain
  arma::mat phi1  = arma::trans(usephi1); // npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma  = GetGamma(npar, gammamult);
  m = nchain / ngroup;
  subchains = arma::linspace<arma::uvec>(0, m - 1, m);

  for (size_t i = 0; i < nchain; i++) {
    start = i*m;
    end   = ( (1 + i)*m -1);
    chains= arma::linspace<arma::uvec>(start, end, m); // subindex, e.g., 4,5,6,7
    for (size_t j = 0; j < m; j++) {
      dchains = PickChains(j, 2, subchains); // Pick 2 within 0, 1, 2, 3, e.g.
      k0 = chains(j); k1 = chains(dchains(0)); k2 = chains(dchains(1));
      tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + R::runif(-rp, rp);

      tmp_lp = sumlogprior(tmp, dists, phi0.col(k0), phi1.col(k0), lower,
        upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model,
        type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
        npda, bw, ncore, gpuid, false);
      tmp_logpos = tmp_ll + tmp_lp;
      cur_logpos = uselp(k0) + usell(k0);

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
  chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  // arma::uvec chains   = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));

  for (size_t i = 0; i < nchain; i++) {
    subchains = PickChains(chains(i), 2, chains); // (b-a) * R::runif(1) + a;
    k0 = chains(i); k1 = subchains(0);  k2 = subchains(1);

    // noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    // update mu and sigma
    tmp_loc = useloc.col(k0) + (hgamma % (useloc.col(k1) - useloc.col(k2))) +
      R::runif(-rp, rp);
    tmp_sca = usesca.col(k0) + (hgamma % (usesca.col(k1) - usesca.col(k2))) +
      R::runif(-rp, rp);

    // Update use.loglike for new theta; nsub x npar x nchain
    usehll(k0) = sumloghlike(theta.slice(k0), pdists, useloc.col(k0),
      usesca.col(k0), plower, pupper, plog);
    cur_logpos = usehlp(k0) + usehll(k0); // Calcualte use.post

    // theta: nsub x npar x nchain == ps: nchain x nsub x npar
    tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
      sp2, llower, slower, lupper, supper, llog, slog);
    tmp_hll = sumloghlike(theta.slice(k0), pdists, tmp_loc, tmp_sca, plower,
      pupper, plog);
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
    subchains = PickChains(chains(i), 2, chains); // exclude chain i and pick two other chains
    k0 = chains(i);  k1 = subchains(0); k2 = subchains(1);
    // noise = 2.0*rp*arma::randu<arma::vec>(npar) + rp;
    tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + R::runif(-rp, rp);

    cur_logpos = uselp(k0) + usell(k0);
    tmp_lp = sumlogprior(tmp, dists, phi0.col(k0), phi1.col(k0), lower, upper, islog);
    tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1, dim2,
      dim3, n1idx, ise, cellidx, RT, matchcell, isr1, npda, bw, ncore, gpuid, debug);
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
  double bw, unsigned int ncore, unsigned int gpuid, bool debug) {
  double tmp_logpos, cur_logpos;
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8; could be just 1 chain
  unsigned int nsubchain = subchains.n_elem;

  arma::mat tmp(npar, nsubchain);
  arma::vec cur_lp(nsubchain), cur_ll(nsubchain), noise;
  arma::vec tmp_lp(nsubchain), tmp_ll(nsubchain);

  for(size_t i = 0; i < nsubchain; i++) {
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    tmp.col(i) = theta.col(subchains(i)) + noise; // proposal
    cur_lp(i) = uselp(subchains(i));
    cur_ll(i) = usell(subchains(i));
    tmp_lp(i) = sumlogprior(tmp.col(i), dists, p1, p2, lower, upper, islog);
    tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames,
      model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT,
      matchcell, isr1, nsim, bw, ncore, gpuid, debug);
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

void MigrateDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll,  arma::cube& thetas, // nsub x npar x nchain
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, arma::uvec subgroups,
  unsigned int ngroup, double rp) {

  unsigned int start, end, next;
  double tmp_hlp, tmp_hll, tmp_logpos, cur_hlp, cur_hll, cur_logpos;
  arma::vec noise, tmp_loc, tmp_sca;
  arma::mat useloc = arma::trans(usephi[0]); // nchain x npar to npar x nchain
  arma::mat usesca = arma::trans(usephi[1]);
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  unsigned int l = subgroups.n_elem; // the number of groups + current group to migrate
  unsigned int m = nchain/ngroup;    // number of chains in each group
  arma::uvec subchains;

  for (size_t i = 0; i < l; i++) { // loop through subgroups, eg 0, 2, 3
    start = i*m; // where to start in overall chain index (eg); 0, 10, 15
    end   = ((1 + i)*m -1); // where to end in overall index (eg): 4, 14, 19
    subchains   =  arma::linspace<arma::uvec>(start, end, m);
    for (size_t j = 0; j < m; j++) { // loop through all chains in a subgroup
      // Infer the reference chain index in the next (cont' or noncont') subgroup;
      // If i is the last subgroup, its reference chain is in the first
      // subgroup; otherwise, its reference chain is in the next subgroup.
      next  = (i == l-1) ? subgroups(0) * m + j : subgroups(i+1) * m + j;
      noise = (2.0 * rp * arma::randu<arma::vec>(npar) - rp);
      tmp_loc = useloc.col(subchains(j)) + noise;
      tmp_sca = usesca.col(subchains(j)) + noise;
      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1,
        lp2, sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(thetas.slice(subchains(j)), pdists, tmp_loc,
        tmp_sca, plower, pupper, plog);
      tmp_logpos = tmp_hlp + tmp_hll;
      cur_logpos = usehlp(next) + usehll(next);

      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0,1) < std::exp(tmp_logpos - cur_logpos)) {
        useloc.col(next) = tmp_loc;
        usesca.col(next) = tmp_sca;
        usehlp(next)    = tmp_hlp;
        usehll(next)    = tmp_hll;
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
  arma::uvec matchcell, arma::uvec isr1, arma::uvec subgroups,
  unsigned int ngroup, unsigned int npda = 16384, double bw = .01,
  unsigned int ncore = 1, unsigned int gpuid = 0, bool debug = false,
  double rp = .001)
{
  arma::mat phi0 = arma::trans(usephi0);
  arma::mat phi1 = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;

  unsigned int l, m, start, end, next;
  double tmp_logpos, cur_logpos, tmp_lp, tmp_ll, cur_lp, cur_ll;
  arma::vec noise, tmp;
  l = subgroups.n_elem;
  m = nchain / ngroup;
  arma::uvec subchains;

  for (size_t i = 0; i < l; i++) {
    start = i*m;
    end   = ( (1 + i)*m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);
    for (size_t j = 0; j < m; j++) {
      next  = (i == l-1) ? subgroups(0) * m + j : subgroups(i+1) * m + j;
      noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
      tmp   = theta.col(subchains(j)) + noise;

      tmp_lp = sumlogprior(tmp, dists, phi0.col(subchains(j)),
        phi1.col(subchains(j)), lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames,
        model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
        npda, bw, ncore, gpuid, debug);
      tmp_logpos = tmp_lp + tmp_ll;

      // Current posterior; Same location in the next group
      // theta_cur = theta.col(next);
      // cur_lp = sumlogprior(theta_cur, dists, p1_.col(next),
      //   p2_.col(next), lower, upper, islog);
      // cur_ll = sumloglike(theta_cur, pnames, allpar, parnames, model,
      //   type, dim1, dim2, dim3, n1idx, ise, cellidx, RT,matchcell, isr1, npda,
      //   bw, ncore, gpuid, debug);
      cur_logpos = uselp(next) + usell(next);
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(next) = tmp;
        uselp(next) = tmp_lp;
        usell(next)  = tmp_ll;
      }
    }
  }
  usetheta = arma::trans(theta);
}

void MigrateDMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll,  arma::cube& theta,
  std::vector<std::string> pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, std::vector<std::string> ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  std::vector<std::string> sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp) {

  arma::vec noise; // nchain x npar to npar x nchain; /* h.migrate */
  arma::mat useloc = arma::trans(usephi(0));
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  arma::uvec subchains = GetSubchains(nchain);
  unsigned int nsubchain = subchains.n_elem;

  arma::mat tmp_loc(npar, nsubchain), tmp_sca(npar, nsubchain);
  tmp_loc.fill(NA_REAL); tmp_sca.fill(NA_REAL);
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain);
  arma::vec tmp_hlp(nsubchain), tmp_hll(nsubchain);

  for(size_t i = 0; i < nsubchain; i++) { // 0, 1, 2, 5, ...
    cur_hlp(i) = usehlp(subchains(i));
    cur_hll(i) = usehll(subchains(i));
    noise = (2.0 * rp * arma::randu<arma::vec>(npar) + rp);
    tmp_loc.col(i) = useloc.col(subchains(i)) + R::runif(-rp, rp);
    tmp_sca.col(i) = usesca.col(subchains(i)) + R::runif(-rp, rp);
    // thetas: nsub x npar x nchain
    tmp_hlp(i) = sumloghprior(tmp_loc.col(i), tmp_sca.col(i), ldists, sdists,
      lp1, sp1, lp2, sp2, llower, slower, lupper, supper, llog, slog);
    tmp_hll(i) = sumloghlike(theta.slice(subchains(i)), pdists, tmp_loc.col(i),
      tmp_sca.col(i), plower, pupper, plog);
  }

  /* Compare the first current chain to the last subchain ------------ */
  double tmp_logpos = tmp_hlp(nsubchain - 1) + tmp_hll(nsubchain - 1);
  double cur_logpos = cur_hlp(0) + cur_hll(0);
  if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  if (R::runif(0, 1) < std::exp( tmp_logpos - cur_logpos)) {
    useloc.col(subchains(0)) = tmp_loc.col(nsubchain - 1);
    usesca.col(subchains(0)) = tmp_sca.col(nsubchain - 1);
    usehlp(subchains(0)) = tmp_hlp(nsubchain - 1);
    usehll(subchains(0)) = tmp_hll(nsubchain - 1);
  }

  if (nsubchain != 1) {
     for (size_t j = 0; j < (nsubchain - 2); j++) {
        tmp_logpos = tmp_hlp(j) + tmp_hll(j);
        cur_logpos = cur_hlp(j+1) + cur_hll(j+1);
        if (std::isnan(tmp_logpos)) { tmp_logpos = -INFINITY; }
        if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
          useloc.col(subchains(j+1)) = tmp_loc.col(j);
          usesca.col(subchains(j+1)) = tmp_sca.col(j);
          usehlp(subchains(j+1)) = tmp_hlp(j);
          usehll(subchains(j+1)) = tmp_hll(j);
        }
      }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void MigrateDMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselogprior, arma::vec& useloglike,
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
  double gammamult = 2.38)
{
  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar    = theta.n_rows;
  unsigned int nchain  = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain); // eg. 0, 1, 3, 4, 8
  unsigned int nsubchain = subchains.n_elem;

  arma::mat tmp(npar, nsubchain);
  arma::vec cur_lp(nsubchain), cur_ll(nsubchain);
  arma::vec tmp_lp(nsubchain), tmp_ll(nsubchain), noise;

  for (size_t i = 0; i < nsubchain; i++) {
    noise = 2.0 * rp * arma::randu<arma::vec>(npar) + rp;
    tmp.col(i) = theta.col(subchains(i)) + noise;
    cur_lp(i)  = uselogprior(subchains(i));
    cur_ll(i)  = useloglike(subchains(i));

    tmp_lp(i) = sumlogprior(tmp.col(i), dists, phi0.col(subchains(i)),
      phi1.col(subchains(i)), lower, upper, islog);
    tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames,
      model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
      nsim, bw, ncore, gpuid, debug);
  }

  // Metropolis step 1; last chain vs first chain
  double tmp_logpos = tmp_lp(nsubchain - 1) + tmp_ll(nsubchain - 1);
  double cur_logpos = cur_lp(0) + cur_ll(0);
  if (std::isnan(tmp_logpos)) { tmp_logpos = -INFINITY; }

  if (R::runif(0,1) < std::exp(tmp_logpos - cur_logpos)) {
    theta.col(subchains(0))   = tmp.col(nsubchain - 1);
    uselogprior(subchains(0)) = tmp_lp(nsubchain - 1);
    useloglike(subchains(0))  = tmp_ll(nsubchain - 1);
  }

  if (nsubchain != 1) {
    for (size_t j = 0; j < (nsubchain - 2); j++) {
      tmp_logpos = tmp_lp(j) + tmp_ll(j);
      cur_logpos = cur_lp(j+1) + cur_ll(j+1);
      if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (R::runif(0,1) < std::exp(tmp_logpos - cur_logpos)) {
        theta.col(subchains(j + 1))   = tmp.col(j);
        uselogprior(subchains(j + 1)) = tmp_lp(j);
        useloglike(subchains(j + 1))  = tmp_ll(j);
      }
    }
  }
  usetheta = arma::trans(theta);
}

//' Fit a Bayesian Model to Data
//'
//' run_dgmc runs Bayesian modelling with DGMC sampler.
//'
//' @param samples a initialized sample
//' @param force PDA re-calculate interval
//' @param report progress report intervel
//' @param ncore number of CPU cores
//' @param pm probability of migration
//' @param qm probability of mutation
//' @param ngroup number of independent groups
//' @return Bayesian samples
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

    for (size_t j = 0; j < ngroup; j++) {

      // start = i*m;
      // end   = ((1 + i)*m - 1);
      // subchains = arma::linspace<arma::uvec>(start, end, m); // overall index
      subgroups = SelectEmigrants(ngroup, j);

      if (R::runif(0.0, 1.0) < pm) {
        MigrateDGMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
          plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
          n1idx, ise, cellidx, RT, mc, isr1, subgroups, ngroup, rp, force(i),
          npda, bw, ncore, gpuid, debug);

      } else if (R::runif(0.0, 1.0) <= qm) {
        MutateDGMCChains(usetheta, uselp, usell, pnames, pdists, pp1,
          pp2, plower, pupper, plog, allpar, parnames, model, type, dim1, dim2,
          dim3, n1idx, ise, cellidx, RT, mc, isr1, j, ngroup, rp, force(i), npda,
          bw, ncore, gpuid, debug);
      } else {
        CrossoverDGMCChains(usetheta, uselp, usell, pnames, pdists,
          pp1, pp2, plower, pupper, plog, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT, mc, isr1, j, ngroup, rp, force(i),
          npda, bw, ncore, gpuid, gammamult, debug);
      }
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
  double gammamult, unsigned int ncore) {

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
      MigrateDMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2, plower,
        pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3, n1idx,
        ise, cellidx, RT, mc, isr1, rp, gammamult, force(i), npda, bw, ncore,
        gpuid, debug);
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
List run_hyper_dmc(List samples, unsigned int report = 100, double pm = 0,
  double gammamult = 2.38, unsigned int ncore = 1, bool debug = false) {

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
  List ppprior  = hyper["pp.prior"]; /* Extract pprior & ppprior */
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
  arma::field<arma::cube> subtheta(nsub); // nchains x npar x nmc; nmc x nchain
  arma::field<arma::mat> usetheta(nsub), lp(nsub), ll(nsub);
  arma::field<arma::vec> uselp(nsub), usell(nsub), allpars(nsub), RTs(nsub);
  arma::field<arma::umat> n1idx(nsub), cellidx(nsub);
  arma::field<arma::uvec> matchcells(nsub), emptycells(nsub), isr1(nsub);
  arma::field<std::vector<std::string>> parnames(nsub), dim1s(nsub), dim2s(nsub), dim3s(nsub);
  arma::field<arma::ucube> models(nsub);

  SpreadSubs(samples_in, subtheta, usetheta, lp, uselp, ll, usell, substore,
    types, allpars, n1idx, matchcells, emptycells, cellidx, parnames, dim1s,
    dim2s, dim3s, isr1, models, npdas, bws, gpuids, RTs);

   for (size_t i = 1; i < nsamp; i++) {
       if (R::runif(0, 1) < pm) {
          MigrateDMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
            plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog, sdists,
            sp1, sp2, slower, supper, slog, rp);
       } else {
          CrossoverDMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
            plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
            sdists, sp1, sp2, slower, supper, slog, rp, gammamult);
       }

      /* usephi(0) and usephi(1): nchain x npar */
      for (size_t j = 0; j < nsub; j++) { // usethetas(j): nchain x npar

          // Because usephi may change in the CrossoverHyperchains or MigrateHyperchains
          uselp(j) = UpdatePriors(usetheta(j), pdists, usephi(0), usephi(1),
             plower, pupper, plog);

          if (R::runif(0, 1) < pm) {
            MigrateDMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j), isr1(j));
          } else {
            CrossoverDMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j), isr1(j));
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
List run_hyper_dgmc(List samples, arma::uvec force, unsigned int report,
  double pm, double qm, double gammamult, unsigned int ncore,
  unsigned int ngroup) {

  List samples_in(clone(samples));
  CheckHyperPnames(samples_in);

  List hyper = samples_in.attr("hyper");
  unsigned int npar = hyper["n.pars"];
  unsigned int nmc  = hyper["nmc"];
  unsigned int thin = hyper["thin"];
  unsigned int start= hyper["start"]; // start_R == 1;
  double rp = hyper["rp"];            // rp is defined in initialise
  unsigned int nsub = samples.size();
  unsigned int start_C = start - 1;   // start_C == 0;
  unsigned int store_i = start_C ;    // store_i == 0;
  unsigned int nsamp = 1 + (nmc - start) * thin;

  /* Hyperparameters, hyper logprior and hyper loglikelihood */
  std::vector<std::string> pnames = hyper("p.names");
  List phi      = hyper["phi"];
  arma::mat hlp = hyper["h_summed_log_prior"]; // nmc x nchain
  arma::mat hll = hyper["h_log_likelihoods"];
  arma::cube location = phi[0]; // nchain x npar x nmc
  arma::cube scale    = phi[1];
  arma::cube theta0  = GetTheta0(samples_in); // nsub x npar x nchain
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

  SpreadSubs(samples_in, subtheta, usetheta, lp, uselp, ll, usell, substore,
    types, allpars, n1idx, matchcells, emptycells, cellidx, parnames, dim1s,
    dim2s, dim3s, isr1, models, npdas, bws, gpuids, RTs);

  // DGMC specifics
  unsigned int l;
  arma::vec useLogPrior_, useLogLike_, current_theta;
  arma::uvec idx, subgroups;
  unsigned int nchain = theta0.n_slices;

  for (size_t i = 1; i < nsamp; i++) {

    for (size_t cur_gp_i = 0; cur_gp_i < ngroup; cur_gp_i++) {

      if (R::runif(0, 1) < pm) {
           l = (unsigned int)std::ceil( (double)(ngroup - 1) * R::runif(0.0, 1.0));
           // not necessary all groups are selected
           subgroups = SelectEmigrants(ngroup, cur_gp_i);
           MigrateDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
              plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
              sdists, sp1, sp2, slower, supper, slog, subgroups, ngroup, rp);
      } else if (R::runif(0, 1) <= qm) {
        MutateDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists, plower,
          pupper, plog, ldists, lp1, lp2, llower, lupper, llog, sdists, sp1,
          sp2, slower, supper, slog, ngroup, rp);

      } else {
        CrossoverDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists, plower,
          pupper, plog, ldists, lp1, lp2, llower, lupper, llog, sdists, sp1,
          sp2, slower, supper, slog, ngroup, rp, gammamult);
      }

    }

      /*  usephi: nchain x npar; Update data level and hyper data */
      for (size_t j = 0; j < nsub; j++) {
        uselp(j) = UpdatePriors(usetheta(j), pdists, usephi(0), usephi(1),
            plower, pupper, plog);

        if (R::runif(0, 1) < pm) {
            for (size_t cur_gp_j = 0; cur_gp_j < ngroup; cur_gp_j++) {
               subgroups = SelectEmigrants(ngroup, cur_gp_j);
               MigrateDGMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
               pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
               parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
               n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
               isr1(j), subgroups, ngroup);
            }

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

        if ( i % thin == 0 ) {
            substore(j)++;
            lp(j).row(substore(j)) = uselp(j).t();
            ll(j).row(substore(j)) = usell(j).t();
            subtheta(j).slice(substore(j)) = usetheta(j);
         }

         for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
            theta0.slice(k).row(j) = usetheta(j).row(k);
         }

      }

      // theta0 = GetUsethetas(usetheta);
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


// void mutation(arma::mat& usetheta, arma::vec& uselogprior, arma::vec& useloglike,  // nchain x 1
//                     std::vector<std::string> pnames, std::vector<std::string> dists,
//                     arma::vec p1, arma::vec p2, arma::vec lower,
//                     arma::vec upper, arma::uvec islog, arma::vec allpar,
//                     std::vector<std::string> parnames, arma::ucube model,
//                     std::string type,
//                     std::vector<std::string> dim1,
//                     std::vector<std::string> dim2,
//                     std::vector<std::string> dim3,
//                     arma::umat n1idx, arma::uvec ise, arma::umat cellidx,
//                     arma::vec RT, arma::uvec matchcell, arma::uvec isr1,
//                     unsigned int ngroup, double rp, bool force, int nsim,
//                     double bw, int ncore, int gpuid, bool debug) {
//
//   unsigned int npar, nchain, m, start, end;
//   double tmp_logprior, tmp_loglike, tmp_logpos, cur_logpos;
//   arma::vec noise, proposal;
//   arma::mat theta;
//   arma::uvec idx;
//   theta  = arma::trans(usetheta);    // theta: npar x nchain
//   npar   = theta.n_rows;
//   nchain = theta.n_cols;
//   m = nchain / ngroup;
//
//   for (size_t i = 0; i < ngroup; i++) {
//     start = i * m;
//     end   = ( (1 + i)*m - 1);
//     idx   = arma::linspace<arma::uvec>(start, end, m);
//     for (size_t j = 0; j < m; j++) {
//       // (b-a) * R::runif(1) + a;
//       proposal = theta.col(idx(j)) + (2.0 * rp * arma::randu<arma::vec>(npar) - rp);
//       tmp_logprior = sumlogprior(proposal, dists, p1, p2, lower, upper, islog);
//       tmp_loglike = sumloglike(proposal, pnames, allpar, parnames, model,
//         type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim,
//         bw, ncore, gpuid, debug);
//       tmp_logpos = tmp_logprior + tmp_loglike;
//       if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
//
//       cur_logpos = useloglike(idx(j)) + uselogprior(idx(j));
//       if (R::runif(0, 1) < std::exp(tmp_logpos - cur_logpos)) {
//         theta.col(idx(j)) = proposal;
//         uselogprior(idx(j)) = tmp_logprior;
//         useloglike(idx(j)) = tmp_loglike;
//       }
//     }
//   }
//   usetheta = arma::trans(theta);
// }
