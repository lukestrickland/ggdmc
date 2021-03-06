#include <ggdmc.hpp>
#include <boost/foreach.hpp>

using namespace Rcpp;

#define foreach BOOST_FOREACH

// [[Rcpp::depends(BH)]]

//' @export
// [[Rcpp::export]]
arma::vec pnormP(arma::vec x, double mean=0, double sd=1, double lt=true,
  double lg=false) {
  arma::vec out(x.n_elem);   // protected pnorm & dnorm
  for (size_t i = 0; i < x.n_elem; i++) {
    if (std::abs(x(i)) < 7) {
      out(i) = R::pnorm(x(i), mean, sd, lt, lg);
    } else if (x(i) < 0) {
      out(i) = 0;
    } else {
      out(i) = 1;
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec dnormP(arma::vec x, double mean=0, double sd=1, double lg=false) {
  arma::vec out(x.n_elem);
  for(size_t i = 0; i < x.n_elem; i++) {
    if (std::abs(x(i)) < 7) {
      out(i) = R::dnorm(x(i), mean, sd, lg);
    } else {
      out(i) = 0;
    }
  }
  return out;
}

//' @rdname removet0
//' @export
// [[Rcpp::export]]
arma::vec remove_t0(arma::vec rt, double t0) {
  foreach( double& elem, rt ) {
    elem -= t0;
    if (elem < 0) elem = 0;
  }
  return rt;
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::mat ttf(arma::mat drifts, arma::vec A, arma::vec b, arma::vec t0,
              double st0, unsigned int nmean_v, unsigned int n) {
  // drifts is n x nmean_v matrix
  arma::mat out(n, nmean_v); // time to finish
  // Generate uniform random numbers: a + (b-a)*randu(n)

  unsigned int nA  = A.n_elem;
  unsigned int nb  = b.n_elem;
  unsigned int nt0 = t0.n_elem;
  if (nA == 1) A = arma::repmat(A, nmean_v, 1);
  if (nb == 1) b = arma::repmat(b, nmean_v, 1);
  if (nt0 == 1) t0 = arma::repmat(t0, nmean_v, 1);

  for (size_t i = 0; i < nmean_v; i++) {
    out.col(i) = t0(i) + (b(i) - A(i) * arma::randu(n)) / drifts.col(i);
  }
  return out.t();
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::vec fptcdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
  double t0 = 0, bool posdrift = true) {

  double x, tv, ts, term1, term2, tmp, denom;
  arma::vec out(rt.n_elem), dt;

  dt = remove_t0(rt, t0);
  denom = !posdrift ? 1.0 : std::max(R::pnorm(mean_v/sd_v, 0, 1, true, false), 1e-10);
  // denom = std::max(R::pnorm(mean_v/sd_v, 0, 1, true, false), 1e-10);

  for (size_t i = 0; i < rt.n_elem; i++) {
    x  = b / dt(i);         // SB's terminology
    tv = dt(i) * mean_v;    // zu; chiminuszu=b-tv; HS's xx=b-A-tv
    ts = dt(i) * sd_v;      // zs; chizu=(b-tv)/ts; chizumax=(b-A-tv)/ts
    term1 = (b-A-tv)*R::pnorm((b - A - tv)/ts, 0, 1, true, false) -
      (b-tv)*R::pnorm((b - tv)/ts, 0, 1, true, false); // HS's tmp2
    term2 = ts * (R::dnorm((b - A - tv)/ts, 0, 1, false) -
      R::dnorm((b-tv)/ts, 0, 1, false)); // HS's tmp1
    if (A < 1e-10) {
      tmp = R::pnorm(x, mean_v, sd_v, false, false) / denom;

      out(i) = std::min(1.0, std::max(0.0, tmp));

    } else {

      tmp = (1.0 + (term1 + term2)/A) / denom;
      out(i) = std::min(1.0, std::max(0.0, tmp));
    }
  }
  return out;
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::vec fptpdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
  double t0 = 0, bool posdrift = true) {

  arma::vec dt = remove_t0(rt, t0);
  double x, tv, ts, term1, term2, denom, tmp;
  arma::vec out(dt.n_elem);
  denom = !posdrift ? 1.0 : std::max(R::pnorm(mean_v/sd_v, 0, 1, true, false), 1e-10);

  // if (b < A) {Rcpp::stop("b must be greater than A. in fptpdf");}
  for (size_t i = 0; i < dt.n_elem; i ++) {
    x  = b / dt(i);       // SB's terminology
    tv = dt(i) * mean_v;  // zu; chiminuszu=b-tv
    ts = dt(i) * sd_v;    // zs; chizu=(b-tv)/ts; chizumax=(b-A-tv)/ts
    term1 = mean_v *(R::pnorm((b-tv)/ts, 0, 1, 1, 0) - R::pnorm((b-A-tv)/ts,0,1,1, 0));
    term2 = sd_v*(R::dnorm((b-A-tv)/ts, 0, 1, 0)    - R::dnorm((b-tv)/ts,  0,1,0));

    if (A < 1e-10) {
      tmp = ( b / (dt(i)*dt(i)) ) * R::dnorm(x, mean_v, sd_v, 0) / denom;
      out(i) = std::max(0.0, tmp);
    } else {
      tmp = (term1 + term2)/ (A*denom); // HS's out_o
      out(i) = std::max(0.0, tmp);
    }
  }
  return out;
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b, arma::vec mean_v,
                       arma::vec sd_v, arma::vec t0) {

  unsigned int nmean_v = mean_v.n_elem;  // Number of accumulators/responses.
  unsigned int n       = rt.n_elem;       // Number of trials
  unsigned int nsd_v   = sd_v.n_elem; // Check for matrix operations
  unsigned int nA      = A.n_elem;
  unsigned int nb      = b.n_elem;
  unsigned int nt0     = t0.n_elem;

  if (nsd_v == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
  if (nA == 1) A = arma::repmat(A, nmean_v, 1);
  if (nb == 1) b = arma::repmat(b, nmean_v, 1);
  if (nt0 == 1) t0 = arma::repmat(t0, nmean_v, 1);

  arma::vec onevec = arma::ones<arma::vec>(n);
  arma::vec node1den = fptpdf(rt, A(0), b(0), mean_v(0), sd_v(0), t0(0));

  if (nmean_v > 1) {
    for (size_t i = 1; i < nmean_v; i++) {
      node1den = node1den % (onevec - fptcdf(rt, A(i), b(i), mean_v(i), sd_v(i), t0(i)));
      // if (posdrift) {
      //   node1den = node1den % fptcdf(rt, A(i), b(i), mean_v(i), sd_v(i), t0(i));
      // } else {
      //   node1den = node1den % (onevec - fptcdf(rt, A(i), b(i), mean_v(i), sd_v(i), t0(i)));
      // }
    }
  }
  return node1den;
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::vec n1PDFfixedt0_pda(arma::vec rt, double A, double b,
                           arma::mat mean_v, arma::vec sd_v, double t0,
                           unsigned int n, double h, bool debug)
{
  // Convert double to arma::vec
  arma::vec A_in(1), b_in(1), t0_in(1);
  A_in = A;
  b_in = b;
  t0_in = t0;
  // mean_v must be a n_acc x 1 matrix or n_acc x n matrix
  arma::mat sim = rlba_norm(n, A_in, b_in, mean_v, sd_v, t0_in, 0, true, false, debug);
  arma::vec sRT = sim.col(0);
  arma::vec sR  = sim.col(1);
  arma::vec RT0 = sRT.elem(arma::find(sR==1));  // return only R == 1
  arma::vec out = spdf(rt, RT0, n, h, debug);
  return out;
}


//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::mat make_r(arma::mat drifts, arma::vec A, arma::vec b, arma::vec t0,
                 double st0,
                 unsigned int nmean_v, unsigned int n, bool return_ttf,
                 bool debug = false) {
  // internal_rt = nmean_v x n matrix;
  if (drifts.n_cols != nmean_v) stop("ncol in drifts != nmean_v");

  unsigned int nA  = A.n_elem;
  unsigned int nb  = b.n_elem;
  unsigned int nt0 = t0.n_elem;
  if (nA == 1) A = arma::repmat(A, nmean_v, 1);
  if (nb == 1) b = arma::repmat(b, nmean_v, 1);
  if (nt0 == 1) t0 = arma::repmat(t0, nmean_v, 1);

  arma::mat internal_rt = ttf(drifts, A, b, t0, st0, nmean_v, n);

  arma::vec RT(n), R(n);
  for (size_t i = 0; i < n; i++) {
       R(i)  = 1 + internal_rt.col(i).index_min();
       RT(i) = internal_rt.col(i).min(); // which accumulator is the minimal RT
       if (st0 > 0) RT(i) = RT(i) + R::runif(0.0, st0);
  }

  arma::mat out;
  if (RT.has_inf()) {
     arma::vec finite_rt, finite_r;
     arma::uvec finite_idx = arma::find_finite(RT);
     finite_rt = RT.elem(finite_idx);
     finite_r  = 1 + R.elem(finite_idx);
     if (debug) Rcout << "Inf/NaN found. Only " << finite_idx.n_elem << " RTs are returned\n";
     out = arma::join_horiz(finite_rt, finite_r);
  } else {
     out = arma::join_horiz(RT, R);
  }
  if (return_ttf) { return internal_rt; } else { return out; }
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::mat make_v(unsigned int n, arma::mat mean_v, arma::vec sd_v,
                 bool posdrift = true)
{
  // posdrift = positive drift rates
  // drifts = n x n_v matrix; n_v is number of accumulators/responses
  // DMC's drifts is n_v x n
  unsigned int nmean_v = mean_v.n_rows;
  arma::mat drifts(n, nmean_v);
  if (sd_v.n_elem == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
  if (sd_v.n_elem != nmean_v) stop("sd_v must either be 1 or nmean_v elements.");

  for (size_t i = 0; i < nmean_v; i++) {
    if (posdrift) {
      drifts.col(i) = rtnorm(n,  mean_v(i), sd_v(i), 0, INFINITY);
    } else {
      drifts.col(i) = sd_v(i) * arma::randn(n) + mean_v(i);
    }
  }

  return drifts;
}

//' @rdname makeR
//' @export
// [[Rcpp::export]]
arma::mat rlba_norm(unsigned int n, arma::vec A, arma::vec b, arma::mat mean_v,
                    arma::vec sd_v, arma::vec t0, double st0 = 0,
                    bool posdrift = true, bool return_ttf = false,
                    bool debug = false)
{
  //if (b  < A && debug) Rcout << "b smaller than A!\t";
  //if (t0 < 0 && debug) Rcout << "t0 is negative!\t";
  unsigned int nmean_v = mean_v.n_rows;
  unsigned int nA  = A.n_elem;
  unsigned int nb  = b.n_elem;
  unsigned int nt0 = t0.n_elem;
  // unsigned int nsd_v   = sd_v.n_elem; // Check for matrix operations

  // drifts is n x nmean_v matrix
  if (nA == 1) A = arma::repmat(A, nmean_v, 1);
  if (nb == 1) b = arma::repmat(b, nmean_v, 1);
  if (nt0 == 1) t0 = arma::repmat(t0, nmean_v, 1);
  if (sd_v.n_elem == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);

  arma::mat drifts, out;
  drifts = make_v(n, mean_v, sd_v, posdrift);
  out    = make_r(drifts, A, b, t0, st0, nmean_v, n, return_ttf, debug);
  return out;
}

// arma::vec n1PDFfixedt0(arma::vec x, arma::vec A, arma::vec b, arma::vec mean_v,
//                        arma::vec sd_v, arma::vec t0) {
//
//   unsigned int nmean_v = mean_v.n_elem;  // Number of accumulators/responses.
//   unsigned int n       = x.n_elem;       // Number of trials
//   unsigned int nsd_v   = sd_v.n_elem; // Check for matrix operations
//   unsigned int nA      = A.n_elem;
//   unsigned int nb      = b.n_elem;
//   unsigned int nt0     = t0.n_elem;
//
//   if (nmean_v < 2) Rcpp::stop("Minimal 2 accumulators/responses!");
//   if (nsd_v == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
//   if (nA == 1)       A = arma::repmat(A, nmean_v, 1);
//   if (nb == 1)       b = arma::repmat(b, nmean_v, 1);
//   if (nt0 == 1)     t0 = arma::repmat(t0, nmean_v, 1);
//
//   arma::vec G(n);
//   arma::mat tmp(n, nmean_v - 1);
//   arma::vec onevec = arma::ones<arma::vec>(n);
//
//   if (nmean_v > 2) {
//     // three or more accumulators; calculate starting from 2nd, 3rd etc.
//     for(size_t i = 1; i < nmean_v; i++)
//       tmp.col(i-1) = onevec - fptcdf(x, A(i), b(i), mean_v(i), sd_v(i), t0(i));
//     for(size_t j = 0; j < n; j++) G(j) = arma::prod(tmp.row(j));
//
//   } else {
//     // calculate 2nd accumulator in a 2-accumulator model
//     G = onevec - fptcdf(x, A(1), b(1), mean_v(1), sd_v(1), t0(1));
//   }
//
//   return G % fptpdf(x, A(0), b(0), mean_v(0), sd_v(0), t0(0));
// }
