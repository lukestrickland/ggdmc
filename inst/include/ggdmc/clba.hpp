#include <RcppArmadillo.h>

arma::vec remove_t0(arma::vec rt, double t0);

arma::vec fptcdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
  double t0, bool posdrift);

arma::vec fptpdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
  double t0, bool posdrift);

arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b, arma::vec mean_v,
                       arma::vec sd_v, arma::vec t0);

arma::vec n1PDFfixedt0_pda(arma::vec rt, double A, double b,
                           arma::mat mean_v, arma::vec sd_v, double t0,
                           unsigned int n, double h, bool debug);

arma::mat rlba_norm(unsigned int n, arma::vec A, arma::vec b, arma::mat mean_v,
                    arma::vec sd_v, arma::vec t0, double st0, bool posdrift,
                    bool return_ttf, bool debug);

// arma::mat rlba(unsigned int n, double A, double b, double t0, arma::vec mean_v,
//   arma::vec sd_v);
