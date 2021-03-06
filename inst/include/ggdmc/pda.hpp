#include <RcppArmadillo.h>

arma::vec spdf(arma::vec x, arma::vec RT, int n, double h_in, bool debug);

arma::mat rplba1(int n, arma::vec A, arma::vec b, double t0, arma::vec mean_v,
  arma::vec mean_w, arma::vec sd_v, double rD, double swt, int ncore,
  bool debug);

arma::mat rplba1_test(int n, arma::vec A, arma::vec b, double t0, arma::vec mean_v,
  arma::vec mean_w, arma::vec sd_v, double rD, double swt, bool debug);

arma::mat rplba2(int n, arma::vec A, arma::vec b, double t0, arma::vec mean_v,
  arma::vec mean_w, arma::vec sd_v, arma::vec sd_w, double rD,
  double swt, int ncore, bool debug);

arma::mat rplba2_test(int n, arma::vec A, arma::vec b, double t0, arma::vec mean_v,
  arma::vec mean_w, arma::vec sd_v, arma::vec sd_w, double rD,
  double swt, int ncore, bool debug);

arma::vec n1PDF_plba1(arma::vec x, arma::vec A, arma::vec b, arma::vec mean_v,
  arma::vec sd_v, double t0, arma::vec mean_w, double rD, double swt, int n,
  double h, int ncore, bool debug);

arma::vec n1PDF_plba2(arma::vec x, arma::vec A, arma::vec b, arma::vec mean_v,
  arma::vec sd_v, double t0, arma::vec mean_w, arma::vec sd_w,
  double rD, double swt, int n, double h, int ncore, bool debug);

arma::vec n1PDF_plba3(arma::vec x, int n, arma::vec A, arma::vec B,
  arma::vec C, arma::vec mean_v, arma::vec sd_v, arma::vec mean_w,
  arma::vec sd_w, double rD, double tD, double swt, double t0, double h);

