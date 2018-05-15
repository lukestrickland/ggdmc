#include <RcppArmadillo.h>

bool checkDDM(std::vector<double> pVec) ;

Rcpp::NumericMatrix getAccumulatorMatrix(Rcpp::NumericVector pVec,
  std::string cell, Rcpp::NumericVector model, bool n1order);

arma::vec density_rd(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, double precision);

arma::vec density_norm(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1);

arma::vec density_norm_pda(arma::vec pVec, std::vector<std::string> pnames,
                           arma::vec allpar, std::vector<std::string> parnames,
                           arma::ucube model,
                           std::string type,
                           std::vector<std::string> dim1,
                           std::vector<std::string> dim2,
                           std::vector<std::string> dim3,
                           arma::umat n1idx, arma::uvec ise, arma::umat cellidx,
                           arma::vec RT, arma::uvec matchcell, arma::uvec isr1,
                           int nsim, double bw, bool debug);

arma::vec density_plba1(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, int nsim, double bw, int ncore, bool debug);


arma::vec density_plba1_gpu(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, int nsim,
  double bw, int ncore, int gpuid, int nthread, bool debug);

double sumloglike(arma::vec pvec,
                     std::vector<std::string> pnames,
                     arma::vec allpar,
                     std::vector<std::string> parnames,
                     arma::ucube model,
                     std::string type,
                     std::vector<std::string> dim1,
                     std::vector<std::string> dim2,
                     std::vector<std::string> dim3,
                     arma::umat n1idx, arma::uvec ise,
                     arma::umat cellidx, arma::vec RT,
                     arma::uvec matchcell, arma::uvec isr1, int nsim, double bw,
                     int ncore, int gpuid, bool debug);

