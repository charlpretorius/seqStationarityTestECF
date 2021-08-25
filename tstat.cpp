#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat tstats(arma::vec x, int T, int m, arma::vec a) {
  int nY = x.n_elem - m + 1;
  int na = a.n_elem;
  arma::mat Y(nY, m);
  arma::mat Wjk(nY, nY);
  
  for (int i = 0; i < nY; i++) {
    Y.row(i) = x.subvec(i, i + m - 1).t();
  }
  
  for (int j = 0; j < nY; j++) {
    for (int k = 0; k < nY; k++) {
      Wjk(j,k) = exp(-sum(pow(Y.row(j) - Y.row(k), 2))/4);
    }
  }
  
  int Tm = T - m + 1;
  int Tmt = 0;
  arma::vec tstat(nY - T);
  arma::mat tstat1(na, nY - T);
  arma::mat Wjka(nY, nY);
  for (int k = 0; k < na; k++) {
    Wjka = pow(Wjk, 1/a(k));
    for (int tt = 0; tt < nY - T; tt++) {
      Tmt = Tm + tt + 1;
      
      tstat(tt) = sum(sum(Wjka.submat(0,0,Tm-1,Tm-1)))/pow(Tm,2) +
        sum(sum(Wjka.submat(0,0,Tmt-1,Tmt-1)))/pow(Tmt,2) -
        2*sum(sum(Wjka.submat(0,0,Tmt-1,Tm-1)))/Tm/Tmt;
      
      tstat(tt) = pow((double) Tmt/(T + tt + 1),2)*tstat(tt);
    }
    tstat1.row(k) = pow(datum::pi/a(k), m/2)*pow(T,2)/Tm*tstat.t();
  }
  
  return tstat1;
}


