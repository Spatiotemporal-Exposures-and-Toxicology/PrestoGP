#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <cmath>

double matern1(const double &x, const double &smooth, const double &alpha)
{
  if (smooth == 0.5)
  {
    return exp(-x * alpha);
  }
  else if (smooth == 1.5)
  {
    double normcon = (pow(alpha, 3) / sqrt(M_PI)) * (tgamma(2.0) / tgamma(1.5));
    return normcon * 0.5 * M_PI * pow(alpha, -3) * exp(-x * alpha) * (1 + x * alpha);
  }
  else if (smooth == 2.5)
  {
    double normcon = (pow(alpha, 5) / sqrt(M_PI)) * (tgamma(3.0) / tgamma(2.5));
    double scaled_x = x * alpha;
    return normcon * 0.125 * M_PI * pow(alpha, -5) * exp(-scaled_x) * (3 + 3 * scaled_x + scaled_x * scaled_x);
  }
  else if (x == 0.0)
  {
    return 1.0;
  }
  else
  {
    double normcon = pow(2.0, 1.0 - smooth) / tgamma(smooth);
    double scaled_x = x * alpha;
    return normcon * pow(scaled_x, smooth) * R::bessel_k(scaled_x, smooth, 1);
  }
}

// computes Euclidean distance between 2 vectors ([slightly slower] helper function for rdist in C++)
double dist1(const arma::rowvec &a, const arma::rowvec &b)
{
  double temp = 0.0;
  // loop over each element and square the distance, add to the placeholder
  for (arma::uword i = 0; i < a.n_elem; i++)
  {
    temp += (a(i) - b(i)) * (a(i) - b(i));
  }
  return sqrt(temp);
}

// [[Rcpp::export]]
arma::vec na_omit_c(arma::vec x)
{
  // make placeholder vector r
  arma::vec r(x.size());
  // make counter for number of non-NAs
  int k = 0;
  for (unsigned j = 0; j < x.size(); j++)
  {
    // if not NA, add 1 to counter and set r(j) = x(j)
    if (x(j) == x(j))
    {
      r(j) = x(j);
      k++;
    }
  }
  // resize r to non-NA size (removing NAs)
  r.resize(k);
  // subtract one because C++ indexes from 0 while inputs from R index from 1
  // r -= 1.0;
  return r;
}

// [[Rcpp::export]]
arma::sp_mat createU_helper_mat(const arma::mat &olocs, const arma::vec &ondx,
                                const arma::mat &curqys, const arma::mat &curqzs,
                                const arma::mat &vijs, const arma::mat &aijs,
                                const arma::mat &full_const, const arma::vec &nugget,
                                const arma::vec &sig2, const arma::vec &U_beginning)
{
  int n = arma::as_scalar(ondx.n_elem);
  int m = arma::as_scalar(curqys.n_rows);
  int n_inds = 2 * n * (m + 3);
  n_inds -= sum(arma::linspace(0, m, m + 1));
  // arma::mat feeder(3, n_inds);
  arma::umat ndx(2, n_inds);
  arma::vec vals(n_inds);
  // vals.rows(0, 6) = U_beginning;
  // feeder.cols(0, 6) = {{1, 1, 2, 1, 3, 3, 4},
  //                     {1, 2, 2, 3, 3, 4, 4}};
  ndx.cols(0, 6) = {{0, 0, 1, 0, 2, 2, 3},
                    {0, 1, 1, 2, 2, 3, 3}};
  vals(arma::span(0, 6)) = U_beginning;
  // feeder.cols(0,6) = {{0, 0, 1, 0, 2, 2, 3},
  //                     {0, 1, 1, 2, 2, 3, 3}, {1,1,1,1,1,1,1}};
  // feeder(2, arma::span(0,6)) = U_beginning.t();
  int ind = 7;
  // int ind = 0;
  for (arma::uword i = 3; i < (ondx.n_elem + 1); i++)
  {
    double temppi = ondx(i - 1);
    arma::vec cy = na_omit_c(curqys.col(i - 1));
    arma::vec cz = (curqzs.col(i - 1));
    arma::vec cq = na_omit_c(arma::join_cols(cy, cz));
    arma::uword nq = cq.n_elem;
    arma::vec k1(nq);
    arma::mat k2(nq, nq);
    // Rcpp::Rcout << i << " ";
    // Rcpp::Rcout << cq;
    for (arma::uword j = 0; j < nq; j++)
    {
      double temppj = ondx(cq(j) - 1);
      k1(j) = full_const(temppi - 1, temppj - 1) *
              matern1(dist1(olocs.row(i - 1), olocs.row(cq(j) - 1)),
                      vijs(temppi - 1, temppj - 1), aijs(temppi - 1, temppj - 1));
      for (arma::uword k = j; k < nq; k++)
      {
        double temppk = ondx(cq(k) - 1);
        k2(j, k) = full_const(temppj - 1, temppk - 1) *
                   matern1(dist1(olocs.row(cq(j) - 1), olocs.row(cq(k) - 1)),
                           vijs(temppj - 1, temppk - 1), aijs(temppj - 1, temppk - 1));
      }
      // // add nugget if neighbor is not latent
      if (j >= cy.n_elem)
      {
        k2(j, j) += nugget(temppj - 1);
      }
    }
    // //make symmetric
    k2 += k2.t();
    // //but diagonal is double counted
    k2.diag() /= 2.0;
    // Rcpp::Rcout << k2;
    arma::vec bi(nq);
    bool status = arma::solve(bi, k2, k1);
    double epsilon = 1e-9;
    while (!status)
    {
      epsilon *= 2;
      k2.diag() += epsilon;
      status = arma::solve(bi, k2, k1);
    }
    double ritemp = arma::as_scalar(sig2(temppi - 1) - bi.t() * k1);
    epsilon = 1e-9;
    while (ritemp <= 0)
    {
      epsilon *= 2;
      k2.diag() += epsilon;
      arma::solve(bi, k2, k1);
      ritemp = arma::as_scalar(sig2(temppi - 1) - bi.t() * k1);
      // Rcpp::Rcout << ritemp << "\n";
    }
    double ri = sqrt(ritemp);
    bi *= (-1.0 / ri);
    // arma::vec bi = arma::solve(k2, k1);
    // uhat(2*i - 1, 2*i - 1) = pow(nugget(temppi - 1), -0.5);
    // uhat(2*i - 2, 2*i - 1) = -1.0 * uhat(2*i - 1, 2*i - 1);

    // uhat(2*i - 2, 2*i - 2) = 1.0 / ri;
    for (arma::uword j = 0; j < cq.n_elem; j++)
    {
      arma::uword xind = 2 * cq(j) - 2;
      if (j >= cy.n_elem)
      {
        xind += 1;
      }
      // feeder.col(ind) = {(double)xind,(double) 2*i - 2, bi(j)};
      ndx.col(ind) = {xind, 2 * i - 2};
      vals(ind) = bi(j);
      // vals[ind] = bi(j);
      ind++;
      // uhat(xind, 2.0*i - 2.0) = bi(j);
    }
    // vals[ind] = 1.0 / ri;
    double feed_temp = pow(nugget(temppi - 1), -0.5); // used twice so temp storage to avoid recomputation
    // feeder.col(ind) = {(double)2*i - 2,(double)2*i - 2, 1.0 / ri};
    // feeder.col(ind + 2) = {(double)2*i - 1,(double) 2*i - 1, feed_temp};
    ndx.col(ind) = {2 * i - 2, 2 * i - 2};
    vals(ind) = 1.0 / ri;
    ndx.col(ind + 2) = {2 * i - 1, 2 * i - 1};
    vals(ind + 2) = feed_temp;
    // vals[ind + 2] = pow(nugget(temppi - 1), -0.5);
    // vals[ind + 1] = -1.0 * vals[ind + 2];
    // feeder.col(ind + 1) = {(double)2*i - 2,(double) 2*i - 1, -1.0 * feed_temp};
    ndx.col(ind + 1) = {2 * i - 2, 2 * i - 1};
    vals(ind + 1) = -1.0 * feed_temp;
    ind += 3;
  }
  arma::sp_mat uhat(ndx, vals);
  // 0 indexing -> 1 indexing for R
  // feeder.row(0) += 1.0;
  // feeder.row(1) += 1.0;
  // return feeder;
  return uhat;
}
