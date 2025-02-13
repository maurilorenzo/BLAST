#include <RcppArmadillo.h>
#include <Rcpp.h>



// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List sample_Lambda(const arma::mat& Y, 
                   const arma::mat& mu_js, 
                   double tau_2, 
                   const arma::mat& M, 
                   double gamma_0 = 1, 
                   double delta_0 = 1, 
                   int N_mc = 1000, 
                   double rho = 1, 
                   arma::uvec subsample_index = arma::regspace<arma::uvec>(0, 99),
                   bool sample_outer_product=true) {
  
  int n = M.n_rows;
  int p = mu_js.n_rows;
  int k = M.n_cols;
  int p_sample = subsample_index.n_elem;
  
  arma::cube Lambda_samples(p, k, N_mc, arma::fill::zeros);
  arma::cube Lambda_outer_samples(p_sample, p_sample, N_mc, arma::fill::zeros);
  arma::mat Lambda_outer_mean(p, p, arma::fill::zeros);
  
  arma::mat Lambda_t(p, k, arma::fill::zeros);
  arma::mat Sigmas_samples(N_mc, p, arma::fill::zeros);
  arma::vec residuals(p, arma::fill::zeros);
  arma::mat residual_matrix = arma::eye(n, n) - M * M.t() / n;

  for(int j = 0; j < p; ++j) {
    residuals[j] = as_scalar(Y.col(j).t() * residual_matrix * Y.col(j));
  }
  
  if (residuals.has_nan()) {
    Rcpp::Rcout << "NA in residuals" << std::endl;
    return List::create();
  }
  
  if (arma::any(residuals < 0)) {
    Rcpp::Rcout << "negative vals in residuals" << std::endl;
    return List::create();
  }
  
  arma::mat D = arma::diagmat((residuals + gamma_0 * delta_0) / (gamma_0 + n - 2)) 
    * k * rho * rho / (n + 1 / tau_2);
  Lambda_outer_mean = mu_js * mu_js.t() + D;
  

  
  // Sampling process
  for (int it = 0; it < N_mc; ++it) {
    for (int j = 0; j < p; ++j) {
      Sigmas_samples(it, j) = 1.0 / R::rgamma((gamma_0 + n) / 2, 2.0 / (gamma_0 * delta_0 + residuals[j]));
    }
  
    
    arma::mat eps_Lambda = rho * arma::randn<arma::mat>(p, k);
    Lambda_t = mu_js + eps_Lambda.each_col() % arma::sqrt(Sigmas_samples.row(it).t() / (n + 1 / tau_2));
    
    if(sample_outer_product){
      arma::mat Lambda_outer = Lambda_t.rows(subsample_index) * Lambda_t.rows(subsample_index).t();
      Lambda_outer_samples.slice(it) = Lambda_outer;
    }

    Lambda_samples.slice(it) = Lambda_t;
  }
  
  return List::create(
    Named("Lambda_outer_samples") = Lambda_outer_samples,
    Named("Lambda_samples") = Lambda_samples,
    Named("Sigmas_samples") = Sigmas_samples,
    Named("Lambda_outer_mean") = Lambda_outer_mean
  );
}


// [[Rcpp::export]]
arma::cube sample_Lambda_outer(arma::cube Lambda_samples){
  
  int n_MC = Lambda_samples.n_slices;
  int k = Lambda_samples.n_cols;
  int p = Lambda_samples.n_rows;
  
  arma::mat Lambda(p, k);
  arma::mat Lambda_outer(p, p);
  arma::cube Lambda_outer_samples(p, p, n_MC);

  for(int s=0; s<n_MC; ++s) {
    Lambda = Lambda_samples.slice(s);
    Lambda_outer = Lambda * Lambda.t();
    Lambda_outer_samples.slice(s) = Lambda_outer;
  }
  
  return Lambda_outer_samples;
}




// [[Rcpp::export]]
List sample_Gamma(const arma::mat& Y, 
                  double tau_2, 
                  const arma::mat& F_hat, 
                  const arma::cube& Lambda_samples, 
                  const arma::mat& M, 
                  const arma::mat& Sigmas_samples, 
                  double gamma_0 = 1, 
                  double delta_0 = 1, 
                  int N_mc = 100, 
                  double rho = 1, 
                  arma::uvec subsample_index = arma::regspace<arma::uvec>(0, 99)) {
  
  int n = M.n_rows;
  int p = Y.n_cols;
  int q = F_hat.n_cols;
  int k = M.n_cols;
  int p_sample = subsample_index.n_elem;
  
  arma::cube Gamma_samples(p , q, N_mc);
  arma::cube Gamma_outer_samples(p_sample, p_sample, N_mc, arma::fill::zeros);
  arma::mat Gamma_outer_mean(p, p, arma::fill::zeros);
  arma::mat Gamma_t(p, q, arma::fill::zeros);
  arma::mat Gamma_t_outer(p_sample, p_sample);
  
  for(int it = 0; it < N_mc; ++it) {
    arma::vec Sigmas_sample = Sigmas_samples.row(it).t();
    arma::mat Lambda_t = Lambda_samples.slice(it);
    arma::mat Y_res = Y - M * Lambda_t.t();
    arma::mat mean_Gamma = Y_res.t() * F_hat / (n + 1 / tau_2);
    arma::mat eps_Gamma = rho*arma::randn<arma::mat>(p, q);
    Gamma_t = mean_Gamma + eps_Gamma.each_col() % arma::sqrt(Sigmas_sample / (n + 1 / tau_2));
    if(p_sample < p){
     Gamma_t_outer = Gamma_t.rows(subsample_index) * Gamma_t.rows(subsample_index).t();
    } else{
     Gamma_t_outer = Gamma_t * Gamma_t.t();
    }
    Gamma_outer_samples.slice(it) = Gamma_t_outer;
    Gamma_samples.slice(it) = Gamma_t;
  }
  
  return List::create(
    Named("Gamma_outer_samples") = Gamma_outer_samples,
    Named("Gamma_samples") = Gamma_samples,
    Named("Gamma_outer_mean") = Gamma_outer_mean
  );
}


// [[Rcpp::export]]
arma::mat compute_B_Lambda(const arma::mat& Lambda, const arma::vec& Sigma_2s) {
  arma::mat Lambda_outer = Lambda * Lambda.t(); 
  int p = Lambda.n_rows; 
  arma::mat B = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      B(j, l) = sqrt(1 + (Lambda_outer(j, j) * Lambda_outer(l, l) + Lambda_outer(j, l) * Lambda_outer(j, l)) / 
        (Sigma_2s(j) * Lambda_outer(l, l) + Sigma_2s(l) * Lambda_outer(j, j)));
      B(l, j) = B(j, l); 
    }
    B(j, j) = sqrt(1 + Lambda_outer(j, j) / (2 * Sigma_2s(j)));
  }
  B(p - 1, p - 1) = sqrt(1 + Lambda_outer(p - 1, p - 1) / (2 * Sigma_2s(p - 1)));
  
  return B;
}

// [[Rcpp::export]]
arma::mat compute_B_Gamma(const arma::mat& Lambda, const arma::mat& Gamma, const arma::vec& Sigma_2s) {
  arma::mat Gamma_outer = Gamma * Gamma.t(); 
  arma::mat Lambda_outer = Lambda * Lambda.t();
  int p = Gamma.n_rows;
  arma::mat B = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      B(j, l) = sqrt(1 + (
        Gamma_outer(j, j) * Gamma_outer(l, l) + Gamma_outer(j, l) * Gamma_outer(j, l) +
          Lambda_outer(j, j)*Gamma_outer(l, l) + Lambda_outer(l, l)*Gamma_outer(j, j) + 
          2 * Lambda_outer(l, j)*Gamma_outer(j, l)) / 
          (Sigma_2s(j) * Gamma_outer(l, l) + Sigma_2s(l) * Gamma_outer(j, j)));
      B(l, j) = B(j, l); 
    }
    B(j, j) = sqrt(1 + (Gamma_outer(j, j) + 2*Lambda_outer(j, j)) / (2 * Sigma_2s(j)));
  }
  B(p - 1, p - 1) = sqrt(1 + (Gamma_outer(p - 1, p - 1) + 2*Lambda_outer(p - 1, p - 1)) 
                           / (2 * Sigma_2s(p - 1))); 
  
  return B;
}


// [[Rcpp::export]]
List compute_svd_cpp(arma::mat Y, int flag=1){
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  int r = std::min(n,p);
  
  arma::mat U(n, n, fill::zeros);
  arma::mat V(p, r, fill::zeros);
  arma::vec d(r, fill::zeros);
  
  if(flag == 1){
    svd(U, d, V, Y, "std");
  } else if(flag == 2){
    svd(U, d, V, Y, "dc");
  } else if(flag == 3){
    svd_econ(U, d, V, Y, "both", "dc");
  } else if(flag == 4){
    svd_econ(U, d, V, Y, "left", "dc");
  }
  
  return List::create(
    Named("u") = U, Named("d") = d, Named("v") = V
  );
  
}