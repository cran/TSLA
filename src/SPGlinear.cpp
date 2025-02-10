#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// the projection operator of the generalized lasso
// [[Rcpp::export]]
arma::vec hardThresholdCpp(arma::mat u, double thres) {
  arma::mat thresvec = ones(u.n_elem) * thres;
  arma::vec out = min(u, thresvec);
  out = max(- thresvec, out);
  return(out);
}

// the projection operator of the group lasso
// [[Rcpp::export]]
arma::vec shrinkGroupCpp(arma::vec u, arma::mat g_idx) {
  arma::vec R = zeros(u.n_elem);
  int V = g_idx.n_rows;
  for(int v = 0; v < V; v ++) {
    arma::vec idx = g_idx.row(v).cols(0, 1).t();
    arma::vec unorm = max(ones(1), sqrt(sum(pow(u.rows(idx(0) - 1, idx(1) - 1), 2))));
    R.rows(idx(0) - 1, idx(1) - 1) = u.rows(idx(0) - 1, idx(1) - 1) / as_scalar(unorm);
  }
  return(R);
}

// calculate group norm for each gamma coefficient

//' Calculate group norms
//'
//' Function to output group norms on the gamma coefficients
//' based on \code{g_idx} and \code{C2} matrix.
//'
//' @param u \code{C2*gamma.coef}, \code{gamma.coef} is the estimated node coefficient vector,
//' \code{C2} matrix is the output from
//' function \code{get_tree_object()}, which gives the weights of the groups.
//' @param g_idx Group structure matrix defined by the \code{C2} matrix.
//' See details in \code{get_tree_object()}.
//' @param type If \code{type == 1}, return sum of group norms; else return individual
//' norm for each group.
//' @return Sum of group norms or individual group norms.
//' @export
// [[Rcpp::export]]
arma::vec cal2norm(arma::vec u, arma::mat g_idx, int type) {
  int V = g_idx.n_rows;
  arma::vec unorm = zeros(V);
  arma::vec sumnorm = zeros(1);
  arma::vec idx = zeros(2);
  for(int v = 0; v<V; v++) {
    idx = g_idx.row(v).cols(0, 1).t();
    unorm.row(v) = sqrt(sum(pow(u.rows(idx(0) - 1, idx(1) - 1), 2)));
  }
  if(type == 1){
    sumnorm = sum(unorm);
    return(sumnorm);
  }else{
    return(unorm);
  }
}

// [[Rcpp::export]]
arma::vec SPGlinear(arma::vec y, arma::mat X_1, arma::mat X_2,
           arma::mat C_1, double C_1norm,
           arma::mat C_2, double C_2norm,
           arma::mat g_idx, arma::vec gamma_est,
           double lambda, double alpha, double Mu, int maxit, double tol,
           bool verbose) {
  // y: response
  // X1: design matrix for un-penalized covariates, including the intercept column
  // X2: design matrix for binary covariates times expand matrix A
  // gamma_est: initial values for coefficients
  // C_1, C_2, C_1norm, C_2norm: intermediate values for regularization
  // g_idx: group structure matrix
  // lambda, alpha: tuning parameter
  // Mu: positive smoothness parameter
  // maxit: maximum number of iterations
  // tol: tolerance for convergencec check
  // verbose: display objective function value

  C_1 = lambda * alpha * C_1;  // generalized lasso type penalty
  C_2 = lambda * (1 - alpha) * C_2;  // group lasso type penalty
  int pall = X_1.n_cols + X_2.n_cols;
  int p1 = X_1.n_cols;


  // set the w value and theta
  arma::vec beta = gamma_est;
  arma::vec w = beta;
  arma::vec w2 = w.rows(p1, pall - 1);
  arma::vec u1_opt = zeros(C_1.n_rows);
  arma::vec u2_opt = zeros(C_2.n_rows);
  arma::vec grad = zeros(gamma_est.n_elem);
  arma::vec grad2 = zeros(gamma_est.n_elem);
  arma::vec beta_new = zeros(gamma_est.n_elem);
  arma::vec beta2_new = zeros(pall - p1);
  arma::vec theta = ones(1);
  arma::vec theta_new = zeros(1);
  arma::vec obj = zeros(maxit);

  // find maximum eigen value of X'X matrix
  arma::mat X = join_rows(X_1, X_2);
  arma::mat xx = trans(X)*X;
  arma::mat xy = trans(X)*y;
  // Lipschitz constant
  double eigenmax = max(eig_sym(xx));
  double L = eigenmax +
    (lambda * alpha * lambda * alpha * C_1norm +
    lambda * (1 - alpha) * lambda * (1 - alpha) * C_2norm) / Mu;

  // start iteration
  for(int iter = 0; iter < maxit; iter++) {
    // compute the optimal solution for u1
    u1_opt = hardThresholdCpp(C_1 * w2 / Mu, 1);
    // compute the optimal solution for u2
    u2_opt = shrinkGroupCpp(C_2 * w2 / Mu, g_idx);
    // gradient
    grad2 = zeros(gamma_est.n_elem);
    grad2.rows(p1, pall - 1) = trans(C_1) * u1_opt + trans(C_2) * u2_opt;
    grad = xx * w - xy + grad2;

    beta_new = w - 1 / L * grad;
    beta2_new = beta_new.rows(p1, pall - 1);
    theta_new = (sqrt(pow(theta, 4) + 4 * pow(theta, 2)) - pow(theta, 2)) / 2;

    // calculate objective function value
    obj.row(iter) = sum(pow(y - X * beta_new, 2)) / 2 + sum(abs(C_1 * beta2_new)) +
      cal2norm(C_2 * beta2_new, g_idx, 1);
    // print out objective function value for debugging
    //if(verbose == 1) {
    //  printf("Iter %i, obj: %f \n", iter+1, obj[iter]);
    //  fflush(stdout);
    //}

    // check convergence
    if(sum(abs(beta_new - w)) / sum(abs(w)) < tol) {
      break;
    }
    //if(iter >= 50){
    //  if(sum(abs((obj.row(iter) - obj.row(iter - 1)) / obj.row(iter - 1))) < tol) {
    //    break;
    //  }
    //}



    // update the w
    w = beta_new +  (beta_new - beta) * as_scalar((1 - theta) / theta * theta_new);
    w2 = w.rows(p1, pall - 1);

    // update the beta and theta
    beta = beta_new;
    theta = theta_new;
  }
  // return the results on gamma
  return(beta_new);
}



// [[Rcpp::export]]
arma::vec SPGlogistic(arma::vec y, arma::mat X_1, arma::mat X_2,
                      arma::mat C_1, double C_1norm,
                      arma::mat C_2, double C_2norm,
                      arma::mat g_idx, arma::vec gamma_est, arma::vec weight,
                      double lambda, double alpha, double Mu, int maxit, double tol,
                      bool verbose) {
  // y: response
  // X1: design matrix for un-penalized covariates, including the intercept column
  // X2: design matrix for binary covariates times expand matrix A
  // gamma_est: initial values for coefficients
  // C_1, C_2, C_1norm, C_2norm: intermediate values for regularization
  // g_idx: group structure matrix
  // weight: weight vector for both classes in logistic regression
  // lambda, alpha: tuning parameter
  // Mu: positive smoothness parameter
  // maxit: maximum number of iterations
  // tol: tolerance for convergencec check
  // verbose: display objective function value

  C_1 = lambda * alpha * C_1;  // generalized lasso type penalty
  C_2 = lambda * (1 - alpha) * C_2;  // group lasso type penalty
  int pall = X_1.n_cols + X_2.n_cols;
  int p1 = X_1.n_cols;


  // set the w value and theta
  arma::vec beta = gamma_est;
  arma::vec w = beta;
  arma::vec w2 = w.rows(p1, pall - 1);
  arma::vec u1_opt = zeros(C_1.n_rows);
  arma::vec u2_opt = zeros(C_2.n_rows);
  arma::vec grad = zeros(gamma_est.n_elem);
  arma::vec grad2 = zeros(gamma_est.n_elem);
  arma::vec v = zeros(y.n_rows);
  arma::vec beta_new = zeros(gamma_est.n_elem);
  arma::vec beta2_new = zeros(pall - p1);
  arma::vec theta = ones(1);
  arma::vec theta_new = zeros(1);
  arma::vec obj = zeros(maxit);

  // find maximum eigen value of X'X matrix
  arma::mat X = join_rows(X_1, X_2);
  arma::mat xx = trans(X)*X;
  arma::mat xy = trans(X)*y;
  // Lipschitz constant
  double eigenmax = max(eig_sym(xx));
  double L = eigenmax / 2 +
    (lambda * alpha * lambda * alpha * C_1norm +
    lambda * (1 - alpha) * lambda * (1 - alpha) * C_2norm) / Mu;

  // start iteration
  for(int iter = 0; iter < maxit; iter++) {
    // compute the optimal solution for u2
    u1_opt = hardThresholdCpp(C_1 * w2 / Mu, 1);
    // compute the optimal solution for u1
    u2_opt = shrinkGroupCpp(C_2 * w2 / Mu, g_idx);
    // gradient
    grad2 = zeros(gamma_est.n_elem);
    grad2.rows(p1, pall - 1) = trans(C_1) * u1_opt + trans(C_2) * u2_opt;
    v = exp(X * w) / (1 + exp(X * w));
    grad = -trans(X) * diagmat(1 - v) * weight[0] * y +
      trans(X) * diagmat(v) * weight[1] * (1-y) + grad2;

    beta_new = w - 1 / L * grad;
    beta2_new = beta_new.rows(p1, pall - 1);
    theta_new = (sqrt(pow(theta, 4) + 4 * pow(theta, 2)) - pow(theta, 2)) / 2;

    // calculate objective function value
    obj.row(iter) = -weight[0] * trans(y) * log(exp(X * beta_new) / (1 + exp(X * beta_new))) +
      weight[1] * trans(1-y) * log(1+exp(X * beta_new)) +
      sum(abs(C_1 * beta2_new)) +
      cal2norm(C_2 * beta2_new, g_idx, 1);
    // print out objective function value for debugging
    //if(verbose == 1) {
    //  printf("Iter %i, obj: %f \n", iter+1, obj[iter]);
    //  fflush(stdout);
    //}

    // check convergence
    if(sum(abs(beta_new - w)) / sum(abs(w)) < tol) {
      break;
    }
    //if(iter >= 50){
    //  if(sum(abs((obj.row(iter) - obj.row(iter - 1)) / obj.row(iter - 1))) < tol) {
    //    break;
    //  }
    //}

    // update the w
    w = beta_new +  (beta_new - beta) * as_scalar((1 - theta) / theta * theta_new);
    w2 = w.rows(p1, pall - 1);

    // update the beta and theta
    beta = beta_new;
    theta = theta_new;
  }
  // return the results on gamma
  return(beta_new);
}


// warm start solver for a sequence of lambda at a fixed alpha value
// [[Rcpp::export]]
arma::mat warm_start_ls(arma::vec y, arma::mat X_1, arma::mat X_2,
                       arma::mat C_1, double C_1norm,
                       arma::mat C_2, double C_2norm,
                       arma::mat g_idx,
                       arma::vec gamma_init,
                       arma::vec lambda, double alpha, double Mu, int maxit, double tol,
                       bool verbose){

  int pall = X_1.n_cols + X_2.n_cols;
  int nlam = lambda.size();
  arma::mat gamma_mat(pall, nlam);
  arma::vec init = gamma_init;
  arma::vec gamma_est = zeros(gamma_init.n_elem);

  for (int k = 0; k < nlam; ++k) {
      gamma_est = SPGlinear(y, X_1, X_2, C_1, C_1norm, C_2, C_2norm,
                            g_idx, init, lambda[k], alpha, Mu, maxit, tol, verbose);
      gamma_mat.col(k) = gamma_est;
      init = gamma_est;
  }
  return(gamma_mat);
}


// warm start solver for a sequence of lambda at a fixed alpha value
// [[Rcpp::export]]
arma::mat warm_start_logistic(arma::vec y, arma::mat X_1, arma::mat X_2,
                        arma::mat C_1, double C_1norm,
                        arma::mat C_2, double C_2norm,
                        arma::mat g_idx,
                        arma::vec gamma_init, arma::vec weight,
                        arma::vec lambda, double alpha, double Mu, int maxit, double tol,
                        bool verbose){

  int pall = X_1.n_cols + X_2.n_cols;
  int nlam = lambda.size();
  arma::mat gamma_mat(pall, nlam);
  arma::vec init = gamma_init;
  arma::vec gamma_est = zeros(gamma_init.n_elem);

  for (int k = 0; k < nlam; ++k) {
    gamma_est = SPGlogistic(y, X_1, X_2, C_1, C_1norm, C_2, C_2norm,
                          g_idx, init, weight, lambda[k], alpha, Mu, maxit, tol, verbose);
    gamma_mat.col(k) = gamma_est;
    init = gamma_est;
  }
  return(gamma_mat);
}


