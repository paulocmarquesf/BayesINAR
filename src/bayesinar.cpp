// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace std;
using namespace Rcpp;

void set_seed(unsigned int seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);
}

double lfactorial(double x) {
    return lgamma(x + 1);
}

double rtbeta(double a0, double b0, double max) {
    double upper_limit = R::pbeta(max, a0, b0, TRUE, FALSE);
    double u = R::runif(0, upper_limit);
    return R::qbeta(u, a0, b0, TRUE, FALSE);
}

IntegerVector table_each_count(NumericVector x) {
    int t = x.size(); 
    IntegerVector n(t); 
    std::map<int, int> cnt;
    for (int j = 0; j < t; j++) cnt[x[j]]++;
    for (int j = 0; j < t; j++) n[j] = cnt[x[j]];
    return n;
}

IntegerVector minus_m_ij(IntegerMatrix x, int i, int j) {
    int n = x.ncol();
    IntegerVector minus_matrix(n-1);
    for (int t = 0; t < n; t++) {
        if (t < j) minus_matrix[t] =  x(i, t); 
        else if (t > j) minus_matrix[t - 1] = x(i, t); 
    }
    return minus_matrix; 
}

NumericVector minus_alpha_ij(NumericMatrix x, int i, int j) {
    int n = x.ncol();
    NumericVector minus_vector(n - 1);
    for (int t = 0; t < n; t++) {
        if (t < j) minus_vector[t] =  x(i, t) ; 
        else if (t > j) minus_vector[t - 1] = x(i - 1, t) ; 
    }
    return minus_vector; 
}

// [[Rcpp::export(.generalized_median)]]
int generalized_median(NumericVector pred) {
    int median = 0;
    double cum_sum = pred[0];
    double least_abs_diff = abs(0.5 - cum_sum);
    for (int i = 1; i < pred.size(); i++) {
        cum_sum += pred[i];
        double abs_diff = abs(0.5 - cum_sum);
        if (abs_diff < least_abs_diff) {
            median = i;
            least_abs_diff = abs_diff;
        }
    }
    return median;
}

// INAR functions

int m_full_conditional_inar(double alpha, double lambda, int yt_i, int yt, IntegerVector m_minus_t) {        
    int M = m_minus_t.size(), sum_m_minus_t; 
    
    sum_m_minus_t = 0; for(int j = 0; j < M; j++) sum_m_minus_t += m_minus_t[j]; 
    
    if (yt_i == 0 || yt - sum_m_minus_t == 0) return 0;

    NumericVector log_weights(min(yt_i, yt - sum_m_minus_t) + 1), weights(min(yt_i, yt - sum_m_minus_t) + 1);
    int m_t;
    double total = 0;

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        log_weights[m_t] = m_t*(log(alpha) - log(lambda) - log(1 - alpha))
                            - lfactorial(m_t) - lfactorial(yt_i - m_t) - lfactorial(yt - m_t - sum_m_minus_t);
    }
    
    double max_log_weights = max(log_weights);
    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if (R_finite(log_weights[m_t])) weights[m_t] = exp(log_weights[m_t] - max_log_weights);
        else                            weights[m_t] = 0;
        total += weights[m_t];
    }

  double u = R::runif(0, 1);

  for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
      if ((u -= (weights[m_t] / total)) <= 0) return m_t;
  }
}

// [[Rcpp::export(.posterior_inar)]]
List posterior_inar(IntegerVector y, 
                    int p, 
                    List prior, 
                    int burn_in, 
                    int N, 
                    unsigned int random_seed,
                    bool verbose)
{
    int T = y.size(), sum_m;
    double sum_minus_alpha;
    NumericVector lambda(burn_in + N), minus_alpha(p - 1), sum1(p), sum2(p);
    NumericMatrix alpha(burn_in + N, p);
    IntegerVector minus_m(T - p - 1);
    IntegerMatrix m(T - p, p);
    Progress pb(burn_in + N + 1, true);
    List model;

    set_seed(random_seed); 
    
    double a_alpha = prior["a_alpha"]; 
    double b_alpha = prior["b_alpha"]; 
    double a_lambda = prior["a_lambda"]; 
    double b_lambda = prior["b_lambda"]; 
    
    for (int j = 0; j < p; j++) alpha(0, j) = 0;
    lambda[0] = 1; 
    
    pb.increment();
    for (int i = 1; i < burn_in + N; i++) {
        if (Progress::check_abort()) return model;

        for(int j = 0; j < p; j++) {
            sum1[j] = 0; for (int t = 0; t < T-p; t++) sum1[j] += m(t, j);
            sum2[j] = 0; for (int t = 0; t < T-p; t++) sum2[j] += y[t + p - j - 1] - m(t, j);
            minus_alpha = minus_alpha_ij(alpha, i, j); 
            sum_minus_alpha = 0; for (int ell = 0; ell < p-1; ell++) sum_minus_alpha += minus_alpha[ell];
            alpha(i, j) = rtbeta(a_alpha + sum1[j], b_alpha + sum2[j], 1 - sum_minus_alpha);        
        }
        
        double sum_lambda = 0; 
        for (int t = 0; t < T-p; t++){
            sum_m = 0; for(int j = 0; j < p; j++)  sum_m += m(t, j); 
            sum_lambda += y[t+p] - sum_m; 
        }
        lambda[i] = R::rgamma(a_lambda + sum_lambda, 1 / (b_lambda + T- p));
        
        for(int j = 0; j < p; j++){
            for (int t = 0; t < T-p; t++) {
                minus_m = minus_m_ij(m, t, j); 
                m(t, j) = m_full_conditional_inar(alpha(i, j), lambda[i], y[t + p - j - 1], y[t + p], minus_m);
            }
        }
        
        pb.increment();
    }
    
    model["time_series"] = y; 
    model["p"] = p; 
    
    List burn_in_pars; 
    burn_in_pars["alpha"] = alpha(Range(0, burn_in - 1), _);
    burn_in_pars["lambda"] = lambda[Range(0, burn_in - 1)]; 
    model["burn_in"] = burn_in_pars; 

    List chain_pars; 
    chain_pars["alpha"] = alpha(Range(burn_in, burn_in + N - 1), _);
    chain_pars["lambda"] = lambda[Range(burn_in, burn_in + N - 1)]; 
    model["chain"] = chain_pars; 
    
    return model;
}

// [[Rcpp::export(.predictive_distribution_inar)]]
NumericVector predictive_distribution_inar(List model, int h, int replications) {
    IntegerVector y = model["time_series"];
    int T = y.size();
    int p = 1;
    List chain = model["chain"];
    NumericMatrix alpha = chain["alpha"];
    NumericVector lambda = chain["lambda"];
    int N = chain["length"];
    
    int y_next = 0;
    IntegerVector y_prev(T + h); 
    IntegerVector yt_plus_h(N * replications); 
    
    for(int i = 0; i < T; i++) y_prev[i] = y[i]; 
    
    for (int i = 0; i < N; i++) {
        for (int ell = 0; ell < replications; ell++) { 
            for (int j = 1; j <= h; j++) {
                y_next = 0;
                for(int r = 0; r < p; r++) y_next = R::rbinom(y_prev[T - p - 1 + r + j], alpha(i, p - r - 1) );
                y_next += R::rpois(lambda[i]);
                y_prev[T + j - 1] = y_next; 
            }
            yt_plus_h[i*replications + ell] = y_next;
        }
    }
    
    NumericVector count = Rcpp::as<NumericVector>(table(yt_plus_h));
    int c = count.size();
    NumericVector pred(c);
    for (int i = 0; i < c; i++) pred[i] = count[i] / (N * replications);
    
    return pred;
}

// AdINAR functions

int m_full_conditional_adinar(double alpha, double theta, double lambda, int yt_i, int yt, int u_t, IntegerVector m_minus_t) {
    int M = m_minus_t.size();
    
    int sum_m_minus_t = 0;
    for(int j = 0; j < M; j++) sum_m_minus_t += m_minus_t[j]; 
    
    if (yt_i == 0 || yt - sum_m_minus_t == 0) return 0;

    NumericVector log_weights(min(yt_i, yt - sum_m_minus_t) + 1), weights(min(yt_i, yt - sum_m_minus_t) + 1);
    int m_t;
    double total = 0;

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if (u_t == 1) // Geometric
            log_weights[m_t] = m_t*(log(alpha) - log(1 - theta) - log(1 - alpha)) - lfactorial(m_t) - lfactorial(yt_i - m_t);
        else          // Poisson
            log_weights[m_t] = m_t*(log(alpha) - log(lambda) - log(1 - alpha)) - lfactorial(m_t) - lfactorial(yt_i - m_t) - lfactorial(yt - m_t - sum_m_minus_t);
    }

    double max_log_weights = max(log_weights);
    
    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if (R_finite(log_weights[m_t])) weights[m_t] = exp(log_weights[m_t] - max_log_weights);
        else                            weights[m_t] = 0;
        total += weights[m_t];
    }

    double u = R::runif(0, 1);

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if ((u -= (weights[m_t] / total)) <= 0) return m_t;
    }
}

// [[Rcpp::export(.posterior_adinar)]]
List posterior_adinar(IntegerVector y, 
                      int p, 
                      List prior, 
                      int burn_in, 
                      int N, 
                      unsigned int random_seed,
                      bool verbose) 
{
    int T = y.size(), sum_m;
    double sum_minus_alpha;
    NumericVector theta(burn_in + N), lambda(burn_in + N), w(burn_in + N), minus_alpha(p - 1), sum1(p), sum2(p);
    IntegerVector minus_m(T - p - 1);
    NumericMatrix alpha(burn_in + N, p);
    IntegerMatrix u(burn_in + N, T - p), m(T - p, p);
    Progress pb(burn_in + N + 1, true);
    List model;

    set_seed(random_seed); 
    
    double a_alpha = prior["a_alpha"]; 
    double b_alpha = prior["b_alpha"]; 
    double a_lambda = prior["a_lambda"]; 
    double b_lambda = prior["b_lambda"]; 
    double a_theta = prior["a_theta"]; 
    double b_theta = prior["b_theta"]; 
    double a_w = prior["a_w"];
    double b_w = prior["b_w"];
    
    for (int j = 0; j < p; j++) alpha(0, j) = 0;
    lambda[0] = 1; theta[0] = 0; w[0] = 0.5;
    
    pb.increment(); 
    for (int i = 1; i < burn_in + N; i++) {
        if (Progress::check_abort()) return model;
        
        for(int j = 0; j < p; j++) {
            sum1[j] = 0; for (int t = 0; t < T-p; t++) sum1[j] += m(t, j);
            sum2[j] = 0; for (int t = 0; t < T-p; t++) sum2[j] += y[t + p - j - 1] - m(t, j);
            minus_alpha = minus_alpha_ij(alpha, i, j); 
            sum_minus_alpha = 0; for (int ell = 0; ell < p-1; ell++) sum_minus_alpha += minus_alpha[ell];
            alpha(i, j) = rtbeta(a_alpha + sum1[j], b_alpha + sum2[j], 1 - sum_minus_alpha);        
        }
        
        int sum_ut = 0;
        double sum_lambda = 0, sum_theta = 0;
        for (int t = 0; t < T - p; t++) {
            sum_m = 0; for(int j = 0; j < p; j++)  sum_m += m(t, j); 
            double w_geo = exp(log(w[i-1]) + log(theta[i-1]) + (y[t + p] - sum_m)*log(1 - theta[i-1]));
            double w_poi = exp(log(1 - w[i-1]) - lambda[i-1] + (y[t + p] - sum_m)*log(lambda[i-1]) - lfactorial(y[t + p] - sum_m));
            double pr_ut_equal_one = w_geo / (w_geo + w_poi);
            u(i, t) = R::rbinom(1, pr_ut_equal_one);
            if (u(i, t) == 1) sum_theta += y[t + p] - sum_m;
            else sum_lambda += y[t + p] - sum_m;
            sum_ut += u(i, t);
        }
        
        theta[i] = R::rbeta(a_theta + sum_ut, b_theta + sum_theta);
        lambda[i] = R::rgamma(a_lambda + sum_lambda, 1 / (b_lambda + T - p - sum_ut));
        w[i] = R::rbeta(a_w + sum_ut, b_w + T - p - sum_ut);
        
        for(int j = 0; j < p; j++){
            for (int t = 0; t < T - p; t++) {
                minus_m = minus_m_ij(m, t, j); 
                m(t, j) = m_full_conditional_adinar(alpha(i, j), theta[i], lambda[i], y[t + p - j - 1], y[t + p], u(i, t), minus_m);
            }
        }

        pb.increment();
    }

    model["time_series"] = y; 
    model["p"] = p; 
    
    List burn_in_pars; 
    burn_in_pars["alpha"] = alpha(Range(0, burn_in - 1), _);
    burn_in_pars["lambda"] = lambda[Range(0, burn_in - 1)]; 
    burn_in_pars["theta"] = theta[Range(0, burn_in - 1)]; 
    burn_in_pars["w"] = w[Range(0, burn_in - 1)]; 
    model["burn_in"] = burn_in_pars; 

    List chain_pars; 
    chain_pars["alpha"] = alpha(Range(burn_in, burn_in + N - 1), _);
    chain_pars["lambda"] = lambda[Range(burn_in, burn_in + N - 1)]; 
    chain_pars["theta"] = theta[Range(burn_in, burn_in + N - 1)]; 
    chain_pars["w"] = w[Range(burn_in, burn_in + N - 1)]; 
    
    model["chain"] = chain_pars; 
    
    return model;
}

// [[Rcpp::export(.predictive_distribution_adinar)]]
NumericVector predictive_distribution_adinar(List model, int h, int replications) {
    IntegerVector y = model["time_series"];
    int T = y.size();
    int p = 1;
    List chain = model["chain"];
    NumericMatrix alpha = chain["alpha"];
    NumericVector lambda = chain["lambda"];
    NumericVector theta = chain["theta"];
    NumericVector w = chain["w"];
    int N = chain["length"];
    
    int y_next = 0;
    IntegerVector y_prev(T + h); 
    IntegerVector yt_plus_h(N * replications); 
    
    for(int i = 0; i < T; i++) y_prev[i] = y[i]; 
    
    for (int i = 0; i < N; i++) {
        for (int ell = 0; ell < replications; ell++) { 
            for (int j = 1; j <= h; j++) {
                y_next = 0;
                for(int r = 0; r < p; r++) y_next = R::rbinom(y_prev[T - p - 1 + r + j], alpha(i, p - r - 1) );
                double u = R::runif(0, 1);
                if (u <= w[i]) y_next += R::rgeom(theta[i]);
                else           y_next += R::rpois(lambda[i]);
                y_prev[T + j - 1] = y_next; 
            }
            yt_plus_h[i*replications + ell] = y_next;
        }
    }
    
    NumericVector count = Rcpp::as<NumericVector>(table(yt_plus_h));
    int c = count.size();
    NumericVector pred(c);
    for (int i = 0; i < c; i++) pred[i] = count[i] / (N * replications);
    
    return pred;
}

// DP-INAR functions

int m_full_conditional_dpinar(double alpha, double lambda, int yt_i, int yt, IntegerVector m_minus_t) {
    int M = m_minus_t.size(), sum_m_minus_t;

    sum_m_minus_t = 0; for(int j = 0; j < M; j++) sum_m_minus_t += m_minus_t[j];

    if (yt_i == 0 || yt - sum_m_minus_t == 0) return 0;

    NumericVector log_weights(min(yt_i, yt - sum_m_minus_t) + 1), weights(min(yt_i, yt - sum_m_minus_t) + 1);
    int m_t;
    double total = 0;

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        log_weights[m_t] = m_t*(log(alpha) - log(lambda) - log(1 - alpha))
                           - lfactorial(m_t) - lfactorial(yt_i - m_t) - lfactorial(yt - m_t - sum_m_minus_t);
    }

    double max_log_weights = max(log_weights);

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if (R_finite(log_weights[m_t])) weights[m_t] = exp(log_weights[m_t] - max_log_weights);
        else                            weights[m_t] = 0;
        total += weights[m_t];
    }

    double u = R::runif(0, 1);

    for (m_t = 0; m_t <= min(yt_i, yt - sum_m_minus_t); m_t++) {
        if ((u -= (weights[m_t] / total)) <= 0) return m_t;
    }

    stop("m_full_conditional: should never get to this point!\n");
}

double lambda_full_conditional_dpinar(NumericVector minus_lambda_t, double a0, double b0, double tau, int y_t, int sum_m) {
    int r, J = minus_lambda_t.size(), j;
    NumericVector log_weights(J), weights(J);
    double weight_gamma, total = 0;

    for (j = 0; j < J; j++) {
        log_weights[j] = - minus_lambda_t[j] + (y_t - sum_m)*log(minus_lambda_t[j]);
    }

    weight_gamma = exp(log(tau) + a0*log(b0) + lgamma(y_t - sum_m + a0) - lgamma(a0) - (y_t - sum_m + a0)*log(b0 + 1));
    total += weight_gamma;

    for (r = 0; r < J; r++) {
        if (R_finite(log_weights[r])) weights[r] = exp(log_weights[r]);
        else                          weights[r] = 0;
        total += weights[r];
     }

     double u = R::runif(0, 1);

     if ((u -= (weight_gamma / total)) <= 0) {
         return R::rgamma(a0 + y_t - sum_m, 1 / (b0 + 1) );
     } else {
          for (r = 0; r < J; r++) {
              if ((u -= (weights[r] / total)) <= 0) return minus_lambda_t[r];
         }
     }

     stop("lambda_full_conditional: should never get to this point!\n");
}

double lambda_star_full_conditional(int j, IntegerVector y, IntegerMatrix m, IntegerVector c, double a0, double b0) {
    int sum = 0, n_j = 0, p = m.ncol(), T = y.size();

    for (int t = p; t < T; t++) {
        int sum_m = 0;
        if (c[t - p] == j + 1) {
            for (int r = 0; r < p; r++) sum_m += m(t - p , r);
            sum += y[t] - sum_m;
            n_j++;
        }
    }

    return R::rgamma(a0 + sum, 1 / (b0 + n_j) );
}

double tau_full_conditional(int m, int k, double a_tau, double b_tau, double u) {
    double tau;

    double w1 = exp( lgamma(a_tau + k) - (a_tau + k)*log(b_tau - log(u)) );
    double w2 = exp( log(m) + lgamma(a_tau + k - 1) - (a_tau + k - 1)*log(b_tau - log(u)) );

    double v = R::runif(0, 1);

    if (v <= (w1 / (w1 + w2))) tau = R::rgamma(a_tau + k, 1 / (b_tau - log(u)));
    else                       tau = R::rgamma(a_tau + k - 1, 1 / (b_tau - log(u)));

    return tau;
}

// [[Rcpp::export(.posterior_dpinar)]]
List posterior_dpinar(IntegerVector y,
                      int p,
                      List prior,
                      int burn_in,
                      int N,
                      unsigned int random_seed,
                      bool verbose)
{
    int T = y.size(), sum_m;
    double sum_minus_alpha;
    NumericVector tau(burn_in + N), minus_lambda_t(T - p - 1), minus_alpha(p - 1), sum_1(p), sum_2(p);
    IntegerVector clusters(burn_in + N), minus_m(T - p - 1);
    NumericMatrix alpha(burn_in + N, p), lambda(burn_in + N, T - p);
    IntegerMatrix m(T - p, p);
    Progress pb(burn_in + N + 1, true);
    List model;

    set_seed(random_seed);

    double a_alpha = prior["a_alpha"];
    double b_alpha = prior["b_alpha"]; 
    double a_tau = prior["a_tau"];
    double b_tau = prior["b_tau"];
    double a0 = prior["a0"];
    double b0 = prior["b0"];

    for (int j = 0; j < p; j++) alpha(0, j) = 0;
    for (int t = 0; t < T - p; t++) lambda(0, t) = 1;
    tau[0] = 1;

    pb.increment();
    for (int i = 1; i < burn_in + N; i++) {
        if (Progress::check_abort()) return model;

        for (int j = 0; j < p; j++) {
            sum_1[j] = 0; for (int t = 0; t < T - p; t++) sum_1[j] += m(t, j);
            sum_2[j] = 0; for (int t = 0; t < T - p; t++) sum_2[j] += y[t + p - j - 1] - m(t, j);
            minus_alpha = minus_alpha_ij(alpha, i, j);
            sum_minus_alpha = 0; for (int ell = 0; ell < p - 1; ell++) sum_minus_alpha += minus_alpha[ell];
            alpha(i, j) = rtbeta(a_alpha + sum_1[j], b_alpha + sum_2[j], 1 - sum_minus_alpha);
        }
        for (int t = 0; t < T - p; t++) {
            for (int s = 0; s < T - p; s++) {
                if (s < t) minus_lambda_t[s] = lambda(i, s);
                else if (s > t) minus_lambda_t[s - 1] = lambda(i - 1, s);
            }
            sum_m = 0; for(int j = 0; j < p; j++)  sum_m += m(t, j);
            lambda(i, t) = lambda_full_conditional_dpinar(minus_lambda_t, a0, b0, tau[i - 1], y[t + p], sum_m);
        }

        for (int j = 0; j < p; j++) {
            for (int t = 0; t < T - p; t++) {
                minus_m = minus_m_ij(m, t, j);
                m(t, j) = m_full_conditional_dpinar(alpha(i, j), lambda(i, t), y[t + p - j - 1], y[t + p], minus_m);
            }
        }
        NumericVector lambda_star = unique(lambda.row(i));
        int k = lambda_star.size();
        clusters[i] = k;
        IntegerVector c = match(lambda.row(i), lambda_star);

        for (int j = 0; j < k; j++) lambda_star[j] = lambda_star_full_conditional(j, y, m, c, a0, b0);

        for (int t = 0; t < T-p; t++) lambda(i, t) = lambda_star[c[t] - 1];

        double u = R::rbeta(tau[i - 1] + 1, T - p);
        tau[i] = tau_full_conditional(T - p, k, a_tau, b_tau, u);

        pb.increment();
    }

    model["time_series"] = y;
    model["p"] = p;

    List burn_in_pars;
    burn_in_pars["alpha"] = alpha(Range(0, burn_in - 1), _);
    burn_in_pars["lambda"] = lambda(Range(0, burn_in - 1), _);
    burn_in_pars["tau"] = tau[Range(0, burn_in - 1)];
    burn_in_pars["num_clusters"] = clusters[Range(0, burn_in - 1)];
    model["burn_in"] = burn_in_pars;

    List chain_pars;
    chain_pars["alpha"] = alpha(Range(burn_in, burn_in + N - 1), _);
    chain_pars["lambda"] = lambda(Range(burn_in, burn_in + N - 1), _);
    chain_pars["tau"] = tau[Range(burn_in, burn_in + N - 1)];
    chain_pars["num_clusters"] = clusters[Range(burn_in, burn_in + N - 1)];
    model["chain"] = chain_pars;

    return model;
}

double polya_blackwell_macqueen(double a0, double b0, double tau, NumericVector lambda) {
    int t = lambda.size();
    double u = R::runif(0, 1);
    if ((u -= (tau / (tau + t))) <= 0) return R::rgamma(a0, 1 / b0);
    for (int r = 0; r < t - 1; r++) if ((u -= (1 / (tau + t))) <= 0) return lambda[r];
    return lambda[t - 1];
}

// [[Rcpp::export(.predictive_distribution_dpinar)]]
NumericVector predictive_distribution_dpinar(List model, int h, int replications) {
    IntegerVector y = model["time_series"];
    int T = y.size();
    int p = 1;
    List chain = model["chain"];
    NumericMatrix alpha = chain["alpha"];
    NumericMatrix lambda_prev = chain["lambda"];
    NumericVector tau = chain["tau"];
    int N = chain["length"];
    List prior = model["prior"];
    double a0 = prior["a0"], b0 = prior["b0"];

    int y_next = 0;
    NumericMatrix lambda(N, T - p + h);
    IntegerVector y_prev(T + h);
    IntegerVector yt_plus_h(N * replications);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < T - p; j++) lambda(i, j) = lambda_prev(i, j);
        // Successive Polya-Blackwell-MacQueen urns
        for (int r = 1; r <= h; r++) {
            NumericVector given_lambda(T - p - 1 + r);
            for (int s = 0; s < T - p - 1 + r; s++) given_lambda[s] = lambda(i, s);
            lambda(i, T - p - 1 + r) = polya_blackwell_macqueen(a0, b0, tau[i], given_lambda);
        }
    }

    for(int j = 0; j < T; j++) y_prev[j] = y[j];

    for (int i = 0; i < N; i++) {
        for (int ell = 0; ell < replications; ell++) {
            for (int j = 1; j <= h; j++) {
                y_next = 0;
                for(int r = 0; r < p; r++) y_next += R::rbinom(y_prev[T - p - 1 + r + j], alpha(i, p - r - 1));
                y_next += R::rpois(lambda(i, T - p - 1 + j));
                y_prev[T + j - 1] = y_next;
            }
            yt_plus_h[i*replications + ell] = y_next;
        }
    }

    NumericVector count = Rcpp::as<NumericVector>(table(yt_plus_h));
    int c = count.size();
    NumericVector pred(c);
    for (int i = 0; i < c; i++) pred[i] = count[i] / (N * replications);

    return pred;
}
