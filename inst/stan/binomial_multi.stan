data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  int<lower = 0, upper=1> y_obs[N_obs];
  int<lower = 2> k;
  int<lower=0> m[k];
  int<lower=0> max_m;
  int<lower = 1, upper = max_m> ind[N_obs, k];
  int<lower = 1, upper = max_m> missing_ind[N_mis, k];
  matrix<lower = 0>[max_m-1, k] diff;
  int<lower=0> p;
  matrix[N_obs, p] X;
  int<lower=0, upper=1> has_covs;
  real<lower = 1> local_dof_stan;
  real<lower = 1> global_dof_stan; // dof of pi(tau), = 1
  real<lower = 0> alpha_scale_stan; // c, Section 2// dof of pi(lambda), = 1
  int<lower=0,upper=1> is_monotone;
  int<lower=0,upper=1> exp_approach;
  int<lower=0,upper=k> q;
  int<lower=1,upper=k> monotone_inds[q];
  int<lower=1,upper=k> nonmonotone_inds[k-q];
  real samp_mean;
  real samp_sd;
  vector[p] beta_mean;
  vector[p] beta_sd;
  vector<lower=1.0>[p] beta_df;
  int<lower=0,upper=1> cauchy_beta;
  int<lower=0> N_new;
  int<lower = 1, upper = max_m> new_ind[N_new, k];
  matrix[N_new, p] new_X;
}

transformed data {
    matrix[max_m-1, k] sqrt_diff = sqrt(diff);
}

parameters {
  real alpha;
  vector[p] beta_sub1;
  vector<lower = 0.0>[p] beta_sub2;
  vector<lower = 0.0>[k] tau_glob1;
  vector<lower = 0.0>[k] tau_glob2;
  matrix[max_m-1,k] gamma_aux;
  matrix<lower = 0.0>[max_m-1,k] lambda_loc1;
  matrix<lower = 0.0>[max_m-1,k] lambda_loc2;
}

transformed parameters {
  vector[k] tau_glob;
  matrix[max_m-1,k] lambda;
  matrix[max_m,k] eta;
  matrix[max_m,k] path;
  matrix[N_obs,k] indiv_path;
  matrix[N_mis,k] miss_path;
  vector[N_obs] finalpath;
  vector[N_obs] f;
  vector[p] beta;
    if (cauchy_beta){
      beta = ((5 * beta_sd .* beta_sub1) ./ sqrt(beta_sub2)) + beta_mean;
    } else{
      beta = (5 * beta_sd .* beta_sub1) + beta_mean;
    }
    eta[1,] = rep_row_vector(0,k);
    lambda = lambda_loc1 .* sqrt(lambda_loc2);
    tau_glob = tau_glob1 .* sqrt(tau_glob2) * alpha_scale_stan;
    if (is_monotone){
      if (exp_approach){
        for (i in 1:q){
        eta[2:m[monotone_inds[i]], monotone_inds[i]] = exp(tau_glob[monotone_inds[i]] * lambda[,monotone_inds[i]] .* gamma_aux[,monotone_inds[i]] .* sqrt_diff[,monotone_inds[i]]);
        eta[m[monotone_inds[i]]+1:max_m, monotone_inds[i]] = rep_vector(0, max_m-monotone_inds[i]);
      } for (j in 1:(k-q)){
        eta[2:m[nonmonotone_inds[j]], nonmonotone_inds[j]] = tau_glob[nonmonotone_inds[j]] * lambda[,nonmonotone_inds[j]] .* gamma_aux[,nonmonotone_inds[j]] .* sqrt_diff[,nonmonotone_inds[j]];
        eta[m[nonmonotone_inds[j]]+1:max_m, nonmonotone_inds[j]] = rep_vector(0, max_m-nonmonotone_inds[j]);
      }
      } else{
        for (i in 1:q){
        eta[2:m[monotone_inds[i]], monotone_inds[i]] = fabs(tau_glob[monotone_inds[i]] * lambda[,monotone_inds[i]] .* gamma_aux[,monotone_inds[i]] .* sqrt_diff[,monotone_inds[i]]);
        eta[m[monotone_inds[i]]+1:max_m, monotone_inds[i]] = rep_vector(0, max_m-monotone_inds[i]);
      } for (j in 1:(k-q)){
        eta[2:m[nonmonotone_inds[j]], nonmonotone_inds[j]] = tau_glob[nonmonotone_inds[j]] * lambda[,nonmonotone_inds[j]] .* gamma_aux[,nonmonotone_inds[j]] .* sqrt_diff[,nonmonotone_inds[j]];
        eta[m[nonmonotone_inds[j]]+1:max_m, nonmonotone_inds[j]] = rep_vector(0, max_m-nonmonotone_inds[j]);
      }
      }
    } else{
      for (i in 1:k){
        eta[2:m[i],i] = tau_glob[i] * lambda[,i] .* gamma_aux[,i] .* sqrt_diff[,i];
      }
    }
    for (i in 1:k){
      path[,i] = cumulative_sum(eta[,i]);
      indiv_path[,i] = path[ind[,i],i];
      if (N_mis > 0){
        miss_path[,i] = path[missing_ind[,i],i];
      }
    }
    finalpath = alpha + indiv_path*rep_vector(1,k);
    if (has_covs==0){
      f = finalpath;
    } else{
      f = X*beta + finalpath;
    }
}

model {
    alpha ~ normal(samp_mean, 5*samp_sd);
    beta_sub1 ~ std_normal();
    beta_sub2 ~ gamma(0.5*beta_df, 0.5*beta_df);
    tau_glob1 ~ std_normal();
    tau_glob2 ~ inv_gamma(0.5*global_dof_stan, 0.5*global_dof_stan);
    for (i in 1:k){
      lambda_loc1[,i] ~ std_normal();
      lambda_loc2[,i] ~ inv_gamma(0.5*local_dof_stan, 0.5*local_dof_stan);
      gamma_aux[,i] ~ std_normal();
    }
    y_obs ~ bernoulli_logit(f);
}
