data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  real y_obs[N_obs];
  int<lower=0> m;
  int<lower = 1, upper = m> ind[N_obs];
  int<lower = 1, upper = m> missing_ind[N_mis];
  vector<lower = 0>[m-1] diff;
  int<lower=0> p;
  matrix[N_obs, p] X;
  int<lower=0, upper=1> has_covs;
  real<lower = 1> local_dof_stan;
  real<lower = 1> global_dof_stan; // dof of pi(tau), = 1
  real<lower = 0> alpha_scale_stan; // c, Section 2// dof of pi(lambda), = 1
  int<lower=0,upper=1> is_monotone;
  int<lower=0,upper=1> exp_approach;
  real samp_mean;
  real samp_sd;
  vector[p] beta_mean;
  vector[p] beta_sd;
  vector<lower=1.0>[p] beta_df;
  int<lower=0,upper=1> cauchy_beta;
  int<lower=0> N_new;
  int<lower=1, upper=m> new_ind[N_new];
  matrix[N_new, p] new_X;
}

transformed data {
    vector[m-1] sqrt_diff = sqrt(diff);
}

parameters {
  real alpha;
  vector[p] beta_sub1;
  vector<lower = 0.0>[p] beta_sub2;
  real<lower = 0.0> tau_glob1;
  real<lower = 0.0> tau_glob2;
  real<lower=0> sig1;
  real<lower = 0.0> sig2;
  vector[m-1] gamma_aux;
  vector<lower = 0.0>[m-1] lambda_loc1;
  vector<lower = 0.0>[m-1] lambda_loc2;
}

transformed parameters {
  real tau_glob;
  real sig;
  vector[m-1] lambda;
  vector[m] eta;
  vector[m] path;
  vector[N_obs] f;
  vector[p] beta;
    if (cauchy_beta){
      beta = ((5 * beta_sd .* beta_sub1) ./ sqrt(beta_sub2)) + beta_mean;
    } else{
      beta = (5 * beta_sd .* beta_sub1) + beta_mean;
    }
    eta[1] = 0;
    lambda = lambda_loc1 .* sqrt(lambda_loc2);
    tau_glob = tau_glob1 * sqrt(tau_glob2) * alpha_scale_stan;
    sig = sig1 * sqrt(sig2) * 5;
    if (is_monotone){
      if (exp_approach){
        eta[2:m] = exp(tau_glob * lambda .* gamma_aux .* sqrt_diff);
      } else{
        eta[2:m] = fabs(tau_glob * lambda .* gamma_aux .* sqrt_diff);
      }
    } else{
      eta[2:m] = tau_glob * lambda .* gamma_aux .* sqrt_diff;
    }
    path = alpha + cumulative_sum(eta);
    if (has_covs==0){
      f = path[ind];
    } else{
      f = X*beta + path[ind];
    }
}

model {
    alpha ~ normal(samp_mean, 5*samp_sd);
    beta_sub1 ~ std_normal();
    beta_sub2 ~ gamma(0.5*beta_df, 0.5*beta_df);
    tau_glob1 ~ std_normal();
    tau_glob2 ~ inv_gamma(0.5*global_dof_stan, 0.5*global_dof_stan);
    lambda_loc1 ~ std_normal();
    lambda_loc2 ~ inv_gamma(0.5*local_dof_stan, 0.5*local_dof_stan);
    sig1 ~ std_normal();
    sig2 ~ inv_gamma(0.5, 0.5);
    gamma_aux ~ std_normal();
    y_obs ~ normal(f, sig);
}

generated quantities{
  real new_pred[N_new];
  if (has_covs==0){
    new_pred = normal_rng(path[new_ind], sig);
  } else{
    new_pred = normal_rng(path[new_ind] + new_X*beta, sig);
  }
}
