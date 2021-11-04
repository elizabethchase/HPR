data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  real y_obs[N_obs];
  int<lower=0> m;
  int<lower = 1, upper = m> ind[N_obs];
  int<lower = 1, upper = m-1> nearest[N_mis];
  vector<lower = 0>[N_mis] miss_diff;
  vector<lower = 0>[m-1] diff;
  int<lower=0> p;
  matrix[N_obs, p] X;
  int<lower=0, upper=1> has_covs;
  real<lower = 1> local_dof_stan;
  real<lower = 1> global_dof_stan; // dof of pi(tau), = 1
  real<lower = 0> alpha_scale_stan; // c, Section 2// dof of pi(lambda), = 1
  real samp_mean;
  real samp_sd;
  int<lower = 0, upper=1> use_mean;
}

transformed data {
    vector[m-1] sqrt_diff = sqrt(diff);
    vector[N_mis] sqrt_miss_diff = sqrt(miss_diff);
}

parameters {
  real alpha;
  vector[p] beta;
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
  vector[m] f;
    eta[1] = 0;
    lambda = lambda_loc1 .* sqrt(lambda_loc2);
    tau_glob = tau_glob1 * sqrt(tau_glob2) * alpha_scale_stan;
    sig = sig1 * sqrt(sig2) * 5;
    eta[2:m] = tau_glob * lambda .* gamma_aux .* sqrt_diff;
    path = cumulative_sum(eta);
    if (has_covs==0){
      f = alpha + path[ind];
    } else{
      f = alpha + X*beta + path[ind];
    }
}

model {
    alpha ~ normal(samp_mean, 5*samp_sd);
    beta ~ normal(0, 5);
    tau_glob1 ~ std_normal();
    tau_glob2 ~ inv_gamma(0.5*global_dof_stan, 0.5*global_dof_stan);
    lambda_loc1 ~ std_normal();
    lambda_loc2 ~ inv_gamma(0.5*local_dof_stan, 0.5*local_dof_stan);
    sig1 ~ std_normal();
    sig2 ~ inv_gamma(0.5, 0.5);
    gamma_aux ~ std_normal();
    y_obs ~ normal(f, sig);
}

generated quantities {
  vector[m] var_vec_obs;
  vector[N_mis] var_vec_mis;
  vector[N_mis] var_mis;
  vector[N_mis] f2;
  real delta = 1e-9;
  matrix[N_obs, N_obs] L_K;
  vector[N_obs] K_div_y1;
  matrix[N_obs, N_mis] k_x1_x2;
  matrix[N_obs, N_mis] v_pred;
  vector[N_mis] f2_mu;
  matrix[N_mis, N_mis] cov_x2;
  matrix[N_mis, N_mis] cov_f2;
  matrix[N_mis, N_mis] diag_delta;
  matrix[N_obs, N_obs] K;
  matrix[N_mis, N_mis] sub_cov;
  int check1[2];
  int check2[2];

  var_vec_obs[1] = 0;
  var_vec_obs[2:m] = square(tau_glob) * diff .* square(lambda);
  for (k in 1:N_mis){
    if (use_mean){
      var_vec_mis[k] = square(tau_glob) * miss_diff[k] * (square(lambda[nearest[k]])*(miss_diff[k]/diff[nearest[k]]) + square(lambda[nearest[k]+1])*(diff[nearest[k]]-miss_diff[k])/diff[nearest[k]]);
    } else{
      var_vec_mis[k] = square(tau_glob) * miss_diff[k] * square(lambda[nearest[k]]);
    }
    var_mis[k] = sum(var_vec_obs[1:nearest[k]]) + var_vec_mis[k];
  }

  for (i in 1:N_obs){
        for (j in 1:N_obs){
          K[i,j] = sum(var_vec_obs[1:min(ind[i], ind[j])]);
        }
  }

  for (n in 1:N_obs) {
        K[n, n] = K[n, n] + square(sig);
  }
  L_K = cholesky_decompose(K);
  K_div_y1 = mdivide_left_tri_low(L_K, to_vector(y_obs));
  K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';

      for (q in 1:N_obs){
        for (r in 1:N_mis){
          if (ind[q] < nearest[r]){
              k_x1_x2[q,r] = sum(var_vec_obs[1:ind[q]]);
          } else{
              k_x1_x2[q,r] = var_mis[r];
          }
        }
      }
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);

      for (l in 1:N_mis){
        for (o in 1:N_mis){
          cov_x2[l,o] = var_mis[min(l,o)];
        }
      }

      sub_cov = v_pred' * v_pred;
}
