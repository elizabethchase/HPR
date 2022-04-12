#' Function to fit a Bayesian horseshoe process regression
#'
#' This function fits a Bayesian horseshoe process regression: it calculates the
#' posterior distribution for a set of q horseshoe process terms on some continuous
#' predictors X, and a set of p linear predictors Z. Outcomes can be Gaussian,
#' binary, or Poisson distributed. Monotonic constraints can be enforced on the
#' horseshoe process terms, if desired. At present, this model is implemented for
#' q = 1 or 2, although future versions may permit higher values for q.
#' If the predictor(s) X are sparsely distributed over their domain, it may be
#' wise to add some additional augmented data to even out the grid of X or to obtain
#' predictions/interpolations where desired. More information can be found in Chase et al. (2022+).
#'
#' @param y A length N numeric vector containing the outcome of interest--if continuous, then this
#' should be real values; if binary, then a vector of 0's and 1's; if count data, then a vector of integers.
#' @param X An N x q matrix of the nonlinear predictors, with q either 1 or 2. All entries of X must be
#' numeric real.
#' @param Z An N x p matrix of linear predictors. All entries of X should be numeric; therefore, any categorical
#' predictors should be converted into dummy variables before calling hpr.
#' @param X_aug An N_aug x q matrix of nonlinear predictors at which interpolations/extrapolations are requested.
#' @param family The family of the outcome; can be either "gaussian", "binomial", or "poisson". The default is "gaussian".
#' @param beta_dist If a matrix Z is provided, this is the prior imposed on the coefficients of the linear predictors included
#' in Z. Can be either "normal" or "cauchy". The default is "normal".
#' @param beta_mean If a matrix Z is provided, then this is a vector of prior location parameters for the coefficients of the
#' linear predictors. The default behavior will set all prior locations to be 0.
#' @param beta_sd If a matrix Z is provided, then this is a vector of prior scale parameters for the coefficients of the
#' linear predictors. The default behavior will set all prior scales to be 5 times the observed standard deviation for that predictor.
#' @param beta_df If a matrix Z is provided and the beta distribution is "cauchy", this is a vector containing the degrees of freedom for each Cauchy
#' prior for each of the coefficients of the linear predictors. The default behavior is to set beta_df to be 1 for all coefficients,
#' which corresponds to a true Cauchy prior. Larger values of beta_df will yield t distributions with that many degrees of freedom.
#' @param monotonic_terms A vector of less than length q containing the indices of the nonlinear predictors that should be
#' constrained to be monotonic. See the vignette for more examples.
#' @param monotonic_approach An indicator of which approach for constraining monotonicity should be used. The default is "abs" (the
#' absolute value) which is recommended; users can also specify "exp" (the exponential function) which may have computational challenges.
#' @param aug_approach An indicator of which approach for data augmentation should be used. The highly recommended default is
#' "MCMC" which draws from the posterior to generate new augmentation values and local shrinkage parameters. Other options are
#' "Mean" or "LVCF", which use the mean value of the two neighboring local shrinkage parameters, or the last value of the local
#' shrinkage parameter to serve as the local shrinkage at the augmentation point, then inputs this value via Gaussian kriging equations.
#' For more details, see Chase et al. (2022+). Again, we highly discourage the use of any approach other than "MCMC".
#' @param new_X An N_new x q matrix of locations at which posterior predictions are requested. Note that all entries in this
#' matrix must already appear in either X or X_aug in order for posterior predictions to be generated.
#' @param new_Z An N_new x p matrix of linear predictor values at which posterior predictions are requested, corresponding to
#' the nonlinear predictor values input in new_X.
#' @param seed A number that will be passed to the Stan MCMC to make it possible to generate identical results on future runs.
#' @param chains The number of chains used for Hamiltonian MCMC; the default is 4. We do not recommend modifying this parameter
#' unless you know what you're doing.
#' @param iter_warmup The number of samples for warmup per chain. These samples will not be included for posterior estimation.
#' The default is 1000 (so 1000 x 4 chains = 4000 warmup samples). We do not recommend modifying this parameter unless you know
#' what you're doing.
#' @param iter_sampling The number of posterior samples per chain. The default is 2000 (so 2000 x 4 chains = 8000 posterior samples).
#' We do not recommend modifying this parameter unless you know what you're doing.
#' @param max_treedepth The max_treedepth parameter passed to the NUTS algorithm within the Hamiltonian MCMC. The default value is
#' 12. We do not recommend modifying this parameter unless you know what you're doing.
#' @param adapt_delta The adapt_delta parameter passed to the Hamiltonian MCMC, which helps to set the step size of the MCMC. The default
#' value is 0.95. We do not recommend modifying this parameter unless you know what you're doing.
#' @param c The scale parameter of the half-Cauchy prior placed on the global shrinkage parameter. The default is 0.01, which we have
#' found yields reasonably good performance. In general, we find that this parameter does not make much difference except in the
#' case of very sparse data, in which case changes in c can make a big difference. We recommend trying several values of c to
#' explore sensitivity of results.
#' @param verbose A logical value indicating whether or not you would like to print updates on Stan sampling progress. The default is
#' verbose = FALSE.
#'
#' @return The function returns a named list with the following contents:
#' \describe{
#'    \item{stan_object}{the finished Stan modeling object}
#'    \item{run_time}{the length of time needed for HMC sampling, in seconds}
#'    \item{model_file}{a character string containing the name of the Stan model file that was used}
#'    \item{seed}{either the value of seed input by user, or the value of a seed that was randomly generated within hpr}
#'    \item{grid}{the unique, ordered values of each column of X and X_aug}
#'    \item{treedepth}{the value of max_treedepth as input by the user}
#'    \item{adapt_delta}{the value of adapt_delta as input by the user}
#' }
#'
#' @examples
#'X <- as.matrix(dat$Day, ncol = 1)
#'y <- dat$Temperature
#'
#'mymodel <- hpr(y = y, X = X, family = "gaussian")
#'
#'mymodel_monotonic <- hpr(y = y, X = X, family = "gaussian",
#'monotonic_terms = c(1), monotonic_approach = "abs")
#'
#'Z <- as.matrix(dat$Aberration, ncol = 1)
#'
#'mymodel2 <- hpr(y = y, X = X, Z = Z, new_X = X, new_Z = Z,
#'family = "gaussian", beta_dist = "cauchy")
#'
#' @importFrom dplyr near case_when
#' @importFrom stats quantile sd plogis qlogis
#' @import cmdstanr
#' @export
hpr <- function(y = NULL,
                X = NULL,
                Z = NULL,
                X_aug = NULL,
                family = "gaussian",
                beta_dist = "normal",
                beta_mean = NULL,
                beta_sd = NULL,
                beta_df = NULL,
                monotonic_terms = NULL,
                monotonic_approach = NULL,
                aug_approach = "MCMC",
                new_X = NULL,
                new_Z = NULL,
                seed = NULL,
                chains = 4,
                iter_warmup = 1000,
                iter_sampling = 2000,
                max_treedepth = 12,
                adapt_delta = 0.95,
                c = 0.01,
                verbose = FALSE
                ){

  if (!is.matrix(X)){stop("X must be a matrix.")}
  if (!is.vector(y)){stop("y must be a vector.")}
  if (is.character(y) | is.logical(y)){stop("y must be numeric or integer.")}
  if (is.null(family)){stop("Please specify a family of model: gaussian, binomial, or poisson.")}
  if (!near(length(y), nrow(X))){stop("y must be the same length as the number of rows in X!")}
  if (family=="gaussian" & !is.numeric(y)){"gaussian data must be of type numeric!"}

  if (!is.null(Z)){
    if (!is.matrix(Z)){stop("Z must be a matrix.")}
    if ((is.character(Z) | is.logical(Z))){stop("Z must be numeric or integer.")}
    if (!near(length(y), nrow(Z))){stop("y must be the same length as the number of rows in Z!")}
    if (!near(nrow(Z), nrow(X))){stop("X and Z must have the same number of rows!")}
    if (!is.null(beta_mean) & !near(length(beta_mean), ncol(Z))){stop("beta_mean must have the
                                                                      same number of entries as the
                                                                      number of columns in Z!")}
    if (!is.null(beta_sd) & !near(length(beta_sd), ncol(Z))){stop("beta_sd must have the
                                                                      same number of entries as the
                                                                      number of columns in Z!")}
    if (!is.null(beta_df) & !near(length(beta_df), ncol(Z))){stop("beta_df must have the
                                                                      same number of entries as the
                                                                      number of columns in Z!")}
    if (is.null(beta_mean)){
      beta_mean <- rep(0.0, ncol(Z))
    }
    if (is.null(beta_sd)){
      beta_sd <- apply(Z, 2, sd)
    }
    if (is.null(beta_df)){
      beta_df <- rep(1.0, ncol(Z))
    }
  }

  if (!is.null(X_aug)){
    if (!near(length(X_aug),ncol(X))){stop("X_aug must be a list with the same number
                                                            of entries as the number of columns in X.")}
  }

  if (!is.null(monotonic_terms)){
    if((min(monotonic_terms) < 0 | max(monotonic_terms) > ncol(X)))
    {stop("monotonic_terms must be a vector containing the indices of the columns of X that are monotonic.")}
    if (is.null(monotonic_approach)){"Please specify an approach for constraining
    monotonicity: either exp or abs."}
  }

  if (aug_approach != "MCMC" & (family != "gaussian" | ncol(X) > 1 | !is.null(monotonic_terms)))
  {"Augmentation approaches other than MCMC are not recommended, and are only implemented for Gaussian outcomes and
    a single, non-monotonic horseshoe smoothing term."}

  if(aug_approach != "MCMC" & !is.null(new_X)){stop("Additional predictions cannot be generated for models that rely on interpolated augmentation!")}
  if(aug_approach != "MCMC" & !is.null(new_Z)){stop("Additional predictions cannot be generated for models that rely on interpolated augmentation!")}

  if(is.null(Z) & !is.null(new_Z)){warning("No Z was inputted, so new_Z will not be used.")}
  if(!is.null(Z) & !is.null(new_X) & is.null(new_Z)){stop("This model has linear predictors--please input a matrix new_Z.")}

  if (!is.null(new_X)){
    if(!near(ncol(new_X), ncol(X))){stop(paste0("new_X has ", ncol(new_X), " columns, but orig_X has ", ncol(orig_X),
                                                " columns. They should have the same number of columns."))}

    for (i in 1:ncol(X)){
      if(!near(length(setdiff(new_X[i], X[i])),0)){stop(paste0("All elements of new_X must already appear in the
                                                                 corresponding column of X. Column ", i, " of
                                                                 new_X has an element that is not in X. To get
                                                                 predictions at unobserved locations, use the X_aug
                                                                 argument of hpr."))}
    }
  }

  if (!is.null(new_Z)){
    if(!near(nrow(new_X), nrow(new_Z))){stop("new_X and new_Z must have the same number of rows!")}
    if(!near(ncol(Z), ncol(new_Z))){stop("new_Z should have ", ncol(Z), " columns, but it has ", ncol(new_Z), " columns.")}
  }

  N_obs <- nrow(X)
  y_obs <- y
  #y <- as.numeric(scale(mydata$Visits_twovar, scale = FALSE))

  if (!is.null(Z)){
    p <- ncol(Z)
    covariates <- as.matrix(scale(Z, scale = FALSE))
    scale_means <- colMeans(Z)
    has_covs <- 1
    if(beta_dist=="cauchy"){
      cauchy_beta <- 1
    } else{
      cauchy_beta <- 0
    }
    if (is.null(new_Z)){
      N_new <- 1
      new_covariates <- as.matrix(covariates[1,], nrow = 1)
    } else{
      N_new <- nrow(new_Z)
      scale_new_Z <- matrix(NA, nrow = nrow(new_Z), ncol = ncol(new_Z))
      for (i in 1:ncol(Z)){
        scale_new_Z[,i] <- new_Z[,i] - scale_means[i]
      }
      new_covariates <- as.matrix(scale_new_Z)
    }
  } else{
    p <- 1
    covariates <- as.matrix(rep(0, N_obs), nrow = N_obs, ncol = 1)
    has_covs <- 0
    cauchy_beta <- 0
    beta_mean <- 0
    beta_sd <- 1
    beta_df <- 1
    if (is.null(new_X)){
      N_new <- 1
      new_covariates <- as.matrix(rep(0, N_new), nrow = N_new, ncol = 1)
    } else{
      N_new <- nrow(new_X)
      new_covariates <- as.matrix(rep(0, N_new), nrow = N_new, ncol = 1)
    }
  }

  if (ncol(X) < 2){
    x_obs <- X[,1]
    if (is.null(X_aug) | aug_approach != "MCMC"){
      x_aug <- NULL
      N_mis <- 0
      missing_ind <- rep(0,0)
    } else{
      x_aug <- X_aug[[1]]
    }
    xdat <- c(x_obs, x_aug)
    xdat_sorted <- sort(xdat)
    grid <- xdat_sorted[c(TRUE, !near(diff(xdat_sorted), 0))]
    m <- length(grid)
    differences <- diff(grid)
    ind <- rep(NA, nrow = N_obs)
    for (i in 1:N_obs){
      ind[i] <- max(which((x_obs[i] >= grid) | near(x_obs[i], grid)))
    }
    if (!is.null(x_aug) & aug_approach=="MCMC"){
     missing_ind <- setdiff(c(1:m), ind)
     N_mis <- length(missing_ind)
    }
    if (!is.null(new_X) & aug_approach=="MCMC"){
      new_ind <- rep(NA, nrow = N_new)
      for (i in 1:N_new){
        new_ind[i] <- max(which((new_X[i,1] >= grid) | near(new_X[i,1], grid)))
      }
    } else{
      new_ind <- 1
    }
    if (aug_approach != "MCMC" & family=="gaussian"){
      x_aug <- X_aug[[1]]
      x_aug2 <- setdiff(x_aug, x_obs)
      if (min(x_aug2) < min(grid)){stop("If using LVCF or Mean imputation, the smallest
                                        augmented value must be greater than or
                                        equal to the smallest observed value.")}
      N_mis <- length(x_aug2)
      nearest <- rep(NA, N_mis)
      miss_diff <- rep(NA, N_mis)
      for (i in 1:N_mis){
        myind <- max(which(grid <= x_aug2[i]))
        nearest[i] <- myind
        miss_diff[i] <- x_aug2[i] - grid[myind]
      }
    }
    if (!is.null(monotonic_terms)){
      is_monotone <- 1
      if (monotonic_approach=="exp"){exp_approach <- 1}
      else {exp_approach <- 0}}
    else {is_monotone <- 0
          exp_approach <- 0
          }

    if (aug_approach=="MCMC"){
      dat <- list(
        N_obs = N_obs,
        N_mis = N_mis,
        y_obs = y_obs,
        m = m,
        ind = ind,
        missing_ind = missing_ind,
        diff = differences,
        p = p,
        X = covariates,
        has_covs = has_covs,
        local_dof_stan = 1,
        global_dof_stan = 1,
        alpha_scale_stan = c,
        is_monotone = is_monotone,
        exp_approach = exp_approach,
        N_new = N_new,
        new_ind = new_ind,
        new_X = new_covariates
      )
    } else if (aug_approach!="MCMC"){
      dat <- list(
        N_obs = N_obs,
        N_mis = N_mis,
        y_obs = y_obs,
        m = m,
        ind = ind,
        nearest = nearest,
        miss_diff = miss_diff,
        diff = differences,
        p = p,
        X = covariates,
        has_covs = has_covs,
        local_dof_stan = 1,
        global_dof_stan = 1,
        alpha_scale_stan = c,
        x_aug = x_aug2
      )
    }
  } else {
    k <- ncol(X)
    m <- rep(NA, k)
    ind <- matrix(NA, nrow = N_obs, ncol = k)
    new_ind <- matrix(NA, nrow = N_new, ncol = k)
    grid <- vector("list", length = k)
    for (i in 1:k){
      x_aug <- X_aug[[i]]
      x_obs <- X[,i]
      xdat <- c(x_obs, x_aug)
      xdat_sorted <- sort(xdat)
      grid[[i]] <- xdat_sorted[c(TRUE, !near(diff(xdat_sorted), 0))]
      m[i] <- length(grid[[i]])
      for (j in 1:length(x_obs)){
        ind[j,i] <- max(which((x_obs[j] >= grid[[i]]) | near(x_obs[j], grid[[i]])))
      }
      if (!is.null(new_X)){
        for (j in 1:N_new){
          new_ind[j,i] <- max(which((new_X[j, i] >= grid[[i]]) | near(new_X[j, i], grid[[i]])))
        }
      } else {
        new_ind <- matrix(rep(1, k), nrow = 1, ncol = k)
      }
    }
    max_m <- max(m)
    differences <- matrix(NA, nrow = max_m-1, ncol = k)
    for (i in 1:k){
      differences[,i] <- c(diff(grid[[i]]), rep(0, max_m - m[i]))
    }

    if (is.null(X_aug)){
      N_mis <- 0
      missing_ind <- matrix(rep(0,0), nrow = 0, ncol = k)
    } else{
      lengths <- rep(NA, k)
      missing_grids <- vector("list", length = k)
      for (i in 1:k){
        missing_grids[[i]] <- setdiff(1:m[i], ind[,i])
        lengths[i] <- length(missing_grids[[i]])
      }
      lengths <- c(0, lengths)
      missing_ind <- matrix(NA, nrow = sum(lengths), ncol = k)
      for (i in 1:k){
        missing_ind[(lengths[i]+1):(lengths[i]+lengths[i+1]),i] <- missing_grids[[i]]
      }
      for (i in 1:k){
        missing_ind[is.na(missing_ind[,i]),i] <- sample(1:m[i], length(which(is.na(missing_ind[,i]))))
      }
    }

    if (!is.null(monotonic_terms)){
      is_monotone <- 1
      if (monotonic_approach=="exp"){exp_approach <- 1}
      else {exp_approach <- 0}
      q <- length(monotonic_terms)
      monotone_inds <- monotonic_terms
      nonmonotone_inds <- setdiff(1:k, monotone_inds)
      }
    else {
      is_monotone <- 0
      q <- 0
      monotone_inds <- rep(0,0)
      nonmonotone_inds <- c(1:k)
      exp_approach <- 0
      }

    dat <- list(N_obs = N_obs,
                 N_mis = N_mis,
                 y_obs = y_obs,
                 k = k,
                 m = m,
                 max_m = max_m,
                 ind = ind,
                 missing_ind = missing_ind,
                 diff = differences,
                 p = p,
                 X = covariates,
                 has_covs = has_covs,
                 local_dof_stan = 1,
                 global_dof_stan = 1,
                 alpha_scale_stan = c,
                 is_monotone = is_monotone,
                 exp_approach = exp_approach,
                 q = q,
                 monotone_inds = monotone_inds,
                 nonmonotone_inds = nonmonotone_inds,
                 N_new = N_new,
                 new_ind = new_ind,
                 new_X = new_covariates
    )
  }

   if (family=="gaussian"){
     dat$samp_mean <- mean(y_obs)
     dat$samp_sd <- sd(y_obs)
     if (ncol(X) > 1){
       mymodel <- cmdstan_model(system.file("stan", "gaussian_multi.stan", package = "HPR"))
       model_file <- "gaussian_multi.stan"
       dat$beta_mean <- beta_mean
       dat$beta_sd <- beta_sd
       dat$beta_df <- beta_df
       dat$cauchy_beta <- cauchy_beta
     } else{
       if (aug_approach=="MCMC"){
         mymodel <- cmdstan_model(system.file("stan", "gaussian_uni.stan", package = "HPR"))
         model_file <- "gaussian_uni.stan"
         dat$beta_mean <- beta_mean
         dat$beta_sd <- beta_sd
         dat$beta_df <- beta_df
         dat$cauchy_beta <- cauchy_beta
       } else if (aug_approach=="LVCF"){
         mymodel <- cmdstan_model(system.file("stan", "gaussian_uni_interp.stan", package = "HPR"))
         model_file <- "gaussian_uni_interp.stan"
         dat$use_mean <- 0
       } else if (aug_approach=="Mean"){
         mymodel <- cmdstan_model(system.file("stan", "gaussian_uni_interp.stan", package = "HPR"))
         model_file <- "gaussian_uni_interp.stan"
         dat$use_mean <- 1
       }
     }
   } else if (family == "binomial"){
     trunc_y <- case_when(
       y_obs==1 ~ 0.995,
       y_obs==0 ~ 0.005
     )

     dat$samp_mean <- qlogis(mean(trunc_y))
     dat$samp_sd = sd(qlogis(trunc_y))
     dat$beta_mean <- beta_mean
     dat$beta_sd <- beta_sd
     dat$beta_df <- beta_df
     dat$cauchy_beta <- cauchy_beta

     if (ncol(X) > 1){
       mymodel <- cmdstan_model(system.file("stan", "binomial_multi.stan", package = "HPR"))
       model_file <- "binomial_multi.stan"
     } else{
       mymodel <- cmdstan_model(system.file("stan", "binomial_uni.stan", package = "HPR"))
       model_file <- "binomial_uni.stan"
     }
   } else if (family=="poisson"){
     dat$samp_mean <- log(mean(y_obs))
     dat$samp_sd = sd(log(y_obs+0.5))
     dat$beta_mean <- beta_mean
     dat$beta_sd <- beta_sd
     dat$beta_df <- beta_df
     dat$cauchy_beta <- cauchy_beta

     if (ncol(X) > 1){
       mymodel <- cmdstan_model(system.file("stan", "poisson_multi.stan", package = "HPR"))
       model_file <- "poisson_multi.stan"
     } else{
       mymodel <- cmdstan_model(system.file("stan", "poisson_uni.stan", package = "HPR"))
       model_file <- "poisson_uni.stan"
     }
   }

  if (verbose){
    refresh <- 200
  } else{
    refresh <- 0
  }

  if (is.null(seed)){
    seed <- sample(1:10000000, size = 1)
  }

  start.time <- Sys.time()
  my_finished_model <- mymodel$sample(
    data = dat,
    chains = chains,
    seed = seed,
    refresh = refresh,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  end.time <- Sys.time()
  time.taken <- as.numeric(difftime(end.time, start.time, units="secs"))

  outputs <- list(
    "stan_object" = my_finished_model,
    "run_time" = time.taken,
    "data" = dat,
    "model_file" = model_file,
    "seed" = seed,
    "grid" = grid,
    "treedepth" = max_treedepth,
    "adapt_delta" = adapt_delta
  )

  return(outputs)
}
