#' Function to fit a Bayesian horseshoe process regression
#'
#' This function fits a Bayesian horseshoe process regression: it calculates the
#' posterior distribution for a set of q horseshoe process terms on some continuous
#' predictors X, and a set of p linear predictors Z. Outcomes can be Gaussian,
#' binary, or Poisson distributed. Monotonic constraints can be enforced on the
#' horseshoe process terms, if desired. At present, this model is implemented for
#' q = 1 or 2, although future versions may permit higher values for q.
#' If the predictor(s) X are sparsely distributed over their domain, it may be
#' wise to add some additional augmented data to even the grid of X or to obtain
#' predictions/interpolations where desired. More information can be found in Chase et al. (2023).
#'
#' @param y
#' @param X
#' @param Z
#' @param X_aug
#' @param family
#' @param monotonic_terms
#' @param monotonic_approach
#' @param aug_approach
#' @param seed
#' @param chains
#' @param iter_warmup
#' @param iter_sampling
#' @param max_treedepth
#' @param adapt_delta
#' @param c
#' @param verbose
#'
#' @return
#'
#' @examples
#'
#' @importFrom dplyr near case_when
#' @importFrom stats quantile sd plogis qlogis
#' @import cmdstanr
#' @export
hpr <- function(y = NULL,
                X = NULL,
                Z = NULL,
                X_aug = NULL,
                family = NULL,
                monotonic_terms = NULL,
                monotonic_approach = NULL,
                aug_approach = "MCMC",
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

  N_obs <- nrow(X)
  y_obs <- y
  #y <- as.numeric(scale(mydata$Visits_twovar, scale = FALSE))

  if (!is.null(Z)){
    p <- ncol(Z)
    covariates <- as.matrix(scale(Z, scale = FALSE))
    has_covs <- 1
  } else{
    p <- 1
    covariates <- as.matrix(rep(0, N_obs), nrow = N_obs, ncol = 1)
    has_covs <- 0
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
        exp_approach = exp_approach
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
                 nonmonotone_inds = nonmonotone_inds
    )
  }

   if (family=="gaussian"){
     dat$samp_mean <- mean(y_obs)
     dat$samp_sd <- sd(y_obs)
     if (ncol(X) > 1){
       mymodel <- cmdstan_model(system.file("stan", "gaussian_multi.stan", package = "HPR"))
       model_file <- "gaussian_multi.stan"
     } else{
       if (aug_approach=="MCMC"){
         mymodel <- cmdstan_model(system.file("stan", "gaussian_uni.stan", package = "HPR"))
         model_file <- "gaussian_uni.stan"
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
