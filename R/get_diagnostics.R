#' Function to extract Stan diagnostics from an HPR model
#'
#' @param object The results object from a run of hpr.
#' @param verbose A logical indicator of whether a full cmdstan diagnostic report should
#' be printed to the console. The default is false.
#'
#' @return A dataframe with columns:
#' \describe{
#'    \item{Divergences}{The number of HMC samples that ended in a divergence.}
#'    \item{Max_Treedepth}{The number of HMC samples that had a max_treedepth warning.}
#'    \item{Rhat}{The number of f parameters that had Rhat greater than 1.1, using the adjusted Rhat
#'    of Vehtari et al. (Bayesian Analysis, 2021).}
#'    \item{Min_Ess_Bulk}{The minimum effective sample size in the bulk of the posterior across
#'    the f parameters. This estimated according to Vehtari et al. (Bayesian Analysis, 2021).}
#'    \item{Min_Ess_Tail}{The minimum effective sample size in the tails of the posterior across
#'    the f parameters. This estimated according to Vehtari et al. (Bayesian Analysis, 2021).}
#'    \item{Num_Param}{The length of the f vector (the systematic component of the model), which
#'    is a function of all other parameters in the model.}
#'    \item{Num_Samples}{The number of HMC posterior samples.}
#'    \item{Time}{The computing time of Stan sampling.}
#' }
#' For more information on these metrics, please see Chase et al. (2022+) or the Stan reference manual.
#'
#' @examples
#' X <- as.matrix(dat$Day, ncol = 1)
#' y <- dat$Temperature
#'
#' mymodel <- hpr(y = y, X = X, family = "gaussian")
#' get_diagnostics(mymodel)
#'
#' @importFrom posterior as_draws_df summarise_draws
#' @import cmdstanr
#' @export
get_diagnostics <- function(object = NULL,
                  verbose = FALSE
){
  if (verbose){
    object$stan_object$cmdstan_diagnose()
  }

  diagnostics_df <- as_draws_df(object$stan_object$sampler_diagnostics())
  f_df <- as_draws_df(object$stan_object$draws("f"))
  sum_stats <- summarise_draws(f_df)

  diagnostics_table <- data.frame("Divergences" = sum(diagnostics_df$divergent__),
                                  "Max_Treedepth" = sum(diagnostics_df$treedepth__ >= object$treedepth),
                                  "RHat" = length(which(sum_stats$rhat > 1.1)),
                                  "Min_Ess_Bulk" = min(sum_stats$ess_bulk),
                                  "Min_Ess_Tail" = min(sum_stats$ess_tail),
                                  "Num_Param" = nrow(sum_stats),
                                  "Num_Samples" = nrow(diagnostics_df),
                                  "Time" = object$run_time
                                  )

  return(diagnostics_table)
}
