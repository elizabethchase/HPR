#' Function to extract Stan diagnostics from an HPR model
#'
#' @param object
#'
#' @return
#'
#' @examples
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
