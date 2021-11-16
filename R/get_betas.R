#' Function to extract the estimates of the linear coefficients from an HPR model
#'
#' @param object
#' @param alpha
#' @param var_names
#'
#' @return
#'
#' @examples
#'
#' @importFrom dplyr near select
#' @importFrom posterior as_draws_df
#' @importFrom stats plogis quantile
#' @import cmdstanr
#' @export
get_betas <- function(object = NULL,
                  alpha = 0.05,
                  var_names = NULL
){
  family <- strsplit(object$model_file, "[_]")[[1]][1]
  has_covs <- object$data$has_covs
  p <- object$data$p
  if (has_covs==0){stop("This model object has no linear coefficients.")}
  lower <- alpha/2
  upper <- 1-(alpha/2)
  if (is.null(var_names)){var_names <- 1:p}

  my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("beta")), -.iteration, -.chain, -.draw)

  if (family=="gaussian"){
    beta_preds <- data.frame("Variable" = var_names,
                             "Median" = apply(my_summary, MARGIN = 2, quantile, 0.5),
                             "Mean" = apply(my_summary, MARGIN = 2, mean),
                             "Lower" = apply(my_summary, MARGIN = 2, quantile, lower),
                             "Upper" = apply(my_summary, MARGIN = 2, quantile, upper))
  } else if (family=="binomial"){
    beta_preds <- data.frame("Variable" = var_names,
                             "Median" = exp(apply(my_summary, MARGIN = 2, quantile, 0.5)),
                             "Mean" = exp(apply(my_summary, MARGIN = 2, mean)),
                             "Lower" = exp(apply(my_summary, MARGIN = 2, quantile, lower)),
                             "Upper" = exp(apply(my_summary, MARGIN = 2, quantile, upper)))
  } else if (family=="poisson"){
    beta_preds <- data.frame("Variable" = var_names,
                             "Median" = exp(apply(my_summary, MARGIN = 2, quantile, 0.5)),
                             "Mean" = exp(apply(my_summary, MARGIN = 2, mean)),
                             "Lower" = exp(apply(my_summary, MARGIN = 2, quantile, lower)),
                             "Upper" = exp(apply(my_summary, MARGIN = 2, quantile, upper)))
  }

 return(beta_preds)

}
