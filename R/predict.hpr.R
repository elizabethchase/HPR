#' Function to extract the predicted horseshoe curve fits from an HPR model
#'
#' @param object
#' @param orig_X
#' @param new_X
#' @param new_Z
#' @param alpha
#'
#' @return
#'
#' @examples
#'
#' @importFrom dplyr near select
#' @importFrom stats plogis quantile
#' @importFrom posterior as_draws_df
#' @importFrom mgcv rmvn
#' @import cmdstanr
#' @export
predict.hpr <- function(object = NULL,
                  orig_X = NULL,
                  new_X = NULL,
                  new_Z = NULL,
                  alpha = 0.05
){
  family <- strsplit(object$model_file, "[_]")[[1]][1]
  num_var <- strsplit(strsplit(object$model_file, "[_]")[[1]][2], "[.]")[[1]][1]
  posthoc_interp <- grepl("interp", object$model_file, fixed = TRUE)
  lower <- alpha/2
  upper <- 1-(alpha/2)
  if (posthoc_interp){stop("Predictions are not generated for posthoc interpolation.")}

  new_model_file <- paste0(family, "_", num_var, "_gq.stan")
  mymodel <- cmdstan_model(system.file("stan", new_model_file, package = "HPR"))

  new_preds <- mymodel$generate_quantities(object$stan_object, data = object$data, seed = object$seed)

  my_summary <- dplyr::select(as_draws_df(new_preds$draws("new_pred")), -.chain, -.iteration, -.draw)

  if (family=="gaussian"){
    preds <- data.frame("Median" = apply(my_summary, MARGIN = 2, FUN = quantile, 0.5),
                        "Mean" = apply(my_summary, MARGIN = 2, FUN = mean),
                        "Lower" = apply(my_summary, MARGIN = 2, FUN = quantile, lower),
                        "Upper" = apply(my_summary, MARGIN = 2, FUN = quantile, upper))
  } else if (family=="poisson"){
    preds <- data.frame("Median" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5)),
                        "Mean" = exp(apply(my_summary, MARGIN = 2, FUN = mean)),
                        "Lower" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, lower)),
                        "Upper" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, upper)))
  } else if (family=="binomial"){
    preds <- data.frame("Median" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5)),
                        "Mean" = plogis(apply(my_summary, MARGIN = 2, FUN = mean)),
                        "Lower" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, lower)),
                        "Upper" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, upper)))
  }

  return(preds)

}
