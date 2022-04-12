#' Function to extract posterior predictions from an HPR model. Note that these are
#' posterior predictions; if you want estimates of the estimated trajectory f
#' (the systematic model component), use get_preds.
#'
#' @param object The results object from a run of hpr, for which a new_X and
#' new_Z matrix were input.
#' @param alpha The uncertainty level for posterior prediction intervals; the
#' default is 0.05 (which corresponds to 95\% prediction intervals).
#'
#' @return A dataframe with N_new rows (corresponding to the values of new_X, new_Z) and columns:
#' \describe{
#'    \item{Mean}{the mean of the posterior prediction samples}
#'    \item{Median}{the median of the posterior prediction samples}
#'    \item{Lower}{the alpha/2 percentile of the posterior prediction samples}
#'    \item{Upper}{the 1-alpha/2 percentile of the posterior prediction samples}
#' }
#'
#' @examples
#' X <- as.matrix(dat$Day, ncol = 1)
#' y <- dat$Temperature
#'
#' mymodel <- hpr(y = y, X = X, family = "gaussian", new_X = X)
#' post_preds <- post_predict(mymodel)
#' post_preds$x <- dat$Day
#'
#' @importFrom dplyr select
#' @importFrom stats plogis quantile
#' @importFrom posterior as_draws_df
#' @import cmdstanr
#' @export
post_predict <- function(object = NULL,
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
