#' Function to extract the nonlinear component of an HPR model, in the presence
#' of augmented data
#'
#' Function to extract the estimated value of the horseshoe component of the model
#' from an HPR model in the presence of augmented data. If there are no augmented data,
#' use get_preds. For posterior predictions, use post_predict.
#'
#' @param object The results object from a run of hpr.
#' @param alpha The uncertainty level for the values of the ho; the
#' default is 0.05 (which corresponds to 95\% credible intervals).
#'
#' @return A dataframe with N rows (corresponding to the rows of X, X_aug) and columns:
#' \describe{
#'    \item{Mean}{the mean of the posterior samples of the horseshoe component}
#'    \item{Median}{the median of the posterior samples of the horseshoe component}
#'    \item{Lower}{the alpha/2 percentile of the posterior samples of the horseshoe component}
#'    \item{Upper}{the 1-alpha/2 percentile of the posterior samples of the horseshoe component}
#' }
#'
#' @examples
#' X <- as.matrix(dat$Day, ncol = 1)
#' y <- dat$Temperature
#'
#' mymodel <- hpr(y = y, X = X, family = "gaussian")
#' my_H <- get_H(mymodel)
#'
#' @importFrom dplyr near select
#' @importFrom posterior as_draws_df
#' @importFrom stats plogis quantile
#' @import cmdstanr
#' @export
get_H <- function(object = NULL,
                      alpha = 0.05
){
  family <- strsplit(object$model_file, "[_]")[[1]][1]
  lower <- alpha/2
  upper <- 1-(alpha/2)

  my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("path")), -.chain, -.iteration, -.draw)

  if (family=="gaussian"){
    preds <- data.frame("Median" = apply(my_summary, MARGIN = 2, FUN = quantile, 0.5),
                        "Mean" = apply(my_summary, MARGIN = 2, FUN = mean),
                        "Lower" = apply(my_summary, MARGIN = 2, FUN = quantile, lower),
                        "Upper" = apply(my_summary, MARGIN = 2, FUN = quantile, upper),
                        "x" = object$grid)
  } else if (family=="poisson"){
    preds <- data.frame("Median" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5)),
                        "Mean" = exp(apply(my_summary, MARGIN = 2, FUN = mean)),
                        "Lower" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, lower)),
                        "Upper" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, upper)),
                        "x" = object$grid)
  } else if (family=="binomial"){
    preds <- data.frame("Median" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5)),
                        "Mean" = plogis(apply(my_summary, MARGIN = 2, FUN = mean)),
                        "Lower" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, lower)),
                        "Upper" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, upper)),
                        "x" = object$grid)
  }

  return(preds)
}
