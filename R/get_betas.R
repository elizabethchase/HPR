#' Function to extract the estimates of the linear coefficients from an HPR model
#'
#' @param object The results object from a run of hpr for which a matrix Z was input.
#' @param alpha The uncertainty level for the betas; the
#' default is 0.05 (which corresponds to 95\% credible intervals).
#' @param var_names An optional character vector containing the variable names of the
#' columns of Z.
#'
#' @return A dataframe with a row for each beta and columns:
#' \describe{
#'    \item{Variable}{the name of the variable, if var_names was used}
#'    \item{Mean}{the mean of the posterior samples of the beta}
#'    \item{Median}{the median of the posterior samples of the beta}
#'    \item{Lower}{the alpha/2 percentile of the posterior samples of the beta}
#'    \item{Upper}{the 1-alpha/2 percentile of the posterior samples of the beta}
#' }
#'
#' @examples
#' X <- as.matrix(dat$Day, ncol = 1)
#' y <- dat$Temperature
#' Z <- as.matrix(dat$Aberration, ncol = 1)
#'
#' mymodel <- hpr(y = y, X = X, Z = Z, family = "gaussian", beta_dist = "cauchy")
#' mybetas <- get_betas(mymodel, var_names = c("Aberration"))
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
  scales <- object$scale_sd
  if (has_covs==0){stop("This model object has no linear coefficients.")}
  lower <- alpha/2
  upper <- 1-(alpha/2)
  if (is.null(var_names)){var_names <- 1:p}

  my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("beta")), -.iteration, -.chain, -.draw)
  for(i in 1:p){
    my_summary[,i] <- my_summary[,i]*scales[i]
  }

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
