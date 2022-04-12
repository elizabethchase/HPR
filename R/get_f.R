#' Function to extract the nonlinear component of a horseshoe model
#'
#' Function to extract only the nonlinear component of a predicted horseshoe curve fit
#' from an HPR model at each unique, ordered value of X. If no Z matrix was used,
#' get_f will return the same results as get_preds. If a Z matrix was used,
#' get_f will return the nonlinear component estimated at the average value of each
#' predictor in Z.
#'
#' @param object The results object from a run of hpr.
#' @param alpha The uncertainty level for the nonlinear component of f; the
#' default is 0.05 (which corresponds to 95\% credible intervals).
#'
#' @return A dataframe with m x q rows (corresponding to the unique, ordered values of each column of X)
#' and columns:
#' \describe{
#'    \item{x}{The value of X}
#'    \item{Median}{the median of the posterior samples of the nonlinear component}
#'    \item{Mean}{the mean of the posterior samples of the nonlinear component}
#'    \item{Lower}{the alpha/2 percentile of the posterior samples of the nonlinear component}
#'    \item{Upper}{the 1-alpha/2 percentile of the posterior samples of the nonlinear component}
#'    \item{Predictor}{if q > 1, then this gives the column index of X for which the nonlinear component
#'    is being estimated}
#' }
#'
#' @examples
#' X <- as.matrix(dat$Day, ncol = 1)
#' y <- dat$Temperature
#'
#' mymodel <- hpr(y = y, X = X, family = "gaussian")
#' my_nonlin <- get_f(mymodel)
#'
#' @importFrom dplyr near select
#' @importFrom stats plogis quantile
#' @importFrom posterior as_draws_df
#' @importFrom mgcv rmvn
#' @import cmdstanr
#' @export
get_f <- function(object = NULL,
                alpha = 0.05
){
  family <- strsplit(object$model_file, "[_]")[[1]][1]
  posthoc_interp <- grepl("interp", object$model_file, fixed = TRUE)
  m <- object$data$m
  k <- length(m)
  lower <- alpha/2
  upper <- 1-(alpha/2)

  if (posthoc_interp){
    N_mis <- object$data$N_mis
    my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("path")), -.iteration, -.chain, -.draw)
    cov_x2 <- as_draws_df(object$stan_object$draws(c("cov_x2")))
    sub_cov <- as_draws_df(object$stan_object$draws(c("sub_cov")))
    f2_mu <- as_draws_df(object$stan_object$draws(c("f2_mu")))
    delta_diag <- diag(x = 1e-9, nrow = N_mis)

    interps <- matrix(NA, nrow = nrow(cov_x2), ncol = N_mis)
    for (i in 1:nrow(cov_x2)){
      cov <- matrix(data = as.numeric(cov_x2[i, 1:N_mis^2]), nrow = N_mis) -
        matrix(data = as.numeric(sub_cov[i, 1:N_mis^2]), nrow = N_mis) + delta_diag
      interps[i,] <- rmvn(n = 1, mu = f2_mu[i, 1:N_mis], V = cov)
    }

    curve_fit <- data.frame("Median" = c(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5),
                                         apply(interps, MARGIN = 2, FUN = quantile, 0.5)),
                            "Mean" = c(apply(my_summary, MARGIN = 2, FUN = mean),
                                         apply(interps, MARGIN = 2, FUN = mean)),
                            "Lower" = c(apply(my_summary, MARGIN = 2, FUN = quantile, lower),
                                        apply(interps, MARGIN = 2, FUN = quantile, lower)),
                            "Upper" = c(apply(my_summary, MARGIN = 2, FUN = quantile, upper),
                                        apply(interps, MARGIN = 2, FUN = quantile, upper)),
                            "x" = c(object$grid, object$data$x_aug))

  } else{
    if (k<2){
      my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("path")), -.chain, -.iteration, -.draw)

      if (family=="gaussian"){
        curve_fit <- data.frame("Median" = apply(my_summary, MARGIN = 2, FUN = quantile, 0.5),
                              "Mean" = apply(my_summary, MARGIN = 2, FUN = mean),
                              "Lower" = apply(my_summary, MARGIN = 2, FUN = quantile, lower),
                              "Upper" = apply(my_summary, MARGIN = 2, FUN = quantile, upper),
                              "x" = object$grid)
      } else if (family=="poisson"){
        curve_fit <- data.frame("Median" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5)),
                              "Mean" = exp(apply(my_summary, MARGIN = 2, FUN = mean)),
                              "Lower" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, lower)),
                              "Upper" = exp(apply(my_summary, MARGIN = 2, FUN = quantile, upper)),
                              "x" = object$grid)
      } else if (family=="binomial"){
        curve_fit <- data.frame("Median" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, 0.5)),
                              "Mean" = plogis(apply(my_summary, MARGIN = 2, FUN = mean)),
                              "Lower" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, lower)),
                              "Upper" = plogis(apply(my_summary, MARGIN = 2, FUN = quantile, upper)),
                              "x" = object$grid)
      }
    } else {
      my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("finalpath")), -.chain, -.iteration, -.draw)
      m <- c(0,m)
      curve_fit <- data.frame("Median" = NA, "Mean" = NA, "Lower" = NA, "Upper" = NA, "x" = NA, "Predictor" = NA)
      if (family=="gaussian"){
        for (i in 1:k){
          subdat <- data.frame("Median" = apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, 0.5),
                               "Mean" = apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = mean),
                               "Lower" = apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, lower),
                               "Upper" = apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, upper),
                               "x" = object$grid[[i]],
                               "Predictor" = i
          )
          curve_fit <- rbind(curve_fit, subdat)
        }
      } else if (family=="binomial"){
        for (i in 1:k){
          subdat <- data.frame("Median" = plogis(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, 0.5)),
                               "Mean" = plogis(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = mean)),
                               "Lower" = plogis(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, lower)),
                               "Upper" = plogis(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, upper)),
                               "x" = object$grid[[i]],
                               "Predictor" = i
          )
          curve_fit <- rbind(curve_fit, subdat)
        }
      } else if (family=="poisson"){
        for (i in 1:k){
          subdat <- data.frame("Median" = exp(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, 0.5)),
                               "Mean" = exp(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = mean)),
                               "Lower" = exp(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, lower)),
                               "Upper" = exp(apply(my_summary[,(m[i]+1):(m[i] + m[i+1])], MARGIN = 2, FUN = quantile, upper)),
                               "x" = object$grid[[i]],
                               "Predictor" = i
          )
          curve_fit <- rbind(curve_fit, subdat)
        }
      }

      curve_fit <- curve_fit[-1,]

    }
  }

  return(curve_fit)
}
