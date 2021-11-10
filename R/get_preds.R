#' Function to extract predicted y values from an HPR model
#'
#' @param object
#' @param alpha
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
get_preds <- function(object = NULL,
                  alpha = 0.05
){
  family <- strsplit(object$model_file, "[_]")[[1]][1]
  lower <- alpha/2
  upper <- 1-(alpha/2)

  my_summary <- dplyr::select(as_draws_df(object$stan_object$draws("f")), -.chain, -.iteration, -.draw)

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
