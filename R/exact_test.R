#' Calculate the power of a one sample exact test for rate
#'
#' @description
#' This function calculate the power of a one sample exact test for rate.
#'
#' @param n A numeric vector of the total sample size, power will be calculated for
#' each n.
#' @param p The rate in alternative hypothesis.
#' @param null_p The rate in null hypothesis.
#' @param alpha The two-sided significance level. Default is 0.05, which is a 0.025
#' one-sided test.
#' @param side The side of the test, "U" and "L" for one-sided test in upper or
#' lower direction. "Both" for a two-sided equivalence test. Default is "U".
#' @param margin The effect to be tested. Default = 0.
#'
#' Recommend to adjust `null_p` directly and not change `margin`.
#'
#' The hypothesis is constructed as:
#'
#' For `side` = "U"
#'
#' H0: p - null_p <= margin
#'
#' H1: p - null_p > margin
#'
#' For `side` = "L"
#'
#' H0: p - null_p >= margin
#'
#' H1: p - null_p < margin
#'
#' For `side` = "Both"
#'
#' H0: p - null_p = margin
#'
#' H1: p - null_p != margin
#'
#' @returns A dataframe of sample size n and power of the test.
#'
#' @importFrom stats qbeta pf
#' @importFrom purrr pmap_dbl
#' @import dplyr
#' @export
#' @examples
#' # example code
#' power_exact_prop(671, p = 0.66, null_p = 0.60)
#'
#' # For multiple n
#' power_exact_prop(600:700, p = 0.66, null_p = 0.60)

power_exact_prop <- function(n, p, null_p, alpha = 0.05, side = c("U", "L", "Both"),
                             margin = 0){
  stopifnot(
    "all n must be integer larger than 1" = all(n == round(n)) && all(n >= 1),
    p >= 0,
    p <= 1,
    null_p >= 0,
    null_p <= 1,
    abs(null_p + margin) <= 1
  )

  side <- match.arg(side)

  if(side != "Both"){

    # The real null hypothesis if input a margin
    null_p <- null_p + margin
  }

  if (side == "U"){
    upper_crit <- purrr::map_dbl(.x = n,
                                 .f = function(x){
                                   above_null <- ceiling(x * null_p):x
                                   valid_cp_l <- qbeta(alpha / 2, above_null, x - above_null + 1) > null_p
                                   if(any(valid_cp_l)){
                                     above_null[which.max(valid_cp_l)]
                                   }else{
                                     NA
                                   }
                                 })
    power <- vector("double", length = length(n))

    v1 <- 2 * na.omit(upper_crit)
    v2 <- 2 * (n - na.omit(upper_crit) + 1)

    power[is.na(upper_crit)] <- 0
    power[!is.na(upper_crit)] <- pf(v2 * p / (v1 * (1 - p)), v1, v2)

  } else if(side == "L"){

    lower_crit <- purrr::map_dbl(.x = n,
                                 .f = function(x){
                                   below_null <- floor(x * null_p):0
                                   valid_cp_u <- qbeta(alpha / 2, below_null + 1, x - below_null, lower.tail = F) < null_p
                                   if(any(valid_cp_u)){
                                     below_null[which.max(valid_cp_u)]
                                   }else{
                                     NA
                                   }
                                 })
    power <- vector("double", length = length(n))

    v1 <- 2 * na.omit(lower_crit)
    v2 <- 2 * (n - na.omit(lower_crit) + 1)

    power[is.na(lower_crit)] <- 0
    power[!is.na(lower_crit)] <- pf(v2 * p / (v1 * (1 - p)), v1, v2, lower.tail = F)

  } else{

    if(margin != 0) message("Note: Margin will be ignored for 2-sided test!")

    upper_crit <- purrr::map_dbl(.x = n,
                                 .f = function(x){
                                   above_null <- ceiling(x * null_p):x
                                   valid_cp_l <- qbeta(alpha / 2, above_null, x - above_null + 1) > null_p
                                   if(any(valid_cp_l)){
                                     above_null[which.max(valid_cp_l)]
                                   }else{
                                     NA
                                   }
                                 })

    lower_crit <- purrr::map_dbl(.x = n,
                                 .f = function(x){
                                   below_null <- floor(x * null_p):0
                                   valid_cp_u <- qbeta(alpha / 2, below_null + 1, x - below_null, lower.tail = F) < null_p
                                   if(any(valid_cp_u)){
                                     below_null[which.max(valid_cp_u)]
                                   }else{
                                     NA
                                   }
                                 })

    power_u <- vector("double", length = length(n))
    power_l <- vector("double", length = length(n))
    u_v1 <- 2 * na.omit(upper_crit)
    u_v2 <- 2 * (n - na.omit(upper_crit) + 1)
    l_v1 <- 2 * na.omit(lower_crit)
    l_v2 <- 2 * (n - na.omit(lower_crit) + 1)

    power_u[is.na(upper_crit)] <- 0
    power_u[!is.na(upper_crit)] <- pf(u_v2 * p / (u_v1 * (1 - p)), u_v1, u_v2)
    power_l[is.na(lower_crit)] <- 0
    power_l[!is.na(lower_crit)] <- pf(l_v2 * p / (l_v1 * (1 - p)), l_v1, l_v2, lower.tail = F)

    power <- power_u + power_l
  }

  data.frame(n = n, power = power)

}

#' Calculate Clopper-Pearson confidence interval
#'
#' @description
#' This function calculate the Clopper-Pearson confidence interval.
#'
#' @param p The observed rate.
#' @param n The sample size.
#' @param alpha The two-sided significance level. Default is 0.05.
#'
#' @returns A dataframe of Clopper-Pearson confidence interval.
#'
#' @importFrom stats qbeta
#' @examples
#' # example code
#' clopper_pearson_ci(0:20/20, 20)

clopper_pearson_ci <- function(p, n, alpha = 0.05){

  stopifnot(
    "all n must be integer larger than 1" = all(n == round(n)) && all(n >= 1),
    all(p >= 0),
    all(p <= 1),
    "n and p must have the same length" = length(n) == 1 || length(p) == 1 || length(n) == length(p)
  )

  UCL <- qbeta(alpha / 2, n * p + 1, n - n * p, lower.tail = F)
  LCL <- qbeta(alpha / 2, n * p, n - n * p + 1)

  data.frame(p = p, cp.low = LCL, cp.up = UCL)
}
