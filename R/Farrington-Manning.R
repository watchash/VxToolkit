#' Calculate the power of a Farrington-Manning test.
#'
#' @description
#' This function calculate the power of a Farrington-Manning test using normal
#' approximation.
#'
#' @param n A numeric vector of the total sample size, power will be calculated for
#' each n.
#' @param r1 The proportion of sample1 in the total sample. Default is 0.5, which
#' lead to equal sample size for sample1 and sample2.
#' @param p1 The population rate where sample1 comes from.
#' @param p2 The population rate where sample2 comes from.
#' @param alpha The two-sided significance level. Default is 0.05, which is a 0.025
#' one-sided test.
#' @param side The side of the test, "U" and "L" for one-sided test in upper or
#' lower direction. "Both" for a two-sided test. Default is "U".
#' @param delta The effect to be tested (|delta| <= 1).
#'
#' The hypothesis is constructed as:
#'
#' For `side` = "U"
#'
#' H0: rate1 - rate2 <= delta
#'
#' H1: rate1 - rate2 > delta
#'
#' For `side` = "L"
#'
#' H0: rate1 - rate2 >= delta
#'
#' H1: rate1 - rate2 < delta
#'
#' For `side` = "Both"
#'
#' H0: rate1 - rate2 = delta
#'
#' H1: rate1 - rate2 != delta

#' @returns A dataframe of sample size n and power of the test.
#'
#' @references Farrington CP, Manning G. Test statistics and sample size formulae
#' for comparative binomial trials with null hypothesis of non-zero risk difference
#' or non-unity relative risk. Stat Med. 1990 Dec;9(12):1447-54.
#' @importFrom stats qnorm pnorm
#' @export
#' @examples
#' # example code
#' # FM test for a NI margin of -0.1, H1: p1 - p2 > -0.1
#' power_Farrington_Manning(n = 500, p1 = 0.8, p2 = 0.8, delta = -0.1)
#'
#' # For multiple n
#' power_Farrington_Manning(n = seq(400, 600, 2), p1 = 0.8, p2 = 0.8, delta = -0.1)

power_Farrington_Manning <- function(n, r1 = 0.5, p1, p2, alpha = 0.05,
                                     side = c("U", "L", "Both"), delta) {
  stopifnot(
    "all n must be integer larger than 2" = all(n == round(n)) && all(n >= 2),
    r1 > 0,
    r1 < 1,
    abs(p1) <= 1,
    abs(p2) <= 1,
    abs(delta) <= 1,
    "When delta = 0, FM test reduce to standard z-test using pooled variance. Use
    z-test instead, or input a very small delta for approximation" = delta != 0
  )

  side <- match.arg(side)

  n1 <- round(n * r1)
  n2 <- n - n1
  theta <- n2 / n1

  if (any(n1 == 0) || any(n2 == 0)) {
    stop("Allocation ratio r1 generates zero sample size in one group!")
  }

  # MLE estimates of two rates under H0 restriction
  a <- 1 + theta
  b <- -(1 + theta + p1 + theta * p2 + delta * (theta + 2))
  c <- delta^2 + delta * (2 * p1 + theta + 1) + p1 + theta * p2
  d <- -p1 * delta * (1 + delta)

  v <- b^3 / (3 * a)^3 - b * c / (6 * a^2) + d / (2 * a)
  u <- sign(v) * sqrt(b^2 / (3 * a)^2 - c / (3 * a))
  w <- (pi + acos(v / u^3)) / 3

  p1_rmle <- 2 * u * cos(w) - b / (3 * a)
  p2_rmle <- p1_rmle - delta

  sd_h1 <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  sd_rmle <- sqrt(p1_rmle * (1 - p1_rmle) / n1 + p2_rmle * (1 - p2_rmle) / n2)

  z <- qnorm(alpha / 2, lower.tail = F)

  if (side == "U") {
    power <- pnorm((p1 - p2 - delta - z * sd_rmle) / sd_h1)
  } else if (side == "L") {
    power <- pnorm((-(p1 - p2 - delta) - z * sd_rmle) / sd_h1)
  } else {
    power <- pnorm((p1 - p2 - delta - z * sd_rmle) / sd_h1) +
      pnorm((-(p1 - p2 - delta) - z * sd_rmle) / sd_h1)
  }

  data.frame(n = n, power = power)
}
