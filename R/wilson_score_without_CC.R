#' Calculating Wilson score confidence interval without continuity correction
#'
#' @description
#' This function construct the (1 - alpha)% Wilson score confidence interval
#' without continuity correction.
#'
#' Either `x` or `n` can use vector as input to conduct vectorized operation for
#' multiple rates. Unequal length of `x` and `n` will throw an error (unless either
#' of the length equals 1).
#'
#' @param x A numeric vector of the numerator of the rate.
#' @param n A numeric vector of the denominator of the rate.
#' @param alpha The significance level. Default is 0.05 which gives 95% CI.
#' @param check Whether to perform input check using [stopifnot()].
#'
#' Default is TRUE, use FALSE to save time.
#'
#' @returns A named list with lower limit (LCI) and upper limit (UCI) of the
#' Wilson score confidence interval without continuity correction.
#'
#' @export
#' @seealso [wilson_rd_ci()] for CI of rate difference.
#' @examples
#' wilson_ci(10, 100)
#'
#' # Several CIs
#' wilson_ci(c(10:15), 100)
#'
#' # Several CIs, different x and n
#' wilson_ci(c(10, 20), c(100, 200))

wilson_ci <- function(x, n, alpha = 0.05, check = TRUE) {
  if (check) {
    stopifnot(
      "x must be integer" = x == round(x),
      "n must be integer" = n == round(n),
      "x and n must have the same length" = length(x) == 1 || length(n) == 1 || length(x) == length(n),
      n > 0,
      x >= 0,
      x <= n,
      alpha > 0,
      alpha < 1,
      is.logical(corrected)
    )
  }

  p <- x / n
  q <- 1 - p
  z <- qnorm(alpha / 2, lower.tail = F)

  UCI <- (2 * n * p + z^2 + z * sqrt(z^2 + 4 * n * p * q)) / (2 * n + 2 * z^2)
  LCI <- (2 * n * p + z^2 - z * sqrt(z^2 + 4 * n * p * q)) / (2 * n + 2 * z^2)

  list(LCI = LCI, UCI = UCI)
}


#' Creating Wilson score confidence interval for rate difference without continuous
#' correction
#'
#' @description
#' This function construct the (1 - alpha)% Wilson score confidence interval for
#' rate difference without continuity correction.
#'
#' Either `x` or `n` (for 1 and 2) can use vector as input to conduct vectorized
#' operation for multiple rate difference. Unequal length of `x` and `n`, or unequal
#' length of rate1 and rate2 will throw an error (unless either of the length equals 1).
#'
#' @details # Note
#' The rate difference is calculated as rate1 - rate2.
#'
#' @param x1 A numeric vector of the numerator of rate1.
#' @param n1 A numeric vector of the denominator of rate1.
#' @param x2 A numeric vector of the numerator of rate2.
#' @param n2 A numeric vector of the denominator of rate2.
#' @param alpha The significance level. Default is 0.05 which gives 95% CI.
#' @param check Whether to perform input check using [stopifnot()].
#'
#' Default is TRUE, use FALSE to save time.
#'
#' @param output Specify whether to output the whole interval.
#'
#' Default is "Both". Use "L" or "U" to only output the lower or upper limits.
#'
#' @returns A named list with lower limit (LCI) and/or upper limit (UCI) of the
#' Wilson score confidence interval of rate1 - rate2.
#'
#' @export
#' @seealso [wilson_ci()] for CI of a single rate.
#' @examples
#' wilson_rd_ci(15, 25, 12, 20)
#'
#' # Several CIs
#' wilson_rd_ci(c(15, 16), 25, c(12, 13), 20)
#'
#' # Several CIs, a series of rate1 compare to a single rate 2
#' wilson_rd_ci(c(15, 16, 17), 25, 12, 20)

wilson_rd_ci <- function(x1, n1, x2, n2, alpha = 0.05, check = TRUE, output = c("Both", "L", "U")) {
  if (check) {
    stopifnot(
      "x1 must be integer" = x1 == round(x1),
      "n1 must be integer" = n1 == round(n1),
      "x2 must be integer" = x2 == round(x2),
      "n2 must be integer" = n2 == round(n2),
      "x1 and n1 must have the same length" = length(x1) == 1 || length(n1) == 1 || length(x1) == length(n1),
      "x2 and n2 must have the same length" = length(x2) == 1 || length(n2) == 1 || length(x2) == length(n2),
      "rate1 and rate2 must have the same length" = (length(x1) == 1 & length(n1) == 1) ||
        (length(x2) == 1 & length(n2) == 1) || max(length(x1), length(n1)) == max(length(x2), length(n2)),
      n1 > 0,
      x1 >= 0,
      x1 <= n1,
      n2 > 0,
      x2 >= 0,
      x2 <= n2,
      alpha > 0,
      alpha < 1
    )
  }

  output <- match.arg(output)

  ci1 <- wilson_ci(x1, n1, alpha, check = FALSE)
  ci2 <- wilson_ci(x2, n2, alpha, check = FALSE)

  p1 <- x1 / n1
  p2 <- x2 / n2

  LCI <- p1 - p2 - sqrt((p1 - ci1[["LCI"]])^2 + (ci2[["UCI"]] - p2)^2)
  UCI <- p1 - p2 + sqrt((ci1[["UCI"]] - p1)^2 + (p2 - ci2[["LCI"]])^2)

  if (output == "Both") {
    list(LCI = LCI, UCI = UCI)
  } else if (output == "L") {
    LCI
  } else {
    UCI
  }
}


#' Calculate the power of a test for rate difference using Wilson score CI
#'
#' @description
#' This function calculate the power of a test for rate difference using Wilson
#' score CI without continuity correction. The function enumerate all possible
#' binary samples of a given sample size. The power is calculated as the average
#' possibility (PMF of binomial distribution) of achieving positive test over all
#' samples.
#'
#' @details # Note
#'
#' The test procedure is actually a one-sided test comparing the lower limit of
#' Wilson score CI of rate difference and delta.
#'
#' @param n A numeric vector of the total sample size, power will be calculated for
#' each n.
#' @param r1 The proportion of sample1 in the total sample. Default is 0.5, which
#' lead to equal sample size for sample1 and sample2.
#' @param p1 The population rate where sample1 comes from.
#' @param p2 The population rate where sample2 comes from.
#' @param alpha The two-sided significance level. Default is 0.05, which is a 0.025
#' one-sided test.
#' @param delta The effect to be tested (|delta| <= 1).
#'
#'    The hypothesis is constructed as:
#'
#'    H0: rate1 - rate2 <= delta
#'
#'    H1: rate1 - rate2 > delta
#'
#' @returns A dataframe of sample size n and power of the test.
#'
#' @importFrom purrr pmap_dbl
#' @importFrom stats qnorm
#' @importFrom magrittr `%>%`
#' @import dplyr
#' @export
#' @seealso [wilson_rd_ci()] for CI of rate difference.
#' @examples
#' # example code
#'
#' power_wilson(n = 792, p1 = 0.76, p2 = 0.76, delta = -0.1)
#'
#' # For multiple n
#' power_wilson(n = seq(300, 500, 10), p1 = 0.99, p2 = 0.99, delta = -0.05)

power_wilson <- function(n, r1 = 0.5, p1, p2, alpha = 0.05, delta) {
  stopifnot(
    "all n must be integer larger than 2" = all(n == round(n)) && all(n >= 2),
    r1 > 0,
    r1 < 1,
    abs(p1) <= 1,
    abs(p2) <= 1,
    abs(delta) <= 1
  )

  n1 <- round(n * r1)
  n2 <- n - n1

  if (any(n1 == 0) || any(n2 == 0)) {
    stop("Allocation ratio r1 generates zero sample size in one group!")
  }

  power <- list(n1, n2) %>% purrr::pmap_dbl(function(n1, n2) {
    sample_x1 <- rep(seq(0, n1, 1), n2 + 1)
    sample_x2 <- rep(seq(0, n2, 1), each = n1 + 1)

    pdf1 <- dbinom(0:n1, n1, p1)
    pdf2 <- dbinom(0:n2, n2, p2)
    pdf <- rep(pdf1, n2 + 1) * rep(pdf2, each = n1 + 1)

    res <- wilson_rd_ci(
      x1 = sample_x1, n1 = n1,
      x2 = sample_x2, n2 = n2,
      alpha = alpha, check = FALSE, output = "L"
    )

    sum(pdf[res > delta])
  })

  data.frame(n = n, power = power)
}
