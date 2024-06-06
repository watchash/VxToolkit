#' Calculate the power of a two sample t-test assuming equal variance for log-normal data
#'
#' @description
#' This function calculate the power of a t-test for log-normal data (e.g., titre).
#' Equal variance is assumed.
#'
#' @param n A numeric vector of the total sample size, power will be calculated for
#' each n.
#' @param r1 The proportion of sample1 in the total sample. Default is 0.5, which
#' lead to equal sample size for sample1 and sample2.
#' @param diff_log10 The true difference of log10 mean.
#' @param sd_log10 The standard deviation of log10 data.
#' @param alpha The two-sided significance level. Default is 0.05, which is a 0.025
#' one-sided test.
#' @param side The side of the test, "U" and "L" for one-sided test in upper or
#' lower direction. "Both" for a two-sided equivalence test.
#'
#' Default is "U".
#'
#' @param margin The effect to be tested.
#'
#' The hypothesis is constructed as:
#'
#' For `side` = "U"
#'
#' H0: geomean1 / geomean2 <= 1 / margin
#'
#' H1: geomean1 / geomean2 > 1 / margin
#'
#' For `side` = "L"
#'
#' H0: geomean1 / geomean2 >= 1 / margin
#'
#' H1: geomean1 / geomean2 < 1 / margin
#'
#' For `side` = "Both"
#'
#' H0: geomean1 / geomean2 = 1 / margin
#'
#' H1: geomean1 / geomean2 != 1 / margin
#'
#' @returns A dataframe of sample size n and power of the test.
#'
#' @importFrom stats qt pt qf pf
#' @export
#' @examples
#' # example code
#' power_ttest_log(400, diff_log10 = 0, sd_log10 = 0.5, margin = 1.5)
#'
#' # For multiple n
#' power_ttest_log(seq(300, 400, 2), diff_log10 = 0, sd_log10 = 0.5, margin = 1.5)

power_ttest_log <- function(n, r1 = 0.5, diff_log10 = 0, sd_log10, alpha = 0.05,
                            side = c("U", "L", "Both"), margin){

  stopifnot(
    "all n must be integer larger than 3" = all(n == round(n)) && all(n >= 3),
    r1 > 0,
    r1 < 1,
    sd_log10 > 0,
    margin > 0
  )

  side <- match.arg(side)

  n1 <- round(n * r1)
  n2 <- n - n1

  delta <- sqrt(n) * sqrt(n1 * n2 / n^2) * (diff_log10 - log10(1 / margin)) / sd_log10
  df <- n - 2

  if (side %in% c("L", "U")) {

    tcrit <- qt(alpha / 2, df = df, lower.tail = F)
    power <-  switch (side,
                      L = pt(q = tcrit, ncp = delta, df = df),
                      U = pt(q = tcrit, ncp = delta, df = df, lower.tail = F))
    } else {

    Fcrit <- qf(alpha, 1, df1 = 1, df2 = df, lower.tail = F)
    power <- pf(q = Fcrit, ncp = delta^2, df1 = 1, df2 = df, lower.tail = F)
  }

  data.frame(n = n, power = power)
}


#' Calculate the power of a one sample t-test for log-normal data
#'
#' @description
#' This function calculate the power of a one sample t-test for log-normal data
#' (e.g., titre).
#'
#' @param n A numeric vector of the sample size for the single group, power will
#' be calculated for each n.
#' @param diff_log10 The true difference of log10 mean and log10 benchmark.
#' @param sd_log10 The standard deviation of log10 data.
#' @param alpha The two-sided significance level. Default is 0.05, which is a 0.025
#' one-sided test.
#' @param side The side of the test, "U" and "L" for one-sided test in upper or
#' lower direction. "Both" for a two-sided test. Default is "U".
#' @param margin The effect to be tested. Default = 1.
#'
#' Recommend to adjust `diff_log10` directly and not change `margin`.
#'
#' The hypothesis is constructed as:
#'
#' For `side` = "U"
#'
#' H0: geomean / benchmark <= 1 / margin
#'
#' H1: geomean / benchmark > 1 / margin
#'
#' For `side` = "L"
#'
#' H0: geomean / benchmark >= 1 / margin
#'
#' H1: geomean / benchmark < 1 / margin
#'
#' For `side` = "Both"
#'
#' H0: geomean / benchmark = 1 / margin
#'
#' H1: geomean / benchmark != 1 / margin
#'
#' @returns A dataframe of sample size n and power of the test.
#'
#' @importFrom stats qt pt qf pf
#' @export
#' @examples
#' # example code
#' power_ttest_log_onesample(200, diff_log10 = 0, sd_log10 = 0.5)
#'
#' # For multiple n
#' power_ttest_log_onesample(seq(150, 200, 2), diff_log10 = 0, sd_log10 = 0.5)

power_ttest_log_onesample <- function(n, diff_log10, sd_log10, alpha = 0.05,
                                      side = c("U", "L", "Both"), margin = 1){

  stopifnot(
    "all n must be integer larger than 2" = all(n == round(n)) && all(n >= 2),
    sd_log10 > 0
  )

  side <- match.arg(side)

  delta <- sqrt(n) * (diff_log10 - log10(1 / margin)) / sd_log10
  df <- n - 1

  if (side %in% c("L", "U")) {

    tcrit <- qt(alpha / 2, df = df, lower.tail = F)
    power <-  switch (side,
                      L = pt(q = tcrit, ncp = delta, df = df),
                      U = pt(q = tcrit, ncp = delta, df = df, lower.tail = F))
  } else {

    Fcrit <- qf(alpha, 1, df1 = 1, df2 = df, lower.tail = F)
    power <- pf(q = Fcrit, ncp = delta^2, df1 = 1, df2 = df, lower.tail = F)
  }

  data.frame(n = n, power = power)
}
