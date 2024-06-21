#' Calculate the power of a Vx efficacy trial using exact method with fixed total
#' case
#'
#' @description
#' This function calculate the power of a Vx efficacy trial.
#'
#' @param case A numeric vector of the total cases, power will be calculated for
#' each case.
#' @param r1 The proportion of sample1 in the total sample. Default is 0.5, which
#' lead to equal sample size for sample1 and sample2.
#' @param ve The assumed vaccine efficacy in alternative hypothesis. e.g., ve = 0.7
#' means 70% efficacy (rate ratio = 1 - 0.7 = 0.3).
#' @param alpha The two-sided significance level. Default is 0.05, which is a 0.025
#' one-sided test.
#' @param incidence0 The incidence rate in sample2(control group), which is used to
#' convert cases to sample size.
#' @param margin The ve margin to be tested.
#'
#' The hypothesis is constructed as:
#'
#' H0: ve <= margin
#'
#' H1: ve > margin
#'
#' @returns A dataframe of number of cases, sample size and power of the test.
#'
#' @importFrom stats qbeta dbinom
#' @export
#' @examples
#' # example code
#' power_exact_ve(case = 80, r1 = 0.5, ve = 0.7, incidence0 = 0.025, margin = 0.3)
#'
#' # For multiple cases
#' power_exact_ve(case = 50:100, r1 = 0.5, ve = 0.7, incidence0 = 0.025, margin = 0.3)

power_exact_ve <- function(case, r1 = 0.5, ve, alpha = 0.05, incidence0, margin){
  stopifnot(
    "all case must be integer larger than 2" = all(case == round(case)) && all(case >= 2),
    r1 > 0,
    r1 < 1,
    incidence0 > 0,
    incidence0 <= 1,
    ve <= 1,
    margin <= 1
  )

  # probability of observe case in sample1 conditional on total case number
  prob_from1 <- (1 - ve) / (1 - ve + (1 - r1) / r1)

  case_all <- unlist(sapply(case, function(x) rep(x, x + 1)))
  case_1 <- unlist(sapply(case, function(x) 0:x))
  UCI_ratio <- qbeta(alpha / 2, case_1 + 1, case_all - case_1, lower.tail = F)
  LCI_ve <- 1 - UCI_ratio / (1 - UCI_ratio)

  pdf <- dbinom(x = case_1, size = case_all, prob = prob_from1)
  pdf_good <- (LCI_ve > margin) * pdf

  power <- tapply(pdf_good, case_all, sum)
  n0 <- ceiling(case * (1 - prob_from1) / incidence0)
  n1 <- ceiling(n0 * r1 / (1 - r1))
  n <- n0 + n1

  data.frame(case = case, power = power, n1 = n1, n0 = n0, n = n)
}


