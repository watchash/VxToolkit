% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttest_logscale.R
\name{power_ttest_log_onesample}
\alias{power_ttest_log_onesample}
\title{Calculate the power of a one sample t-test for log-normal data}
\usage{
power_ttest_log_onesample(
  n,
  diff_log10,
  sd_log10,
  alpha = 0.05,
  side = c("U", "L", "Both"),
  margin = 1
)
}
\arguments{
\item{n}{A numeric vector of the sample size for the single group, power will
be calculated for each n.}

\item{diff_log10}{The true difference of log10 mean and log10 benchmark.}

\item{sd_log10}{The standard deviation of log10 data.}

\item{alpha}{The two-sided significance level. Default is 0.05, which is a 0.025
one-sided test.}

\item{side}{The side of the test, "U" and "L" for one-sided test in upper or
lower direction. "Both" for a two-sided test. Default is "U".}

\item{margin}{The effect to be tested. Default = 1.

Recommend to adjust \code{diff_log10} directly and not change \code{margin}.

The hypothesis is constructed as:

For \code{side} = "U"

H0: geomean / benchmark <= 1 / margin

H1: geomean / benchmark > 1 / margin

For \code{side} = "L"

H0: geomean / benchmark >= 1 / margin

H1: geomean / benchmark < 1 / margin

For \code{side} = "Both"

H0: geomean / benchmark = 1 / margin

H1: geomean / benchmark != 1 / margin}
}
\value{
A dataframe of sample size n and power of the test.
}
\description{
This function calculate the power of a one sample t-test for log-normal data
(e.g., titre).
}
\examples{
# example code
power_ttest_log_onesample(200, diff_log10 = 0, sd_log10 = 0.5)

# For multiple n
power_ttest_log_onesample(seq(150, 200, 2), diff_log10 = 0, sd_log10 = 0.5)
}
