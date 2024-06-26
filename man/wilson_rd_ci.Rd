% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wilson_score_without_CC.R
\name{wilson_rd_ci}
\alias{wilson_rd_ci}
\title{Creating Wilson score confidence interval for rate difference without continuous
correction}
\usage{
wilson_rd_ci(
  x1,
  n1,
  x2,
  n2,
  alpha = 0.05,
  check = TRUE,
  output = c("Both", "L", "U")
)
}
\arguments{
\item{x1}{A numeric vector of the numerator of rate1.}

\item{n1}{A numeric vector of the denominator of rate1.}

\item{x2}{A numeric vector of the numerator of rate2.}

\item{n2}{A numeric vector of the denominator of rate2.}

\item{alpha}{The significance level. Default is 0.05 which gives 95\% CI.}

\item{check}{Whether to perform input check using \code{\link[=stopifnot]{stopifnot()}}.

Default is TRUE, use FALSE to save time.}

\item{output}{Specify whether to output the whole interval.

Default is "Both". Use "L" or "U" to only output the lower or upper limits.}
}
\value{
A named list with lower limit (LCI) and/or upper limit (UCI) of the
Wilson score confidence interval of rate1 - rate2.
}
\description{
This function construct the (1 - alpha)\% Wilson score confidence interval for
rate difference without continuity correction.

Either \code{x} or \code{n} (for 1 and 2) can use vector as input to conduct vectorized
operation for multiple rate difference. Unequal length of \code{x} and \code{n}, or unequal
length of rate1 and rate2 will throw an error (unless either of the length equals 1).
}
\section{Note}{
The rate difference is calculated as rate1 - rate2.
}

\examples{
wilson_rd_ci(15, 25, 12, 20)

# Several CIs
wilson_rd_ci(c(15, 16), 25, c(12, 13), 20)

# Several CIs, a series of rate1 compare to a single rate 2
wilson_rd_ci(c(15, 16, 17), 25, 12, 20)
}
\seealso{
\code{\link[=wilson_ci]{wilson_ci()}} for CI of a single rate.
}
