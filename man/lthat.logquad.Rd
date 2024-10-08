% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lt_model_lq.R
\name{lthat.logquad}
\alias{lthat.logquad}
\title{Estimated life table using the log-quadratic model}
\usage{
lthat.logquad(
  coefs,
  x,
  q0_5,
  k,
  radix,
  axmethod = "pas",
  a0rule = "ak",
  Sex = "m",
  IMR = NA,
  region = "w",
  mod = TRUE,
  SRB = 1.05
)
}
\arguments{
\item{coefs}{Estimated coefficients}

\item{q0_5}{5q0. The probability that a new-born will die during the
subsequent 5 years;}

\item{radix}{Life table radix. Default: 10^5;}

\item{axmethod}{character. Either \code{"pas"} or \code{"un"}.}

\item{a0rule}{character. Either \code{"ak"} (default) or \code{"cd"}.}

\item{Sex}{Choose the sex of the population. This choice defines the use
of a corresponding Log-Quadratic (\code{wilmoth})
model fitted for the whole Human Mortality Database (as of Dec 2019,
there are 968 life tables for each sex).
The following options are available: \itemize{
\item{\code{"b"}} -- Both sex;
\item{\code{"f"}} -- Females;
\item{\code{"m"}} -- Males.
}}

\item{IMR}{numeric. Infant mortality rate \ifelse{html}{\out{q<sub>0}}{\eqn{q_0}}, in case available and \code{nqx} is not specified. Default \code{NA}.}

\item{region}{character. North, East, South, or West: \code{"n"}, \code{"e"}, \code{"s"}, \code{"w"}. Default \code{"w"}.}

\item{mod}{logical. If \code{"un"} specified for \code{axmethod}, whether or not to use Nan Li's modification for ages 5-14. Default \code{TRUE}.}

\item{SRB}{the sex ratio at birth (boys / girls), default 1.05}
}
\description{
Estimated life table using the log-quadratic model
}
\keyword{internal}
