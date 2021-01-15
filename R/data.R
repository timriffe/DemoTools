# Author: IK
# Date: 2018-10-12 UPD 2019-12-04
################################################################################


# #' #' Example data [this is the short explanation to appear in the list of functions]
# #' #'
# #' #' Some more explanation of the data [this will appear on the help page]
# #' #'
# #' #' @format
# #' #'   A data frame with X rows and Y variables:
# #' #'   \describe{
# #' #'     \item{rank}{Rank of the country in a given year.}
# #' #'     \item{country}{Country.}
# #' #'   }
# #' #'
# #' #' @source
# #' #'   The data comes from
# #' #'   \url{https://URL}
# #' "data_object"



#' Indian male population 1991
#'
#' Indian male population 1991
#' @docType data
#' @format
#'   A numeric vector of length 101
#'
#' @source
#'   REGISTRAR GENERAL & CENSUS COMMISSIONER, INDIA, 1976,
#'   Census of India 1971, Series 1, India, Social and Cultural Tables,
#'   Part II-C(ii), New Delhi
"pop1m_ind"


#' Matrix of population over 5 years
#'
#' Matrix of population over 5 years, 5 year age groups, some population 1950-1954
#' @docType data
#' @format
#'   A matrix of 21 rows and 5 columns
#'
#' @source
#'   The data comes from
#'   \url{https://URL}
"pop5_mat"


#' Male population by 5 year age groups
#'
#' Male population by 5 year age groups from PASEX AGESMTH
#' @docType data
#' @format
#'   A numeric vector of length 17
#'
#' @source
#'   The data comes from
#'   \url{https://URL}
"pop5m_pasex"


#' Male population by 1 year age groups
#'
#' Male population by 1 year age groups from PASEX SINGAGE
#' @docType data
#' @format
#'   A numeric vector of length 100
#'
#' @source
#'   The data comes from
#'   \url{https://URL}
"pop1m_pasex"


#' Deaths by 5 year age groups
#'
#' Deaths by 5 year age groups in South Africa 1997, from Feeney Zigzag 2013
#' @docType data
#' @format
#'   A numeric vector of length 20
#'
#' @source
#'   The data comes from
#'   \url{http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/}
"dth5_zigzag"



#' Abridged population from PAS AGEINT -- earlier
#'
#' Abridged population from PAS AGEINT -- earlier, 1980
#' @docType data
#' @format
#'   A numeric vector of length 18
#'
#' @source
#'   The data comes from
#'   \url{http://}
"popA_earlier"


#' Abridged population from PAS AGEINT -- later
#'
#' Abridged population from PAS AGEINT -- later, 1990
#' @docType data
#' @format
#'   A numeric vector of length 18
#'
#' @source
#'   The data comes from
#'   \url{http://}
"popA_later"

#' Russian census 2002 male population by 1 year age groups
#'
#' Male population by 1 year age groups from Russian census help on 2002-10-16
#' @docType data
#' @format
#'   A numeric vector of length 101
#'
#' @source
#'   The data comes from
#'   \url{http://www.demoscope.ru/weekly/ssp/rus2002_01.php}
"pop1m_rus2002"

#' Russian census 2010 male population by 1 year age groups
#'
#' Male population by 1 year age groups from Russian census help on 2010-10-25
#' @docType data
#' @format
#'   A numeric vector of length 101
#'
#' @source
#'   The data comes from
#'   \url{http://www.demoscope.ru/weekly/ssp/rus_age1_10.php}
"pop1m_rus2010"


# model life tables --fitted LogQuad models -------------------------------
 
#' LogQuad model for BOTH SEX fitted for all HMD life tables
#'
#' LogQuad model fitted using `MortalityEstimate::wilmoth` for all BOTH SEX period life tables present in Human Mortality Database (<https://mortality.org>) in December 2019 (968 life tables). Object of class `wilmoth`.
#' @docType data
#' @format
#'   List of 6:
#'   \describe{
#'     \item{input}{List of parameters passed to `MortalityEstimate::wilmoth`.}
#'     \item{call}{R code call line used for fitting.}
#'     \item{coefficients}{Dataframe with 24 rows and 4 variables -- parameters `ax`, `bx`, `cx`, `vx` for each age group.}
#'     \item{k}{`k` fitting parameter of the LogQuag model, a vector of 968 values.}
#'     \item{fitted.values}{Fitted death rates, a matrix 968x24.}
#'     \item{model.info}{Model formula.}
#'  }
#' @source Human Mortality Database. Retrieved 2019-11-28, from <https://mortality.org>
"fitted_logquad_b"

#' LogQuad model for FEMALES fitted for all HMD life tables
#'
#' LogQuad model fitted using \code{MortalityEstimate::wilmoth} for all FEMALE period life tables present in Human Mortality Database (<https://mortality.org>) in December 2019 (968 life tables). Object of class \code{wilmoth}.
#' @docType data
#' @format
#'   List of 6:
#'   \describe{
#'     \item{input}{List of parameters passed to \code{MortalityEstimate::wilmoth}.}
#'     \item{call}{R code call line used for fitting.}
#'     \item{coefficients}{Dataframe with 24 rows and 4 variables -- parameters \code{ax}, \code{bx}, \code{cx}, \code{vx} for each age group.}
#'     \item{k}{\code{k} fitting parameter of the LogQuag model, a vector of 968 values.}
#'     \item{fitted.values}{Fitted death rates, a matrix 968x24.}
#'     \item{model.info}{Model formula.}
#'     }
#' @source Human Mortality Database. Retrieved 2019-11-28, from <https://mortality.org>
"fitted_logquad_f"
#' 
#' LogQuad model for MALES fitted for all HMD life tables
#'
#' LogQuad model fitted using \code{MortalityEstimate::wilmoth} for all MALE period life tables present in Human Mortality Database (\url{https://mortality.org}) in December 2019 (968 life tables). Object of class \code{wilmoth}.
#' @docType data
#' @format
#'   List of 6:
#'   \describe{
#'     \item{input}{List of parameters passed to \code{MortalityEstimate::wilmoth}.}
#'     \item{call}{R code call line used for fitting.}
#'     \item{coefficients}{Dataframe with 24 rows and 4 variables -- parameters \code{ax}, \code{bx}, \code{cx}, \code{vx} for each age group.}
#'     \item{k}{\code{k} fitting parameter of the LogQuag model, a vector of 968 values.}
#'     \item{fitted.values}{Fitted death rates, a matrix 968x24.}
#'     \item{model.info}{Model formula.}
#'     }
#' @source Human Mortality Database. Retrieved 2019-11-28, from <https://mortality.org>
"fitted_logquad_m"