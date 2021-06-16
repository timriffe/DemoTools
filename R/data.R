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

#' Swedish abridged mortality rates
#'
#' Mortality rates in tidy format for each sex in dates 1990-07-01, 2000-07-01, 2010-07-01
#' @docType data
#' @format
#'   A data frame with:
#'   \describe{
#'     \item{Date}{Reference time for the rates estimate.}
#'     \item{Age}{Inferior age for abridged groups. Carefull: last age 100 is not an OAG}
#'     \item{Sex}{Male \code{m} and female \code{m}.}
#'     \item{nMx}{Mortality rates.}
#'     }
#' @source Human Mortality Database. Retrieved 2021-20-01, from <https://mortality.org>
"mA_swe"

#' Swedish life expectancy at birth
#'
#' Life expectancy at birth by sex in tidy format for dates from 1960-07-01 to 2015-07-01 by 5 calendar years.
#' @docType data
#' @format
#'   A data frame with:
#'   \describe{
#'     \item{Date}{Reference time.}
#'     \item{Sex}{Male \code{m} and female \code{m}.}
#'     \item{e0}{Life expectancy at birth.}
#'     }
#' @source Human Mortality Database. Retrieved 2021-20-01, from <https://mortality.org>
"e0_swe"

#' Population matrix for males five year age groups between 1950 and 2050
#'
#' Population matrix for males five year age groups between 1950 and 2050 for
#' unknown country
#' @docType data
#' @format
#'   A matrix of dimensions 21 x 21
#'
#' @source
#'   Migration residual PAS spreadhseet
"pop_m_mat_five"

#' Population matrix for females five year age groups between 1950 and 2050
#'
#' Population matrix for females five year age groups between 1950 and 2050 for
#' unknown country
#' @docType data
#' @format
#'   A matrix of dimensions 21 x 21
#'
#' @source
#'   Migration residual PAS spreadhseet
"pop_f_mat_five"

#' Survival rates matrix for males five year age groups between 1950 and 2045
#'
#' Survival rates matrix for males five year age groups between 1950 and 2045
#' for unknown country
#' @docType data
#' @format
#'   A matrix of dimensions 21 x 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"sr_m_mat_five"

#' Survival rates matrix for females five year age groups between 1950 and 2045
#'
#' Survival rates matrix for females five year age groups between 1950 and 2045
#' for unknown country
#' @docType data
#' @format
#'   A matrix of dimensions 21 x 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"sr_f_mat_five"

#' Age-specific fertility rates for age groups 15 to 45 between 1950 and 2045
#'
#' Age-specific fertility rates for age groups 15 to 45 between 1950 and 2045
#' for unknown country
#' @docType data
#' @format
#'   A matrix of dimensions 7 x 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"asfr_mat_five"

#' Sex ratio at birth between 1950 and 2045
#'
#' Sex ratio at birth between 1950 and 2045 for unknown country
#' @docType data
#' @format
#'   A vector of length 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"srb_vec_five"

#' Ages between 0 and 100 abridged in five year age groups
#'
#' Ages between 0 and 100 abridged in five year age groups for unknown
#' country
#' @docType data
#' @format
#'   A vector of length 21
#'
#' @source
#'   Migration residual PAS spreadhseet
"ages_five"

#' Ages between 15 and 45 in five year age groups
#'
#' Ages between 15 and 45 in five year age groups for unknown
#' country
#' @docType data
#' @format
#'   A vector of length 7
#'
#' @source
#'   Migration residual PAS spreadhseet
"ages_asfr_five"

#' Population matrix for males single ages between 1999 and 2019
#'
#' Population matrix for males single ages between 1999 and 2019 for
#' Sweden
#' @docType data
#' @format
#'   A matrix of dimensions 101 x 21
#'
#' @source
#'   Migration residual PAS spreadhseet
"pop_m_mat_single"

#' Population matrix for females single ages between 1999 and 2019
#'
#' Population matrix for females single ages between 1999 and 2019 for
#' Sweden
#' @docType data
#' @format
#'   A matrix of dimensions 101 x 21
#'
#' @source
#'   Migration residual PAS spreadhseet
"pop_f_mat_single"

#' Survival rates matrix for males single ages between 1999 and 2019
#'
#' Survival rates matrix for males single ages between 1999 and 2019 for
#' Sweden
#'
#' @docType data
#' @format
#'   A matrix of dimensions 101 x 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"sr_m_mat_single"

#' Survival rates matrix for females single ages between 1999 and 2019
#'
#' Survival rates matrix for females single ages between 1999 and 2019 for
#' Sweden
#'
#' @docType data
#' @format
#'   A matrix of dimensions 101 x 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"sr_f_mat_single"

#' Age-specific fertility rates for single ages 15 to 49 between 1999 and 2018
#'
#' Age-specific fertility rates for single ages 15 to 49 between 1999 and 2018
#' for Sweden
#'
#' @docType data
#' @format
#'   A matrix of dimensions 35 x 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"asfr_mat_single"

#' Sex ratio at birth between 1999 and 2019
#'
#' Sex ratio at birth between 1999 and 2019 for Sweden
#'
#' @docType data
#' @format
#'   A vector of length 20
#'
#' @source
#'   Migration residual PAS spreadhseet
"srb_vec_single"

#' Single ages between 0 and 100
#'
#' Single ages between 0 and 100 for Sweden, 1999-2019.
#' @docType data
#' @format
#'   A vector of length 101
#'
#' @source
#'   Migration residual PAS spreadhseet
"ages_single"

#' Single ages between 15 and 49
#'
#' Single ages between 15 and 49 for Sweden
#' @docType data
#' @format
#'   A vector of length 36
#'
#' @source
#'   Migration residual PAS spreadhseet
"ages_asfr_single"

#' Parameters for considered migration profiles
#'
#' Roger-Castro estimated parameters using `mig_estimate_rc` for Pre Working Age and Working Age profiles of migration.
#' @docType data
#' @format
#'   A data frame with:
#'   \describe{
#'     \item{family}{Types Family, Male Labor or Female Labor.}
#'     \item{sex}{Male and Female.}
#'     \item{mig_sign}{Inmigration or Emigration.}
#'     \item{param}{Parameters from Roger-Castro.}
#'     \item{median}{median of posterior distribution using Monte Carlo Markov Chains in `mig_estimate_rc`.}
#'     }
#' @source UN spreadsheet "UNPD_Migration Age Profiles.xlsx"
"mig_un_params"

#' Proportion of net migrants by age and sex for considered migration profiles
#'
#' Roger-Castro estimated proportion of total net migrants using parameters from `mig_un_params` data.
#' @docType data
#' @format
#'   A data frame with:
#'   \describe{
#'     \item{family}{Types Family, Male Labor or Female Labor.}
#'     \item{sex}{Male and Female.}
#'     \item{mig_sign}{Inmigration or Emigration.}
#'     \item{age}{Simple ages from 0 to 80 (OAG).}
#'     \item{prop}{Proportion of net migrants due to that sex and age.}
#'     }
#' @source UN spreadsheet "UNPD_Migration Age Profiles.xlsx"
"mig_un_families"