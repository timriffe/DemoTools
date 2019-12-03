# Function renaming coming soon.
# This page will contain deprecated function names.
# https://ropensci.org/technotes/2017/01/05/package-evolution/


# Ideas: 
# 1) snake case for method/ function families:
# 2) first part denotes family, second function.
# 3) so instead of LTabr() use lt_abridged()
#    lt_single(). For identity functions. mx2qx -> lt_identity_mx_to_qx
# These are more verbose names, but the longer ones are just internal,
# using segmented regular names like this helps with autocomplete in RStudio.

# 4) for the sake of autocomplete, carrier_farrag_smth(), should become smth_carrier_farrag()

# 5) like wise, Sprague -> graduate_sprague, Beers -> graduate_beers()
# 6) Other families: age_util_XXXX, 


#' Renamed functions in DemoTools
#' 
#' These functions still work but will be removed (defunct) in the next version.
#' 
#' \itemize{
#'  \item \code{\link{monoCloseout}}: This function is now called \code{\link{graduate_mono_closeout}}.
#' 
#'  \item \code{\link{splitMono}}: This function is now called \code{\link{graduate_mono}}.
#' 
#'  \item \code{\link{grabill}}: This function is now called \code{\link{graduate_grabill}}.
#' 
#'  \item \code{\link{sprague}}: This function is now called \code{\link{graduate_sprague}}.
#' 
#'  \item \code{\link{beers}}: This function is now called \code{\link{graduate_beers}}.
#'  
#'  \item \code{\link{splitUniform}}: This function is now called \code{\link{graduate_uniform}}.
#'  
#'  \item \code{\link{aomegaMortalityLaws}}: This function is now called \code{\link{lt_ax_closeout}}.
#'   
#'  \item \code{\link{lx2dx}}: This function is now called \code{\link{lt_id_l_d}}.
#'  
#'   \item \code{\link{Lx2Tx}}: This function is now called \code{\link{lt_id_L_T}}.
#' }
#' 
#' 
#' @name DemoTools-renamed
NULL

# graduate family

#' @export
#' @rdname graduate_mono_closeout
monoCloseout <- function(...){
  if (as.character(match.call()[[1]]) == "monoCloseout") {
    warning("please use graduate_mono_closeout() instead of monoCloseout().", call. = FALSE)
  }
  graduate_mono_closeout(...)
}

#' @export
#' @rdname graduate_mono
splitMono <- function(...){
  if (as.character(match.call()[[1]]) == "splitMono") {
    warning("please use graduate_mono() instead of splitMono().", call. = FALSE)
  }
  graduate_mono(...)
}

#' @export
#' @rdname graduate_grabill
grabill <- function(...){
  if (as.character(match.call()[[1]]) == "grabill") {
    warning("please use graduate_grabill() instead of grabill().", call. = FALSE)
  }
  graduate_grabill(...)
}

#' @export
#' @rdname graduate_sprague
sprague <- function(...){
  if (as.character(match.call()[[1]]) == "sprague") {
    warning("please use graduate_sprague() instead of sprague().", call. = FALSE)
  }
  graduate_sprague(...)
}

#' @export
#' @rdname graduate_beers
beers <- function(...){
  if (as.character(match.call()[[1]]) == "beers") {
    warning("please use graduate_beers() instead of beers().", call. = FALSE)
  }
  graduate_beers(...)
}

#' @export
#' @rdname graduate_uniform
splitUniform <- function(...){
  if (as.character(match.call()[[1]]) == "splitUniform") {
    warning("please use graduate_uniform() instead of splitUniform().", call. = FALSE)
  }
  graduate_uniform(...)
}
# --------------------------------------------------- #

# lifetable related functions:
#' @export
#' @rdname lt_ax_closeout
aomegaMortalityLaws <- function(...){
  if (as.character(match.call()[[1]]) == "aomegaMortalityLaws") {
    warning("please use lt_ax_closeout() instead of aomegaMortalityLaws().", call. = FALSE)
  }
  lt_ax_closeout(...)
}

#' @export
#' @rdname lt_id_l_d
lx2dx <- function(...){
  if (as.character(match.call()[[1]]) == "lx2dx") {
    warning("please use lx2dx() instead of lt_id_l_d().", call. = FALSE)
  }
  lt_id_l_d(...)
}

#' @export
#' @rdname lt_id_L_T
Lx2Tx <-function(...){
  if (as.character(match.call()[[1]]) == "Lx2Tx") {
    warning("please use Lx2Tx() instead of lt_id_L_T().", call. = FALSE)
  }
  lt_id_L_T(...)
}