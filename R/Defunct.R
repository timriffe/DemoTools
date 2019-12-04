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
#'  \item \code{\link{splitMono}}: This function is now called \code{\link{graduate_mono}}.
#'  \item \code{\link{grabill}}: This function is now called \code{\link{graduate_grabill}}.
#'  \item \code{\link{sprague}}: This function is now called \code{\link{graduate_sprague}}.
#'  \item \code{\link{beers}}: This function is now called \code{\link{graduate_beers}}.
#'  \item \code{\link{splitUniform}}: This function is now called \code{\link{graduate_uniform}}.
#'  \item \code{\link{aomegaMortalityLaws}}: This function is now called \code{\link{lt_ax_closeout}}.
#'  
#'  \item \code{\link{lx2dx}}: This function is now called \code{\link{lt_id_l_d}}.  
#'  \item \code{\link{Lx2Tx}}: This function is now called \code{\link{lt_id_L_T}}.
#'  \item \code{\link{Lx2Tx}}: This function is now called \code{\link{lt_id_L_T}}.
#'  \item \code{\link{Lxlx2Sx}}: This function is now called \code{\link{lt_id_Ll_S}}.
#'  \item \code{\link{qx2lx}}: This function is now called \code{\link{lt_id_q_l}}.
#'  \item \code{\link{mxax2qx}}: This function is now called \code{\link{lt_id_ma_q}}.
#'  \item \code{\link{qxax2mx}}: This function is now called \code{\link{lt_id_qa_m}}.
#'  \item \code{\link{mxorqx2ax}}: This function is now called \code{\link{lt_id_morq_a}}.
#'  \item \code{\link{mxax2qx_Backstop}}: This function is now called \code{\link{lt_id_ma_q_robust}}.
#' }
#' 
#' 
#' @name DemoTools-renamed
NULL

# graduate family

#' #' @export
#' #' @rdname graduate_mono_closeout
#' monoCloseout <- graduate_mono_closeout

#' #' @export
#' #' @rdname graduate_mono
#' splitMono <- graduate_mono

#' #' @export
#' #' @rdname graduate_grabill
#' grabill <- graduate_grabill

#' #' @export
#' #' @rdname graduate_sprague
#' sprague <- graduate_sprague

#' #' @export
#' #' @rdname graduate_beers
#' beers <- graduate_beers

#' #' @export  
#' #' @rdname graduate_uniform
#' splitUniform <- graduate_uniform
# --------------------------------------------------- #

# lifetable related functions:
#' #' @export
#' #' @rdname lt_ax_closeout
#' aomegaMortalityLaws <- lt_ax_closeout

#' #' @export
#' #' @rdname lt_id_l_d
#' lx2dx <- lt_id_l_d

#' #' @export
#' #' @rdname lt_id_L_T
#' Lx2Tx <- lt_id_L_T
#' 
#' 

