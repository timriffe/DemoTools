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
#'  
#'  \item \code{\link{aomegaMortalityLaws}}: This function is now called \code{\link{lt_a_closeout}}.
#'  \item \code{\link{lx2dx}}: This function is now called \code{\link{lt_id_l_d}}.  
#'  \item \code{\link{Lx2Tx}}: This function is now called \code{\link{lt_id_L_T}}.
#'  \item \code{\link{Lx2Tx}}: This function is now called \code{\link{lt_id_L_T}}.
#'  \item \code{\link{Lxlx2Sx}}: This function is now called \code{\link{lt_id_Ll_S}}.
#'  \item \code{\link{qx2lx}}: This function is now called \code{\link{lt_id_q_l}}.
#'  \item \code{\link{mxax2qx}}: This function is now called \code{\link{lt_id_ma_q}}.
#'  \item \code{\link{mx2qx}}: This function is now called \code{\link{lt_id_ma_q}}.
#'  \item \code{\link{mxax2qx_Backstop}}: This function is now called \code{\link{lt_id_ma_q}}.
#'  \item \code{\link{qxax2mx}}: This function is now called \code{\link{lt_id_qa_m}}.
#'  \item \code{\link{mxorqx2ax}}: This function is now called \code{\link{lt_id_morq_a}}.
#'  \item \code{\link{axUN}}: This function is now called \code{\link{lt_a_un}}.
#'  \item \code{\link{axPAS}}: This function is now called \code{\link{lt_a_pas}}.
#'  \item \code{\link{ax.greville.mortpak}}: This function is now called \code{\link{lt_id_morq_a_greville}}.
#'  \item \code{\link{LTabr}}: This function is now called \code{\link{lt_abridged}}.
#'  \item \code{\link{lt_single_simple}}: This function is now called \code{\link{lt_single_mx}}.
#'  \item \code{\link{extra.mortality}}: This function is now called \code{\link{lt_rule_m_extrapolate}}.
#'  \item \code{\link{geta0CD}}: This function is now called \code{\link{lt_rule_1a0_cd}}.
#'  \item \code{\link{geta1_4CD}}: This function is now called \code{\link{lt_rule_4a1_cd}}.
#'  \item \code{\link{M04_2_D0}}: This function is now called \code{\link{lt_rule_4m0_D0}}.
#'  
#'  \item \code{\link{agesmth}}: This function is now called \code{\link{smooth_age_5}}.
#'  \item \code{\link{T9R5L}}: This function is now called \code{\link{smooth_age_5_feeney}}.
#'  \item \code{\link{carrier_farrag_smth}}: This function is now called \code{\link{smooth_age_5_cf}}.
#'  \item \code{\link{kkn_smth}}: This function is now called \code{\link{smooth_age_5_kkn}}.
#'  \item \code{\link{arriaga_smth}}: This function is now called \code{\link{smooth_age_5_arriaga}}.
#'  \item \code{\link{mav_smth}}: This function is now called \code{\link{smooth_age_5_mav}}.
#'  \item \code{\link{zigzag_smth}}: This function is now called \code{\link{smooth_age_5_zigzag}}.
#'  \item \code{\link{strong_smth}}: This function is now called \code{\link{smooth_age_5_strong}}.
#'  \item \code{\link{united_nations_smth}}: This function is now called \code{\link{smooth_age_5_un}}.
#'  
#'  \item \code{\link{adjustAge}}: This function is now called \code{\link{rescale_vector}}.
#'  
#'  \item \code{\link{Whipple}}: This function is now called \code{\link{check_heaping_whipple}}.
#'  \item \code{\link{Myers}}: This function is now called \code{\link{check_heaping_myers}}.
#'  \item \code{\link{Bachi}}: This function is now called \code{\link{check_heaping_bachi}}.
#'  \item \code{\link{CoaleLi}}: This function is now called \code{\link{check_heaping_coale_li}}.
#'  \item \code{\link{Noumbissi}}: This function is now called \code{\link{check_heaping_noumbissi}}.
#'  \item \code{\link{Spoorenberg}}: This function is now called \code{\link{check_heaping_spoorenberg}}.
#'  \item \code{\link{KannistoHeap}}: This function is now called \code{\link{check_heaping_kannisto}}.
#'  \item \code{\link{Jdanov}}: This function is now called \code{\link{check_heaping_jdanov}}.
#'  \item \code{\link{zero_pref_sawtooth}}: This function is now called \code{\link{check_heaping_sawtooth}}.
#'  \item \code{\link{five_year_roughness}}: This function is now called \code{\link{check_heaping_roughness}}.
#' }
#' 
#' 
#' @name DemoTools-renamed
NULL
