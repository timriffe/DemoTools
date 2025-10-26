#' Construct Single-Year Base Population Using Reverse Survival
#'
#' Constructs a single-year, age-specific population for both sexes using the reverse
#' survival method. This combines life table data, fertility rates, and sex ratio
#' at birth (SRB) to estimate age-specific population counts, optionally using
#' user-supplied demographic data or WPP inputs.
#'
#' @param location Optional character string naming the location (used only for messages).
#' @param refDate Numeric reference date (e.g., census year). Required.
#' @param Age Optional numeric vector of single ages. Must match population vector length if supplied.
#' @param country_code Numeric ISO country code. Required for WPP data retrieval.
#' @param Females_single Optional numeric vector of female single-age population counts.
#' @param Males_single Optional numeric vector of male single-age population counts.
#' @param nLxFemale Optional life table \code{nLx} data for females.
#' @param nLxMale Optional life table \code{nLx} data for males (currently not used directly).
#' @param nLxDatesIn Optional numeric vector of years corresponding to life table data.
#' @param AsfrMat Optional matrix of age-specific fertility rates (ASFR).
#' @param AsfrDatesIn Optional numeric vector of years corresponding to ASFR data.
#' @param SRB Optional numeric vector or single value giving sex ratio at birth (males per female).
#' @param SRBDatesIn Optional numeric vector of years corresponding to SRB data.
#' @param radix Numeric starting radix for life table calculations. Defaults to 1.
#' @param verbose Logical; if \code{TRUE}, prints progress messages. Default = \code{TRUE}.
#' @param ... Additional arguments passed to internal DemoTools functions.
#'
#' @details
#' The function performs the following sequence:
#' \enumerate{
#'   \item Validates input arguments and loads the most recent \code{wpp} dataset installed.
#'   \item Retrieves population, mortality, fertility, and SRB data from WPP if not provided by the user.
#'   \item Constructs survival ratios (\code{SxRev}) from life tables to enable reverse survival.
#'   \item Estimates female and male exposures by age for reproductive ages (15–49).
#'   \item Computes births by year and allocates them to sex using the SRB.
#'   \item Reverse-survives infants and young ages to reconstruct the base population.
#'   \item Returns a list containing reconstructed populations and visual diagnostics.
#' }
#'
#' The reverse survival method is useful when historical population data are unavailable
#' or incomplete but life table and fertility information are known.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{pop_hat_m}{A tibble containing reconstructed single-year population by age, for male counts.}
#'   \item{pop_hat_f}{A tibble containing reconstructed single-year population by age, for female counts.}
#'   \item{figureF}{A \code{ggplot2} line plot comparing estimated and observed female populations (ages 0–9).}
#'   \item{figureM}{A \code{ggplot2} line plot comparing estimated and observed male populations (ages 0–9).}
#' }
#'
#' @importFrom dplyr filter select mutate arrange group_by ungroup right_join left_join summarize
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_line
#' @importFrom magrittr %>%
#' @importFrom DemoTools download_SRB lt_infer_radix_from_1L0 age2int lt_id_q_l lt_rule_1a0_ak lt_id_ma_q lt_id_q_l download_Asfr
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Example: Estimate base population for Germany, year 2000
#' basepop <- basepop_single(
#'   country_code = 276,
#'   refDate = 2000
#' )
#'
#' basepop$pop_hat_m  # Population tibble for men
#' basepop$figureF  # Female population plot
#' basepop$figureM  # Male population plot
#' }
#'
#' @export
#'

basepop_single <- function(location       = NULL,
                           refDate        = NULL,
                           Age            = NULL,
                           country_code   = NULL,
                           Females_single = NULL,
                           Males_single   = NULL,
                           nLxFemale      = NULL,
                           nLxMale        = NULL,
                           nLxDatesIn     = NULL,
                           AsfrMat        = NULL,
                           AsfrDatesIn    = NULL,
                           SRB            = NULL,
                           SRBDatesIn     = NULL,
                           radix          = 1,
                           verbose        = TRUE,
                           ...) {
  
  # --- Global setup ---------------------------------------------------------
  # options from DemoTools version
  options(basepop_verbose = verbose)
  on.exit(options(basepop_verbose = NULL))
  
  # Jan 1 2000 female pop;
  # note:  1999 means dec 31, 1999, so we treat as jan 1, 2000.
  # we need these dates to filter period
  # refDate is treated as Jan 1 of the following year for consistency with WPP
  refDate       <- dec.date(refDate) - 1
  refDate_start <- refDate - 9
  
  # --- Basic validation -----------------------------------------------------
  # we need minimum of date and country to run a function
  if(is.null(refDate) | is.null(country_code)) {
    
    stop("At least refDate and country_code should be provided for function to work")
    
  }
  
  
  # here we attach the latest wpp package available
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  if(length(installed_wpp) == 0) {
    
    stop("No wpp package installed.")
    
  }
  
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  suppressPackageStartupMessages(library(latest_wpp, character.only = TRUE))
  
  if(parse_number(latest_wpp) < 2022) { 
    
    warning("No single ages are availabe in wpp versions earlier that wpp2022. Consider updating the wpp package or change to five year solution.")
    
  }
  
  # --- Load population and mortality data ----------------------------------
  # if user did not provide population we can calculate it from wpp
  if(is.null(Females_single) | is.null(Males_single)) {
    
    data("popAge1dt", package = latest_wpp)
    
  }
  
  # this data is for future calculations
  data("mx1dt", package = latest_wpp)
  
  # --- Prepare user-supplied population vectors ----------------------------
  # Convert numeric vectors to tidy tibbles with cohort tagging
  # if user provided Population vector we turn it into tibble
  if(!is.null(Females_single)) {
    
    Females_single <- tibble(pop    = Females_single,
                             year   = refDate + 1,
                             age    = row_number() - 1,
                             cohort = year - age - 1)
    
    
  }
  
  # same for Males
  if(!is.null(Males_single)) {
    
    Males_single <- tibble(pop    = Males_single,
                           year   = refDate + 1,
                           age    = row_number() - 1,
                           cohort = year - age - 1)
    
  }
  
  # --- Retrieve WPP population if user did not provide one -----------------
  # if user did not provide Population
  if(is.null(Females_single)) {
    
    Females_single <- popAge1dt |>
      filter(country_code == !!country_code, year == refDate) |>
      select(year, age, pop = popF) |>
      mutate(year   = year + 1,
             cohort = year - age - 1)
    
  }
  
  # same operation for males
  if(is.null(Males_single)) {
    
    Males_single <- popAge1dt |>
      filter(country_code == !!country_code, year == refDate) |>
      select(year, age, pop = popM) |>
      mutate(year   = year + 1,
             cohort = year - age - 1)
    
  }
  
  # --- Age validation and setup --------------------------------------------
  # Age setup
  if(!is.null(Age)) {
    
    stopifnot(is_single(Age))
    stopifnot(length(Age) == nrow(Females_single))
    
  } else {
    
    Age <- as.integer(sort(unique(Females_single$age)))
    
  }
  
  # --- Determine life table radix ------------------------------------------
  if(is.null(radix) & !is.null(nLxFemale)) {
    
    radix <- lt_infer_radix_from_1L0(nLxFemale[1, 1])
    
    if(verbose) {
      
      cat(paste0("Setting radix to value of lx: ", radix, 
                 ". Can be overwritten with the `radix` argument"),
          sep = "\n")
      
    }
    
  }
  
  
  if(is.null(radix) & is.null(nLxFemale)) { 
    
    radix <- 1
    cat(paste0("Setting radix to value of lx: ", radix, 
               ". Can be overwritten with the `radix` argument"),
        sep = "\n")
    
  }
  
  # --- Compute reverse survival ratios (SxRev) ------------------------------
  # These represent the probability of surviving backward in time.
  # Female and male SxRev matrices are created either from nLx input or WPP mx data.
  # Inflation factor is cumulative survival ratio (used for exposure reconstruction).

  if(!is.null(nLxFemale)) {
    
    mxF <- nLxFemale %>% 
      group_by(year) %>% 
      mutate(n  = age2int(.data$age, OAvalue = Inf),
             ax = ifelse(age == 0, 
                         solve_a0_from_L0(L0_target = .data$Lxp[.data$age == 0],
                                          l0        = radix,
                                          sex       = "f"), 
                         0.5)) %>% 
      unnest(ax) %>% 
      mutate(qx      = (Lxp - lead(Lxp)) / Lxp,
             lx      = lt_id_q_l(nqx = qx, radix = radix),
             dx      = lx * qx,
             Lx      = lx - dx * (1 - ax),
             mx      = lt_id_qa_m(nqx = qx, nax = ax, AgeInt = n),
             cohort  = year - age - 1,
             qx      = ifelse(age == max(age), 1, qx),
             mx      = ifelse(age == max(age), 1 / ax[age == max(age)], mx)) %>% 
      arrange(cohort, age) |>
      group_by(cohort) |>
      mutate(
        lx    = lt_id_q_l(nqx = qx, radix = radix),
        dx    = lx * qx,
        Lx    = lx - dx * (1 - ax),
        SxRev = Lx / lead(Lx, default = 1)
      ) |>
      ungroup() |>
      arrange(year, age) |>
      group_by(year) |>
      mutate(
        lxp   = lt_id_q_l(nqx = qx, radix = radix),
        dxp   = lxp * qx,
        Lxp   = lxp - dxp * (1 - ax),
        SxRev = ifelse(year == max(year), Lxp / lead(Lxp), SxRev)
      ) |>
      ungroup() |>
      filter(age < 100) |>
      select(cohort, year, age, SxRev) |>
      arrange(cohort, -age) |>
      group_by(cohort) |>
      # inflation for reverse-surviving
      mutate(inflation_factor = cumprod(SxRev))
    
  }
  
  if(!is.null(nLxMale)) {
    
    mxM <- nLxMale %>% 
      group_by(year) %>% 
      mutate(n  = age2int(.data$age, OAvalue = Inf),
             ax = ifelse(age == 0, 
                         solve_a0_from_L0(L0_target = .data$Lxp[.data$age == 0],
                                          l0        = radix,
                                          sex       = "m"), 
                         0.5)) %>% 
      unnest(ax) %>% 
      mutate(qx      = (Lxp - lead(Lxp)) / Lxp,
             lx      = lt_id_q_l(nqx = qx, radix = radix),
             dx      = lx * qx,
             Lx      = lx - dx * (1 - ax),
             mx      = lt_id_qa_m(nqx = qx, nax = ax, AgeInt = n),
             cohort  = year - age - 1,
             qx      = ifelse(age == max(age), 1, qx),
             mx      = ifelse(age == max(age), 1 / ax[age == max(age)], mx)) %>% 
      arrange(cohort, age) |>
      group_by(cohort) |>
      mutate(
        lx    = lt_id_q_l(nqx = qx, radix = radix),
        dx    = lx * qx,
        Lx    = lx - dx * (1 - ax),
        SxRev = Lx / lead(Lx, default = 1)
      ) |>
      ungroup() |>
      arrange(year, age) |>
      group_by(year) |>
      mutate(
        lxp   = lt_id_q_l(nqx = qx, radix = radix),
        dxp   = lxp * qx,
        Lxp   = lxp - dxp * (1 - ax),
        SxRev = ifelse(year == max(year), Lxp / lead(Lxp), SxRev)
      ) |>
      ungroup() |>
      filter(age < 100) |>
      select(cohort, year, age, SxRev) |>
      arrange(cohort, -age) |>
      group_by(cohort) |>
      # inflation for reverse-surviving
      mutate(inflation_factor = cumprod(SxRev))
    
  }
  
  # if not provided by user
  if(is.null(nLxMale)) {
    
    mxM <- mx1dt |>
      filter(country_code == !!country_code,
             between(year, refDate_start, refDate)) |>
      as_tibble() |>
      select(year, age, mx = mxM) |>
      # could try to warp to PC shape here,
      # but uncertain infants. Maybe using
      # an a0 assumption it'd be doable.
      # need cohorts to structure reverse survival
      mutate(
        cohort = year - age - 1,
        age_int = 1,
        ax = if_else(age == 0, lt_rule_1a0_ak(M0 = mx, Sex = "m"), 0.5),
        qx = lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int)
      ) |>
      arrange(cohort, age) |>
      group_by(cohort) |>
      mutate(
        lx = lt_id_q_l(nqx = qx, radix = radix),
        dx = lx * qx,
        Lx = lx - dx * (1 - ax),
        SxRev = Lx / lead(Lx, default = 1)
      ) |>
      ungroup() |>
      arrange(year, age) |>
      group_by(year) |>
      mutate(
        lxp = lt_id_q_l(nqx = qx, radix = radix), 
        dxp = lxp * qx,
        Lxp = lxp - dxp * (1 - ax),
        SxRev = ifelse(year == max(year), Lxp / lead(Lxp), SxRev)
      ) |>
      ungroup() |>
      filter(age < 100) |>
      select(cohort, year, age, SxRev) |>
      arrange(cohort, -age) |>
      group_by(cohort) |>
      # inflation for reverse-surviving
      mutate(inflation_factor = cumprod(SxRev))
  }
  
  # if not provided by user
  if(is.null(nLxFemale)) {
    
    mxF <- mx1dt |>
      filter(country_code == !!country_code,
             between(year, refDate_start, refDate)) |>
      as_tibble() |>
      select(year, age, mx = mxF) |>
      # could try to warp to PC shape here,
      # but uncertain infants. Maybe using
      # an a0 assumption it'd be doable.
      # need cohorts to structure reverse survival
      mutate(
        cohort = year - age - 1,
        age_int = 1,
        ax = if_else(age == 0, lt_rule_1a0_ak(M0 = mx, Sex = "f"), 0.5),
        qx = lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int)
      ) |>
      arrange(cohort, age) |>
      group_by(cohort) |>
      mutate(
        lx = lt_id_q_l(nqx = qx, radix = radix),
        dx = lx * qx,
        Lx = lx - dx * (1 - ax),
        SxRev = Lx / lead(Lx, default = 1)
      ) |>
      ungroup() |>
      arrange(year, age) |>
      group_by(year) |>
      mutate(
        lxp = lt_id_q_l(nqx = qx, radix = radix), 
        dxp = lxp * qx,
        Lxp = lxp - dxp * (1 - ax),
        SxRev = ifelse(year == max(year), Lxp / lead(Lxp), SxRev)
      ) |>
      ungroup() |>
      filter(age < 100) |>
      select(cohort, year, age, SxRev) |>
      arrange(cohort, -age) |>
      group_by(cohort) |>
      # inflation for reverse-surviving
      mutate(inflation_factor = cumprod(SxRev))
  }
  
  
  # --- Compute exposures ---------------------------------------------------
  # Exposures are estimated as mean of adjacent populations (mid-year approx).
  # Restricted to ages 15–49 for fertility-related computations.
  
  # (Exposure calculations for expF and expM unchanged)
  
  expF <- Females_single |>
    select(cohort, pop) |>
    right_join(mxF, by = join_by(cohort)) |>
    mutate(pop_hat = pop * inflation_factor) |>
    select(year, age, pop = pop_hat) |>
    bind_rows(Females_single |> select(year, age, pop)) |>
    filter(between(age, 15, 50)) |>
    arrange(age, year) |>
    mutate(pop_1p1  = lead(pop),
           exposure = (pop + pop_1p1) / 2) |>
    filter(age < 50)
  
  expM <- Males_single |>
    select(cohort, pop) |>
    right_join(mxM, by = join_by(cohort)) |>
    mutate(pop_hat = pop * inflation_factor) |>
    select(year, age, pop = pop_hat) |>
    bind_rows(Females_single |> select(year, age, pop)) |>
    filter(between(age, 15, 50)) |>
    arrange(age, year) |>
    mutate(pop_1p1  = lead(pop),
           exposure = (pop + pop_1p1) / 2) |>
    filter(age < 50)
  
  
  # --- Fertility rates (ASFR) ----------------------------------------------
  # If user supplies AsfrMat, reshape to long format; otherwise, download from WPP.
  if(!is.null(AsfrMat)) {
    
    AsfrMat <- AsfrMat %>%
      pivot_longer(-c(country_code, name, age),
                   names_to  = "year",
                   values_to = "asfr") %>%
      mutate(year = as.numeric(year)) %>%
      select(-name)
    
    
  }
  
  if(is.null(AsfrMat)) {
    
    AsfrMat <- download_Asfr(Asfrmat     = NULL,
                             location    = country_code,
                             AsfrDatesIn = refDate_start:refDate,
                             method      = "linear",
                             output      = "single") %>%
      as.data.frame() %>% 
      rownames_to_column("age") %>% 
      pivot_longer(-c(age),
                   names_to  = "year",
                   values_to = "asfr") %>%
      mutate(year = as.numeric(year),
             age  = as.numeric(age))
    
  }
  
  # --- Birth computations --------------------------------------------------
  # Estimate annual births from ASFR and exposures, then compute age since birth.
  Bt <- left_join(expF, AsfrMat , by = join_by(year, age)) |>
    filter(year < refDate + 1) |>
    mutate(Bx = asfr * exposure) |>
    group_by(year) |>
    summarize(B = sum(Bx)) |>
    mutate(age  = refDate - year)
  
  
  
  # --- Sex ratio at birth (SRB) --------------------------------------------
  # If missing, download from WPP; otherwise, align with cohort years.
  if(!is.null(SRB)) {
    
    SRB <- tibble(SRB    = SRB,
                  cohort = refDate_start:refDate)
    
  }
  
  if(is.null(SRB)) {
    
    SRB <- download_SRB(SRB      = NULL,
                        location = country_code,
                        DatesOut = refDate_start:refDate,
                        verbose  = TRUE)
    
    SRB <- tibble(SRB    = SRB,
                  cohort = refDate_start:refDate)
    
  }
  
  # Separate births into M and F
  Bt <- Bt %>% 
    left_join(SRB, by = c("year" = "cohort")) %>% 
    mutate(Bm = B * SRB  / (1 + SRB),
           Bf = B - Bm)
  
  
  
  # --- Reverse-survival of infants -----------------------------------------
  # Use infant mortality (Lx values) to reverse-survive births and estimate
  # age 0–9 populations for females and males.
  # final part
  pop_hat_f <- mx1dt |>
    filter(country_code == !!country_code) |>
    as_tibble() |>
    select(year, age, mx = mxF) |>
    # could try to warp to PC shape here,
    # but uncertain infants. Maybe using
    # an a0 assumption it'd be doable.
    # need cohorts to structure reverse survival
    mutate(cohort = year - age - 1,
           age_int = 1,
           ax = ifelse(age == 0,
                       lt_rule_1a0_ak(M0 = mx, Sex = "f"),
                       0.5),
           qx = lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int)) |>
    filter(between(cohort, refDate_start,refDate),
           between(year, refDate_start, (refDate + 1)),
           age < 10) |>
    arrange(cohort, age) |>
    group_by(cohort) |>
    mutate(lx = lt_id_q_l(nqx = qx, radix = radix),
           dx = lx * qx,
           Lx = lx - dx * (1-ax)) |>
    filter(year == max(year)) |>
    select(age, Lx, cohort) |>
    left_join(Bt, by = join_by(age)) |>
    left_join(SRB, by = join_by(cohort)) |>
    mutate(pop_hat_f = Lx * Bf) %>%
    ungroup() %>% 
    select("age", "pop_hat_f", "cohort")
  
  # for males
  pop_hat_m <- mx1dt |>
    filter(country_code == !!country_code) |>
    as_tibble() |>
    select(year, age, mx = mxM) |>
    # could try to warp to PC shape here,
    # but uncertain infants. Maybe using
    # an a0 assumption it'd be doable.
    # need cohorts to structure reverse survival
    mutate(cohort = year - age - 1,
           age_int = 1,
           ax = ifelse(age == 0,
                       lt_rule_1a0_ak(M0 = mx, Sex = "f"),
                       0.5),
           qx = lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int)) |>
    filter(between(cohort, refDate_start,refDate),
           between(year, refDate_start, (refDate + 1)),
           age < 10) |>
    arrange(cohort, age) |>
    group_by(cohort) |>
    mutate(lx = lt_id_q_l(nqx = qx, radix = radix),
           dx = lx * qx,
           Lx = lx - dx * (1-ax)) |>
    filter(year == max(year)) |>
    select(age, Lx, cohort) |>
    left_join(Bt, by = join_by(age)) |>
    left_join(SRB, by = join_by(cohort)) |>
    mutate(pop_hat_m = Lx * Bm) %>%
    ungroup() %>% 
    select("age", "pop_hat_m", "cohort")
  
  # --- Visualization -------------------------------------------------------
  # Plots compare estimated infant population (black) to observed (red).
  # Helpful for assessing reverse survival model fit.
  females <- pop_hat_f |>
    select(age, pop =pop_hat_f) |>
    ggplot(aes(x = age, y = pop)) +
    geom_line() +
    geom_line(data = Females_single |> filter(age < 10), color = "red")
  
  males <- pop_hat_m |>
    select(age, pop = pop_hat_m) |>
    ggplot(aes(x = age, y = pop)) +
    geom_line() +
    geom_line(data = Males_single |> filter(age < 10), color = "red")
  
  return(list(pop_hat_m  = pop_hat_m,
              pop_hat_f  = pop_hat_f,
              figureF    = females,
              figureM    = males
  ))
  
}

# example: Source download function first
source("RT_UPDATE/a0_identification.R")
source("RT_UPDATE/download_asfr.R")
source("RT_UPDATE/download_nLx.R")
source("RT_UPDATE/download_SRB.R")

basepop <- basepop_single(country_code = 276, 
                          refDate      = 2000)

# results
basepop$pop_hat_m
basepop$pop_hat_f
basepop$figureM
basepop$figureF

