#' Check and Load the Latest Version of Available `wpp` Package
#'
#' This function checks the PPgp GitHub organization for the latest available
#' `wpp` package (e.g., `wpp2024`, `wpp2022`, etc.), compares it with any
#' installed versions, and optionally updates to and loads the most recent version.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Retrieves all `wpp`-like repositories from the \href{https://github.com/PPgp}{PPgp GitHub}.
#'   \item Identifies the most recent year-based `wpp` package (e.g., `wpp2024`).
#'   \item Checks which `wpp` version, if any, is currently installed.
#'   \item Compares installed and available versions.
#'   \item Prompts the user to update to the latest version, if applicable.
#'   \item Loads the most recent version into the R session.
#' }
#'
#' @return
#' Returns the name of the `wpp` package loaded (e.g., `"wpp2024"`).
#' Also loads the package into the current R session.
#'
#' @note
#' Requires internet access and the following packages: `jsonlite`, `httr`, `remotes`.
#'
#' @examples
#' \dontrun{
#' check_and_load_latest_wpp()
#' }
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @importFrom utils installed.packages compareVersion packageVersion
#' @importFrom remotes install_github
#' @export

# function
check_and_load_latest_wpp <- function() {
  
  # ---------------------------------------------------------------------------- #
  # Step 1: Get all wpp-like repositories from GitHub
  repos     <- jsonlite::fromJSON("https://api.github.com/users/PPgp/repos?per_page=100")
  wpp_repos <- repos$name[grepl("^wpp\\d{4}$", repos$name)]
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 2: Extract years and determine the latest
  wpp_years <- as.integer(sub("wpp", "", wpp_repos))
  if(length(wpp_years) == 0) {
    
    stop("No wpp packages found on GitHub.")
    
  }
  latest_year <- max(wpp_years)
  latest_pkg  <- paste0("wpp", latest_year)
  
  message("Latest available wpp package: ", latest_pkg)
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 3: Check which wpp package (if any) is currently installed
  installed_pkgs <- utils::installed.packages()[, "Package"]
  installed_wpp  <- installed_pkgs[grepl("^wpp\\d{4}$", installed_pkgs)]
  
  if(length(installed_wpp) == 0) {
    
    message("No wpp package currently installed.")
    
    installed_version <- NA_character_
    
  } else {
    
    installed_wpp_years <- as.integer(sub("wpp", "", installed_wpp))
    current_pkg         <- installed_wpp[which.max(installed_wpp_years)]
    installed_version   <- as.character(packageVersion(current_pkg))
    
    message("Installed wpp package: ", current_pkg, " (v", installed_version, ")")
  }
  
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 4: Get version of the latest package from GitHub
  desc_url <- paste0("https://raw.githubusercontent.com/PPgp/", latest_pkg, "/main/DESCRIPTION")
  
  desc <- tryCatch(
    httr::content(httr::GET(desc_url), "text", encoding = "UTF-8"),
    error = function(e)
      stop("Could not read DESCRIPTION for ", latest_pkg)
  )
  
  latest_version <- sub("^Version:\\s+", 
                        "",
                        grep("^Version:", strsplit(desc, "\n")[[1]], value = TRUE))
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 5: Compare versions
  needs_update <- is.na(installed_version) || 
    utils::compareVersion(installed_version, latest_version) < 0
  
  if(needs_update) {
    
    message("Newer version available: ", latest_pkg, " (v", latest_version, ")")
    
    ans <- readline(prompt = "Do you want to update to the latest version? [y/n]: ")
    
    if(tolower(ans) == "y") {
      
      message("Installing ", latest_pkg, "...")
      
      remotes::install_github(paste0("PPgp/", latest_pkg), upgrade = "always")
      
      message("Updated to ", latest_pkg, " (v", latest_version, ")")
      
    } else {
      
      message("Update skipped. Using currently installed version.")
      
    }
    
  } else {
    
    message("Installed version is already up-to-date.")
    
  }
  
  # ---------------------------------------------------------------------------- #
  # Step 6: Load the most recent version available (installed or just updated)
  final_pkg <- if(needs_update && tolower(ans) == "y") { 
    
    latest_pkg
    
  } else {
    
    current_pkg
    
  }
  
  library(final_pkg, character.only = TRUE)
  
  # show the version that was returned
  final_pkg
  # ---------------------------------------------------------------------------- #
  
}

# ---------------------------------------------------------------------------- #
# Example use
check_and_load_latest_wpp()
data(pop1dt)
# ---------------------------------------------------------------------------- #