# remotes::install_github("lorenzwalthert/gitsum")
if (system("whoami",intern=TRUE) == "tim"){
  library(gitsum)
  library(tidyverse)
  library(readr)
  library(lubridate)
  
  get_versions <- function(hash){
    D <- readLines(paste0("https://raw.githubusercontent.com/timriffe/DemoTools/",hash,"/DESCRIPTION") )
    D[grepl(D,pattern = "Version: ")]  %>%  
      gsub(pattern = "Version: ", replacement = "") 
  }
  
  update_lookup <- function(){
    DESC_changes <- 
      parse_log_detailed() %>% 
      unnest_log() %>% 
      dplyr::filter(changed_file == "DESCRIPTION") %>% 
      group_by(lubridate::as_date(date)) %>% 
      slice(n()) %>% 
      ungroup() %>% 
      mutate(date = as_date(date))
        #date = paste(year(date),month(date),day(date),sep="-"))
    # DESC_changes %>% 
    #   mutate(get_versions(hash))
    vers <- list()
    for (i in 1:nrow(DESC_changes)){
      # Sys.sleep(1)
      vers[[i]] <- try(get_versions(DESC_changes$hash[i]))
    }
    closeAllConnections()
    errors <- lapply(vers,class) %>% unlist() 
    errors <- errors == "try-error"
    vers[errors]<-NA
    DESC_changes <- cbind(DESC_changes,version=tibble(version=unlist(vers)))
    
    version_lookup <- 
      DESC_changes %>% 
      filter(!is.na(version)) %>% 
      select(date,version, hash) 
    
    write_csv(version_lookup,"version_lookup.csv")
  
   }
}

get_DemoTools_versions <- function(){
  readr::read_csv("https://raw.githubusercontent.com/timriffe/DemoTools/master/version_lookup.csv")
}

install_DemoTools_version <- function(version = NULL, date = NULL, hash = NULL){
  
  
  
  if (is.null(version) & is.null(date) & is.null(hash)){
    cat("version identifiers include (one of):
        1. Version entries in the DESCRIPTION file
        2. Date (yyyy-mm-dd)
        3. the hash key of the commit\n\n")
    cat("no identifier given, you can view options with 
    get_DemoTools_version()\n\n")
    cat("to install the current head, try:
        remotes::install_github('timriffe/DemoTools'")
    }
  
  versions <- get_DemoTools_versions() %>% 
    mutate(date = lubridate::as_date(date))
  
  out <- NULL
  if (!is.null(hash)){
    if (hash %in% version$hash){
      remotes::install_github("timriffe/DemoTools", ref = hash)
    } else {
      stop("hash not found, try using one from the (incomplete) list returned by:\n get_DemoTools_versions()")
    }
    out <- 1
  }
  
  # try date
  if (is.null(out)){
    if (!is.null(date)){
      if (is.character(date)){
        date = lubridate::ymd(date)
      }
      if (is.na(date)){
        cat("date didn't parse. It should be yyyy-mm-dd format if given as character\n")
        stop()
      }
      dateK <- date
        dateL <- 
          versions %>% 
          mutate(dist = abs(dateK - date)) %>% 
          dplyr::filter(dist == min(dist)) %>% 
          dplyr::pull(date) %>% 
          '['(1)
        hash <- 
          versions %>% 
          dplyr::filter(date == dateL) %>% 
          dplyr::pull(hash)
        remotes::install_github("timriffe/DemoTools", ref = hash)
        out <- 1
        if (date != dateL){
          cat("date not in the (incomplete) set returned by get_DemoTools_versions()
              using the closest date in that subset instead:", as.character(dateL),"\n")
        }
      }
  }
  
  # try version
  if (is.null(out)){
    if (!is.null(version)){
      .version <- version
      hash <- 
        versions %>% 
        dplyr::filter(version == .version) %>% 
        dplyr::pull(hash)
      if (length(hash) == 1){
        remotes::install_github("timriffe/DemoTools", ref = hash)  
        out <- 1
      } else {
        cat("version not in the set returned by get_DemoTools_versions()
            we didn't try approximating. Have a look at the version snapshots available
            in the lookup table.\n")
      }
    }
  }
  # catch-all
  if (is.null(out)){
    cat("Looks like no installation attempted. For the most recent version try:
        remotes::install_github('timriffe/DemoTools')
        otherwise easiest thing to roll back is to pick something 
        out from the lookup table.")
  }
    
  
}
  
  
