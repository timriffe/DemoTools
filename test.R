library(tidyverse)
library(readxl)
library(janitor)

xl_file <- "~/Downloads/Li_2010_WPP_BY1950_MIG_Canada.xlsx"


# INPUTS ------------------------------------------------------------------
# example data

# Survival ratios
sr_m <- read_xlsx(xl_file, sheet = 2)
sr_f <- read_xlsx(xl_file, sheet = 3)

# Age-specific-fertility rate
asfr <- read_xlsx(xl_file, sheet = 4)

# Population
pop_m <- read_xlsx(xl_file, sheet = 5)
pop_f <- read_xlsx(xl_file, sheet = 6)

pop_female <- pop_f
pop_male <- pop_m
survival_ratio_male <- sr_m
survival_ratio_female <- sr_f

# Migration
mig_m <- read_xlsx(xl_file, sheet = 7)
mig_f <- read_xlsx(xl_file, sheet = 8)
mig_mL <- read_xlsx(xl_file, sheet = 9)
mig_mU <- read_xlsx(xl_file, sheet = 10)
mig_fL <- read_xlsx(xl_file, sheet = 11)
mig_fU <- read_xlsx(xl_file, sheet = 12)

# pull out sex ratio at birth
srb <- asfr %>% filter(Age == "SRB") %>% select(-Age)

# only keep asgr
asfr <- asfr %>% filter(Age != "SRB")
srb <- srb %>% pivot_longer(everything(), names_to = "year", values_to = "srb")

age_min <- 0
age_max <- 100
age_interval <- 5
year_min <- 1950
year_max <- 2050
year_interval <- 5


# DEFINE AGES AND YEARS ---------------------------------------------------

ages <- seq(age_min, age_max, by = age_interval)
years <- seq(year_min, year_max, by = year_interval)

age_labels <- c(paste(ages[1:(length(ages)-1)], ages[2:length(ages)]-1, sep = "-"), paste0(age_max, "+"))
year_labels <- paste(years[1:(length(years)-1)], years[2:length(years)], sep = "-")

ages_fertility <- seq(15, 45, by = age_interval)
age_labels_fertility <- c(paste(ages_fertility[1:(length(ages_fertility)-1)],
                                ages_fertility[2:length(ages_fertility)]-1, sep = "-"),
                          paste0(ages_fertility[length(ages_fertility)], "-49"))


# GET IN TIDY FORM --------------------------------------------------------


# assume inputs are of the excel form but we need to make it long and put it all together

names(pop_m) <- tolower(names(pop_m))
names(pop_f) <- tolower(names(pop_f))
names(sr_m) <- tolower(names(sr_m))
names(sr_f) <- tolower(names(sr_f))
names(sr_f) <- tolower(names(sr_f))
names(asfr) <- tolower(names(asfr))

# Female/Male population merged in long format
df_pop <-
  pop_m %>%
  mutate(age = ages) %>%
  pivot_longer(-age, names_to = "year", values_to = "population") %>%
  mutate(sex = "m") %>%
  bind_rows(pop_f %>%
              mutate(age = ages) %>%
              pivot_longer(-age, names_to = "year", values_to = "population") %>%
              mutate(sex = "f")) %>%
  mutate(year = as.numeric(year))

# Survival ratios for Male/Female merged in long format
df_sr <-
  sr_m %>%
  mutate(age = ages) %>%
  pivot_longer(-age, names_to = "year_original", values_to = "sr") %>%
  mutate(sex = "m") %>%
  bind_rows(sr_f %>%
              mutate(age = ages) %>%
              pivot_longer(-age, names_to = "year_original", values_to = "sr") %>%
              mutate(sex = "f")) %>%
  rowwise() %>%
  mutate(year = as.numeric(years[which(year_labels==year_original)])) %>%
  select(age, year, sr, sex)

# Age-specific fertility rate long format
df_asfr <-
  asfr %>%
  mutate(age = ages_fertility) %>%
  pivot_longer(-age, names_to = "year_original", values_to = "asfr") %>%
  rowwise() %>%
  mutate(year = as.numeric(years[which(year_labels==year_original)])) %>%
  select(age, year, asfr)

# sex ratio at birth
df_srb <-
  srb %>%
  ## pivot_longer(everything(), names_to = "year_original", values_to = "srb") %>%
  rename(year_original = year) %>%
  rowwise() %>%
  mutate(year = as.numeric(years[which(year_labels==year_original)])) %>%
  select(year, srb)


# Calculate survivors -----------------------------------------------------

# Calculate survivors and deceased
test <-
  df_pop %>%
  left_join(df_sr) %>%
  mutate(cohort = year - age) %>%
  arrange(sex, year, age) %>%
  mutate(implied_survivors = population * lead(sr)) %>%
  group_by(sex, cohort) %>%
  arrange(sex, cohort, age) %>%
  mutate(diff = population - lag(implied_survivors)) %>%
  arrange(year, sex, age)

# Replicate CohortMigFlow -------------------------------------------------


# Try copying vba code

pop_m_mat <- as.matrix(pop_m %>% select(-age))
pop_f_mat <- as.matrix(pop_f %>% select(-age))
sr_m_mat <- as.matrix(sr_m %>% select(-age))
sr_f_mat <- as.matrix(sr_f %>% select(-age))
asfr_mat <- as.matrix(asfr %>% select(-age))
srb_mat <- as.matrix(srb)

nages <- length(ages)
nyears <- length(years)

net_mig_m <- matrix(NA, nrow = nages, ncol = nyears)
net_mig_f <- matrix(NA, nrow = nages, ncol = nyears)

# 1. net_mig is pop minus people that survived
for (i in 2:(nages - 1)) {
  for (j in 2:(nyears)) {
    previous_year <- j - 1
    previous_age <- i - 1

    # Net migration is pop minus the people that survived from the previous age/cohort
    net_mig_m[i, j] <- pop_m_mat[i, j] -  pop_m_mat[previous_age, previous_year] * sr_m_mat[i, previous_year]
    net_mig_f[i, j] <- pop_f_mat[i, j] -  pop_f_mat[previous_age, previous_year] * sr_f_mat[i, previous_year]
  }
}

# Last age group
for (j in 2:nyears) {
  previous_year <- j - 1
  previous_ageg <- nages - 1

  # For last age group, net migration is:
  # pop for that age group in year j, minus the people from the previous age group
  # the survived
  net_mig_m[nages, j] <- pop_m_mat[nages, j] -  (pop_m_mat[nages, previous_year] + pop_m_mat[previous_ageg, previous_year]) * sr_m_mat[nages, previous_year]
  net_mig_f[nages, j] <- pop_f_mat[nages, j] -  (pop_f_mat[nages, previous_year] + pop_f_mat[previous_ageg, previous_year]) * sr_f_mat[nages, previous_year]
}

## Up until this point, net_mig does not have any value
## for the first age group because we need to know
## the survival rate of the previous age group, which
## is inexistent.

# births
all_births <- rep(NA, length(years))
fertility_index <- which(ages %in% ages_fertility)

for (j in 2:nyears) {
  births_this_year <- 0
  for (i in fertility_index) {
    print(i)

    # Sum female pop from previous year and this year
    f_pop <- pop_f_mat[i, j - 1] + pop_f_mat[i, j]

    ## Get the asfr of the previous year for the current
    ## age group. normalize_index just matches the index
    ## from fertility index to the actual index of
    ## the asfr_mat (so the index 4 in fertility_index
    ## is just 1 in asfr_mat, etc..)
    normalize_index <- i - (min(fertility_index) - 1)
    asfr_previousyear <- asfr_mat[normalize_index, j - 1]

    ## Births that occurred this year j for all age groups
    these_births <- age_interval * (0.5 * (f_pop) * asfr_previousyear) / 1000

    ## Accumulate all the births for the age groups
    ## fir this year
    births_this_year <- births_this_year + these_births
    print(births_this_year)
  }
  # All births that occurred this year. This
  # is a vector same as length as number of years
  # with births
  all_births[j] <- births_this_year
}

# With all_births already calculated, separate between
# female/male births with the sex ratio at birth
srb_vec <- as.numeric(srb_mat[, 2])
births_m <- all_births[2:length(all_births)] * (srb_vec / (1 + srb_vec))
births_f <- all_births[2:length(all_births)] * (1 / (1 + srb_vec))

# 2. So far, net_mig is pop minus people that survived minus births
for (j in 2:nyears) {
  previous_year <- j - 1

  ## pop for the first age group minus the births for the
  ## previous year * by the survival ratio for the previous year.
  ## In other words, tell me all the births that survived and
  ## subtract from the population of the first age group.
  net_mig_m[1, j] <- pop_m_mat[1, j] - births_m[previous_year] * sr_m_mat[1, previous_year]
  net_mig_f[1, j] <- pop_f_mat[1, j] - births_f[previous_year] * sr_f_mat[1, previous_year]

  # Why all of this only for the first age group? Because we didn't
  # calculate the net_mig for the this first age group before. The only
  # thing we can't calculate is the net migration for the first year
  # because we don't info on the *previous year*. Other than that,
  # the age groups x years matrix should be filled entirely.
}

mig_upper_m <- matrix(NA, nrow = nages, ncol = nyears)
mig_lower_m <- matrix(NA, nrow = nages, ncol = nyears)
mig_rectangle_m <- matrix(NA, nrow = nages, ncol = nyears)

mig_upper_f <- matrix(NA, nrow = nages, ncol = nyears)
mig_lower_f <- matrix(NA, nrow = nages, ncol = nyears)
mig_rectangle_f <- matrix(NA, nrow = nages, ncol = nyears)

# 3. Estimate the lower/upper bounds for the net migration
for (i in 2:(nages - 1)) {
  for (j in 2:(nyears)) {
    previous_year <- j - 1
    previous_ageg <- i - 1

    # Upper bound is net mig / 2 times the survival ratio of last year ^ 0.5
    mig_upper_m[i, j] <- net_mig_m[i, j] / (2 * sr_m_mat[i, previous_year]^0.5)
    # Lower bound is net migration minus the upper bound for the
    # previous age group
    mig_lower_m[previous_ageg, j] <- net_mig_m[i, j] - mig_upper_m[i, j]

    mig_upper_f[i, j] <- net_mig_f[i, j] / (2 * sr_f_mat[i, previous_year]^0.5)
    mig_lower_f[previous_ageg, j] <- net_mig_f[i, j] - mig_upper_f[i, j]

  }
}

# Estimate upper bounds for the first age group. Why
# no lower bound for the first age group? because we have
# no previous age group. See above.
for (j in 2:(nyears)) {
  previous_year <- j - 1
  mig_upper_m[1, j] <- net_mig_m[1, j] / (sr_m_mat[1, previous_year]^0.5)
  mig_upper_f[1, j] <- net_mig_f[1, j] / (sr_f_mat[1, previous_year]^0.5)
}

# last age group
for (j in 2:(nyears)) {
  previous_ageg <- nages - 1
  mig_lower_m[previous_ageg, j] <- mig_upper_m[previous_ageg, j]
  mig_lower_f[previous_ageg, j] <- mig_upper_f[previous_ageg, j]

  mig_upper_m[nages, j] <- net_mig_m[nages, j] * 0.5
  mig_upper_f[nages, j] <- net_mig_f[nages, j] * 0.5
  mig_lower_m[nages, j] <- net_mig_m[nages, j] * 0.5
  mig_lower_f[nages, j] <- net_mig_f[nages, j] * 0.5
}

# Combine both upper/lower bound into a single rectangle
for (i in 1:nages) {
  for (j in 1:nyears) {
    mig_rectangle_m[i, j] <- mig_upper_m[i, j] + mig_lower_m[i, j]
    mig_rectangle_f[i, j] <- mig_upper_f[i, j] + mig_lower_f[i, j]
  }
}


# compare -----------------------------------------------------------------

round(as.matrix(mig_mL %>% select(-Age)) - mig_lower_m[,2:nyears], 1)
round(as.matrix(mig_mU %>% select(-Age)) - mig_upper_m[,2:nyears], 2)
round(as.matrix(mig_m %>% select(-Age)) - mig_rectangle_m[,2:nyears], 2)

round(as.matrix(mig_fL %>% select(-Age)) - mig_lower_f[,2:nyears], 2)
round(as.matrix(mig_fU %>% select(-Age)) - mig_upper_f[,2:nyears], 2)
round(as.matrix(mig_f %>% select(-Age)) - mig_rectangle_f[,2:nyears], 2)

mig_m %>%
  pivot_longer(-Age, names_to = "year_original", values_to = "mig") %>%
  rowwise() %>%
  mutate(year = as.numeric(years[which(year_labels==year_original)]),
         age = as.numeric(ages[which(age_labels==Age)])) %>%
  select(age, year, mig) %>%
  mutate(cohort = year - age) %>%
  ggplot(aes(age, mig, color = factor(cohort))) + geom_line()


plot(mig_rectangle_m[2,2:20], type = "o")
lines(as.matrix(mig_m[,-1])[2,2:20] , col = 2)


plot(mig_rectangle_m[2:20,3], type = "o")
lines(as.matrix(mig_m[,-1])[2:20,2] , col = 2)
