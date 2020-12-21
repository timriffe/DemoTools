library(tidyverse)
library(readxl)
library(janitor)

xl_file <- "~/Downloads/canada_migflow.xlsx"

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
mig_m <- read_xlsx(xl_file, sheet = 7) %>% mutate(across(where(is.numeric), round, 5))
mig_f <- read_xlsx(xl_file, sheet = 8) %>% mutate(across(where(is.numeric), round, 5))
mig_mL <- read_xlsx(xl_file, sheet = 9) %>% mutate(across(where(is.numeric), round, 5))
mig_mU <- read_xlsx(xl_file, sheet = 10) %>% mutate(across(where(is.numeric), round, 5))
mig_fL <- read_xlsx(xl_file, sheet = 11) %>% mutate(across(where(is.numeric), round, 5))
mig_fU <- read_xlsx(xl_file, sheet = 12) %>% mutate(across(where(is.numeric), round, 5))

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

pop_m_mat <- as.matrix(pop_m %>% select(-age))
pop_f_mat <- as.matrix(pop_f %>% select(-age))
sr_m_mat <- as.matrix(sr_m %>% select(-age))
sr_f_mat <- as.matrix(sr_f %>% select(-age))
asfr_mat <- as.matrix(asfr %>% select(-age))
srb_vec <- as.numeric(as.matrix(srb[, 2]))

nages <- length(ages)
nyears <- length(years)

# Replicate MigFlow -------------------------------------------------

mig_res <-
  mig_resid_stock(
    pop_m_mat = pop_m_mat,
    pop_f_mat = pop_f_mat,
    sr_m_mat = sr_m_mat,
    sr_f_mat = sr_f_mat,
    asfr_mat = asfr_mat,
    srb_vec = srb_vec,
    ages = ages,
    ages_fertility = ages_fertility
  )

net_mig_m_t <- mig_res$mig_m[1:4, 1:5]
sr_m_mat <- sr_m_mat[1:4, 1:5]

n <- nrow(net_mig_m_t)
p <- ncol(net_mig_m_t)

net_mig_m_t[1, ] <- 2 * net_mig_m_t[1, ]

double_pop <- (2 * net_mig_m_t[2:n, ])
mig_sr <- net_mig_m_t[-n, ] * sr_m_mat[2:n, ]
net_mig_m_t[2:n, ] <-  double_pop - mig_sr

corrected <- round(as.matrix(mig_m %>% select(-Age)), 5)[1:4, 1:5]
res <- round(net_mig_m_t, 5)

corrected - res


## net_mig_f_t <- net_mig_f
## net_mig_f_t[1, ] <- 2 * net_mig_f_t[1, 2]
## net_mig_f_t[next_ageg, j] <- 2 * net_mig_f_t[i, j] - net_mig_f_t[i, j] * sr_f_mat[next_ageg, j]


## for (j in 2:nyears) {
##   print("j is:")
##   print(j)
##   for (i in 2:nages) {
##     print(i)
##     next_ageg <- i
##     net_mig_m_t[next_ageg, j] <- 2 * net_mig_m_t[i, j] - net_mig_m_t[i, j] * sr_m_mat[next_ageg, j]
##     net_mig_f_t[next_ageg, j] <- 2 * net_mig_f_t[i, j] - net_mig_f_t[i, j] * sr_f_mat[next_ageg, j]
##   }
## }


## round(as.matrix(mig_f %>% select(-Age)), 2) - round(net_mig_f_t[, 2:nyears], 5)
