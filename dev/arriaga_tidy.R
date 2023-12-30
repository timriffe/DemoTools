Ages         <- seq(0, 80, by = 5)
AMales       <- smooth_age_5_arriaga(Value = pop5m_pasex, Age = Ages, OAG = TRUE)
# PAS spreadsheet result:
Atest        <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277,
                  161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189, 34729)
all(round(AMales) - Atest == 0, na.rm = TRUE)
plot(Ages, pop5m_pasex)
lines(as.integer(names(AMales)),AMales)
# Before:
smooth_age_5_arriaga <- function(Value,
                                 Age,
                                 OAG = TRUE) {
  
  # these values are not used, it's just for lengths, and to make sure we
  # end on an even 10. Technically we could even provide data in 10-year
  # age groups and it'd still not break.
  Value1     <- graduate_uniform(Value = Value, Age = Age, OAG = OAG)
  Value5     <-
    groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
  N          <- length(Value5)
  Age5       <- as.integer(names(Value5))
  
  # would need to move this up to ensure?
  # or in case of 85+ would we want to keep 80-84, 85+ as-is?
  Value10    <- groupAges(Value, Age = Age, N = 10)
  
  # what OAG is a strange digit? Then take OAG after grouping.
  if (OAG) {
    OAGvalue <- Value5[length(Value5)]
    Value10[length(Value10)] <- NA
    Value5[length(Value5)]   <- NA
  }
  
  # again get staggered vectors
  Value10LL   <- shift.vector(Value10,-2, fill = NA)
  Value10L    <- shift.vector(Value10,-1, fill = NA)
  Value10R    <- shift.vector(Value10, 1, fill = NA)
  Value10RR   <- shift.vector(Value10, 2, fill = NA)
  
  # alternating calc, with differences at tails
  vevens      <- (-Value10R + 11 * Value10 + 2 * Value10L) / 24
  # tails different
  vevens[1]   <-
    (8 * Value10[1] + 5 * Value10L[1] - Value10LL[1]) / 24
  lastind     <- which(is.na(vevens))[1]
  vevens[lastind] <-
    Value10[lastind] - (8 * Value10[lastind] + 5 * Value10R[lastind] - Value10RR[lastind]) / 24
  # odds are complement
  vodds       <- Value10 - vevens
  
  # prepare output
  interleaf  <- c(rbind(vodds, vevens))
  # produce results vector
  out        <- Value5 * NA
  n          <- min(c(length(interleaf), N))
  out[1:n]   <- interleaf[1:n]
  
  # if OA ends in 5, then we can save penultimate value too
  na.i       <- is.na(out)
  out[na.i]  <- Value5[na.i]
  if (OAG) {
    out[N] <- OAGvalue
  }
  
  out
}
library(tidyverse)
data <- tibble(Age = Ages, Pop = pop5m_pasex)

age2single <- function(Age, OAG = TRUE, AgeInt = NULL){
  maxA <- ifelse(OAG, max(Age), max(Age) + AgeInt[which.max(Age)])
  min(Age):maxA
}

library(collapse)
library(data.table)
library(tidyfast)
smooth_age_5_arriaga_tidy <- function(data, variable, OAG = TRUE){
  
  if (OAG){
    OAvalue <- data |> 
      fsubset(Age == max(Age)) |> 
      pull(!!sym(variable))
  }
 
  data <- 
    data |> 
    rename(V = !!variable)
  V  <- data$V
  A  <- data$Age
  V1 <- graduate_uniform(Value = V, 
                         Age = Age, 
                         OAG = OAG)
  A1 <- age2single(A, OAG = OAG)
    # force to single
    # reframe(
    #   V = graduate_uniform(Value = V,
    #                        Age = Age,
    #                        OAG = OAG),
    #   Age = age2single(Age))  
  data.frame(V = V1, Age = A1) |> 
  # data |> 
    # group to 5 (innocuous if 5-year data given)
    fmutate(Age = Age - Age %% 5,
            V = fifelse(Age == max(Age) & OAG, NA_real_, V)) |> 
    fgroup_by(Age) |> 
    fsummarize(V = sum(V))  |> 
    fungroup() |> 
    fmutate(Age10 = Age - Age %% 10) |> 
    fmutate(V10 = sum(V), .by = Age10) |> 
    fmutate(V10LL = flag(V10,-4),
           V10L = flag(V10,-2),
           V10R = flag(V10, 2),
           V10RR = flag(V10,4),
           vevens = dt_case_when(Age10 == min(Age10) ~           (8 * V10 + 5 * V10L - V10LL) / 24,
                          Age10 == (max(Age10) - 10) ~ V10 - (8 * V10 + 5 * V10R - V10RR) / 24,
                          TRUE ~ (-V10R + 11 * V10 + 2 * V10L) / 24),
           vodds = V10 - vevens,
           V5out = dt_case_when(
                        Age == max(Age)~ OAvalue,
                        Age %% 10 == 0 ~ vodds, 
                        Age %% 5 == 0  ~ vevens,
                        TRUE ~ V),
           V5out = fifelse(is.na(V5out), V, V5out)) |> 
    fselect(Age, Age10, V5out) |> 
    rename(!!variable := V5out)
    
  }
# ---------------- #
# test equivalency #
# ---------------- #
# for OAG divisible by 10
smooth_age_5_arriaga_tidy(tibble(Age = Ages, 
                                 Pop = pop5m_pasex), 
                          variable = "Pop", 
                          OAG = TRUE) 
  mutate(V5orig = smooth_age_5_arriaga(pop5m_pasex, Ages, TRUE))


# for OAG divisible by 5
Age75 <- seq(0,75,by=5)
v75 <- pop5m_pasex
v75[16] <- sum(pop5m_pasex[16:17])
v75 <- v75[-17]

data <- tibble(Age = Age75, 
               Pop = v75)
smooth_age_5_arriaga_tidy(data, 
                          variable = "Pop", 
                          OAG = TRUE) |> 
  mutate(V5orig = smooth_age_5_arriaga(v75, Age75, TRUE))

# ----------------------- #
# Compare execution speed #
# ----------------------- #
library(rbenchmark)
benchmark(smooth_age_5_arriaga_tidy(data, 
                                    variable = "Pop", 
                                    OAG = TRUE)) #  .522 | 0.402 if we remove reframe()
benchmark(smooth_age_5_arriaga(v75, Age75, TRUE)) # .054, 10 times faster...

# Conclusion:
# Possibly data.table could do the same a tic faster than collapse/tidyfast, 
# but only marginally do, whereas base-powered DemoTools arithmetic is 10x faster
