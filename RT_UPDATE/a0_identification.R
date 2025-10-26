# Andreev–Kingkade a0
a0_AK <- function(m0, sex = "f") {
  
  lt_rule_ak_m0_a0(M0 = m0, Sex = sex)
  
}


# Solve for a0 given L0 
solve_a0_from_L0 <- function(L0_target, 
                             l0  = 1, 
                             sex = "f") {
  
  # Basic sanity check
  if (L0_target >= l0) {
    
    stop("L0_target must be less than l0 (since deaths must occur in the interval).")
    
  }
  
  # Define function whose root gives the correct l1
  f_root <- function(l1) {
    
    m0      <- (l0 - l1) / L0_target
    a0      <- a0_AK(m0, sex)
    L0_pred <- l1 + a0 * (l0 - l1)
    return(L0_pred - L0_target)
    
  }
  
  # Numerically solve for l1
  sol <- uniroot(f_root, lower = 0, upper = l0)
  l1  <- sol$root
  
  # Compute implied m0 and a0
  m0  <- (l0 - l1) / L0_target
  a0  <- a0_AK(m0, sex)
  
  # Return everything useful
  return(list(
    a0        = a0,
    m0        = m0,
    l1        = l1,
    L0_target = L0_target,
    sex       = sex
  ))
}
