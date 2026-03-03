vector_flip <- function(
    chosen_vector,
    T_start = NULL,
    T_current = NULL,
    p_max = 0.05) {
  # define flip chance based on temperature
  if (!is.null(T_start) && !is.null(T_current)) {
    p <- min(p_max, T_current/T_start)
  } else {
    p <- p_max
  }
  
  # define flips
  n <- length(chosen_vector)
  flip <- as.logical(rbinom(n, 1, p))
  
  # flip at least one
  if (!any(flip)) {
    flip[sample.int(n,1)] <- TRUE
  }
  
  chosen_vector[flip] <- !chosen_vector[flip]
  
  return(chosen_vector)
}

# Based roughly on Ben-Ameur, Walid (2004) to estimate T0 recursively, should be a bit more accurate than using the typical closed form
define_T_start <- function(
    start_vector,
    S,
    n,
    EBIC = TRUE,
    neighbour_proposals = 500,
    uphill_accept = 0.85,
    p = 1,
    max_flip_chance = 0.05
) {
  # Initialize some neighbours to estimate differences of uphill values
  nodes <- ncol(S)
  
  # Collect uphill transitions
  delta_vec <- numeric(0)
  
  current <- start_vector
  Ei <- cost_function(current, S, n, nodes, EBIC)
  
  for (i in 1:neighbour_proposals) {
    neighbour <- vector_flip(current, p_max = 0)  # single flip
    Ej <- cost_function(neighbour, S, n, nodes, EBIC)
    
    if (Ej > Ei) {
      delta_vec <- c(delta_vec, Ej - Ei)
    }
  }
  
  T_n <- -mean(delta_vec) / log(uphill_accept) # starting with closed form
  
  # Safety check for minimum of 1 uphill value
  if (length(delta_vec)==0) {
    stop("No uphill moves found.")
  }
  
  #### Recursive Iteration for T ----
  repeat {
    # chi_hat formula
    x_hat_T <- mean(exp(-delta_vec / T_n))
    
    # break condition
    if (abs(x_hat_T - uphill_accept) < 0.001) break # some arbitrary epsilon
    
    # determine next T_value if not yet satisfactory
    T_n <- T_n * (log(x_hat_T) / log(uphill_accept))^(1/p)
  }
  return(T_n)
}

acceptance_probability <- function(delta_E, T) {
  if (delta_E < 0) { 
    return(1) # immediately accept if better
  } else {
    return(exp(-delta_E / T)) # else Manhattan / Simulated Annealing T
  }
}

bitstring <- function(v) {
  paste0(v, collapse = "")
}


#### Testing -------------------------------------------------------------------
if (sys.nframe() == 0) {
  define_T_start(start_vector = start_state, S = sample_cov, n = sample_size) # works fine
}


