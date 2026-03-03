simulated_annealing <- function(
    S,
    n,
    start_vector = NULL,
    EBIC = TRUE,
    # Algorithm parameters
    max_iterations  = 100000L,
    T_start = NULL,
    T_end = 0.001,
    cooling = 0.95,
    max_flip_chance = 0, # only one flip at a time if set to 0
    k =  1, # k defines the multiplier of the iterations per cooling level, k = 1 => roughly 63% of edges visited from starting neighbourhood, k = 2 => 86%, k = 3 => 95%
    initial_accept = 0.85 # inital acceptance rate for temperature initialization
) {
  t0 <- proc.time()
  # start state initialization
  if (is.null(start_vector)) {
    # Generate theta, turn into vector
    theta_start <- Matrix::forceSymmetric(glasso::glasso(S, rho = 0)$wi)
    theta_vec <- theta_start[lower.tri(theta_start)]
    # initialize empty start vector
    start_vector <- numeric(length(theta_vec))
    # split into top and bottom 50%
    n_vals <- length(theta_vec)
    ranks <- rank(theta_vec, ties.method = "first")
    top_index <- ranks > theta_vec / 2
    bottom_index <- ranks <= theta_vec / 2
    # apply probabilites for initial values based on split as predefined
    start_vector[top_index] <- as.numeric(runif(sum(top_index)) < 0.75)
    start_vector[bottom_index] <- as.numeric(runif(sum(bottom_index)) < 0.25)
  }
  # temperature initialization
  if (is.null(T_start)) {
    T_start = define_T_start(start_vector = start_vector, S = S, n = n, EBIC = EBIC, uphill_accept = initial_accept, max_flip_chance=max_flip_chance)
  }
  # further variable initialization
  nodes <- ncol(S) # number of nodes, needed for some calculations
  edges <- (nodes*(nodes-1))/2 # for undirected graphs its nodes over 2, needed to define Neighbourhood size of a given State, S_i for L_k formula
  L_k <- k*edges # idea of L_k formula from Delahaye (2018)
  counter <- 0
  T_current <- T_start
  current_vector <- start_vector
  current_value <- cost_function(start_vector, S, n, nodes, EBIC)
  best_vector <- current_vector
  start_value <- best_value <- current_value
  no_improve_counter <- 0
  convergence <- TRUE
  T_endstate <- FALSE
  repeats <- 0 # Track how often cache was used for
  max_no_improve = 3*edges # stop if no improvements to best after 3*edges, covering roughly 95% of possible edges in the final temperature state where only better solutions will be accepted
  cache <- fastmap::fastmap()
  
  # Do first run outside of loop to sample the acceptance rate of the first run
  first_acceptance_sum <- 0 # initialize a sum to capture first L_k acceptance rate to check if temperature initialization went correctly
  first_acceptance_counter <- 0
  for (i in 1:L_k) {
    counter <- counter + 1
    candidate <- vector_flip(current_vector, T_start, T_current, max_flip_chance)
    
    # Check with cache before calling cost_function
    candidate_string <- bitstring(candidate)
    if (cache$has(candidate_string)) {
      candidate_value <- cache$get(candidate_string)
      repeats <- repeats + 1
    } else {
      candidate_value <- cost_function(candidate, S, n, nodes, EBIC)
      cache$set(candidate_string, candidate_value)
    }
    
    
    delta_E <- candidate_value - current_value
    # calculate acceptance probability
    accept_prob <- acceptance_probability(delta_E, T_current)
    
    
    # acceptance loop
    if (accept_prob == 1) {
      # always update if better
      current_vector <- candidate
      current_value  <- candidate_value
    } else {
      # add to sum if uphill and thus not automatic 1
      first_acceptance_sum <- first_acceptance_sum + accept_prob
      first_acceptance_counter <- first_acceptance_counter + 1
      if (runif(1) < accept_prob){
        # update even if worse (hill-climb)
        current_vector <- candidate
        current_value  <- candidate_value
      }
    }
    # check if best has changed
    if (current_value < best_value) {
      best_vector <- current_vector
      best_value <- current_value
      # found a new best, therefore reset
      no_improve_counter <- 0
    } else {
      # new current is worse than best
      no_improve_counter <- no_improve_counter + 1
    }
  }
  empirical_first_acceptance_rate <- first_acceptance_sum / first_acceptance_counter
  
  # geometric cooling
  if(!T_endstate) {
    T_current <- T_current * cooling
  }
  
  # perform annealing
  while (counter < max_iterations) {
    for (i in 1:L_k) {
      counter <- counter + 1
      candidate <- vector_flip(current_vector, T_start, T_current, max_flip_chance)
      
      # Check with cache before calling cost_function
      candidate_string <- bitstring(candidate)
      if (cache$has(candidate_string)) {
        candidate_value <- cache$get(candidate_string)
        repeats <- repeats + 1
      } else {
        candidate_value <- cost_function(candidate, S, n, nodes, EBIC)
        cache$set(candidate_string, candidate_value)
      }
      
      delta_E <- candidate_value - current_value
      # calculate acceptance probability
      accept_prob <- acceptance_probability(delta_E, T_current)
      
      # acceptance loop
      if (accept_prob == 1) {
        # always update if better
        current_vector <- candidate
        current_value  <- candidate_value
      } else {
        if (runif(1) < accept_prob){
          # update even if worse (hill-climb)
          current_vector <- candidate
          current_value  <- candidate_value
        }
      }
      # check if best has changed
      if (current_value < best_value) {
        best_vector <- current_vector
        best_value <- current_value
        # found a new best, therefore reset
        no_improve_counter <- 0
      } else {
        # new current is worse than best
        no_improve_counter <- no_improve_counter + 1
      }
    }
    
    # geometric cooling
    if(!T_endstate) {
      T_current <- T_current * cooling
    }
    
    
    if (T_current <= T_end) {
      T_endstate <- TRUE
    }
    # Convergence - ending before max_iterations
    if (T_endstate && no_improve_counter >= max_no_improve) {
      break
    }
  }
  
  convergence <- (counter < max_iterations)
  if (convergence) {
    convergence <- counter
  }
  runtime <- proc.time() - t0 # measure runtime
  
  return(list(
    best_vector = best_vector,
    best_value = best_value,
    start_vector = start_vector,
    start_value = start_value,
    convergence = convergence,
    cache_uses = repeats,
    EBIC = EBIC,
    temperature = T_start,
    first_accept = empirical_first_acceptance_rate,
    runtime = runtime
  ))
}


#### Testing -------------------------------------------------------------------
if (sys.nframe() == 0) {
  # results_1 <- simulated_annealing(start_vector = null_vec, S = sample_cov, n = sample_size, EBIC = TRUE) # works fine
  # print(results_1)
  
  results_1 <- simulated_annealing(S = sample_cov, n = sample_size, EBIC = FALSE) # works fine
  print(results_1)
  
  results_2 <- simulated_annealing(start_vector = start_state, S = sample_cov, n = sample_size, EBIC = FALSE) # works fine
  print(results_2)
  start_vector <- results_2$best_vector
  results_3 <- simulated_annealing(start_vector = start_vector, S = sample_cov, n = sample_size, EBIC = TRUE) # works fine
  print(results_3)
  start_vector <- results_3$best_vector
  
  # results_3 <- simulated_annealing(start_vector = algo_run$best_sol$edge_encoding, S = sample_cov, n = sample_size, EBIC = TRUE) # works fine
  # print(results_3)
  
  # Feel free to modify here! Like trying with EBIC=F
  #results <- simulated_annealing(start_vector = start_state, S = sample_cov, n = sample_size, EBIC = FALSE) # works fine
  #print(results)
}



