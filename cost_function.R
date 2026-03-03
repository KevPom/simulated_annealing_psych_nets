cost_function <- function(given_vector, sample_cov, sample_size, nodes, EBIC=T) {
  
  # transform edge encoding into adjacency matrix
  adjacency <- vector_to_adj(given_vector, nodes)
  
  # reestimate network with given adjaceny matrix
  theta <- tryCatch(
    reestimate(sample_cov, adjacency),
    error = function(e) return(NULL)
  )
  if (is.null(theta)) {
    return(cost = Inf)
  }
  
  # use matrix to calculate (E)-BIC
  if(EBIC) {
    cost <- gaussian_EBIC(
      theta = theta,
      S = sample_cov,
      n = sample_size,
      gamma = 0.5
    )
  } else {
    cost <- gaussian_EBIC(
      theta = theta,
      S = sample_cov,
      n = sample_size,
      gamma = 0
    )
  }
  
  # return value of cost function
  return(cost)
}

#### Testing -------------------------------------------------------------------
if (sys.nframe() == 0) {
  cost_function(start_state, sample_cov, sample_size, ncol(sample_cov)) # works fine
}
