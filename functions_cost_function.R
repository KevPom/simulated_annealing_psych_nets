vector_to_adj <- function(maskVector, p) {
  A <- matrix(0, p, p)
  A[lower.tri(A)] <- maskVector
  A <- A + t(A)
  diag(A) <- 1
  return(A)
}

retrieve_zeros <- function(matrix){
  p <- nrow(matrix)
  
  #which gives one value index for every entry of matrix
  #divide value by p and return the remainder
  row <- which(matrix == 0)%%p
  
  #values without remainder are in row p
  row[row == 0] <- p
  
  #calculate column by division and using next bigger integer value
  col <- ceiling(which(matrix == 0)/p)
  
  zero <- cbind(row, col)
  return(zero)
}

# second function uses first function to reestimate the covariance matrix based on the results

reestimate <- function(covMat, adjMat){
  #find indexes for zero entries
  zero_indices <- retrieve_zeros(matrix = adjMat)
  if (length(zero_indices) == 0) zero_indices <- NULL
  
  #refit network with glasso and edges forced to zero
  inverse <- Matrix::forceSymmetric(glasso::glasso(covMat, rho = 0, zero = zero_indices, penalize.diagonal = FALSE)$wi) # penalize.diagonal not needed with 1 on diagonal of adjacency but doesnt hurt as safety
  inverse <- as.matrix(inverse)
  return(inverse)
}


#### Goal parameter(s) ----
# LogLik
LL_mat <- function(S, theta, n) {
  const1 <- - n/2 * ncol(theta) * log(2 * pi)
  ll <- n/2 * (log(det(theta)) - sum(diag(S %*% theta)))
  return(ll)
}

# EBIC - call with gamma = 0 and you have BIC too!
gaussian_EBIC <- function (indi = FALSE, theta, gamma = 0.5,
                           S = NULL, n = NULL, data = NULL,
                           E = NULL, countDiagonal = FALSE)
{
  
  if(is.null(data) & (is.null(S) | is.null(n))){
    stop("Either the raw data or the sample correlation matrix and the number of observations are required.")
  }
  
  if(!indi){
    if (is.null(S)){S <- cor(data)}
    if (is.null(n)){n <- nrow(data)}
    ll <- LL_mat(S = S, theta = theta, n = n)
  }
  else{
    stop("Individual log-likelihoods have not been implemented yet")
    if(is.null(data)) {stop("Raw data have to be available for individual log-likelihood computation.")}
  }
  
  if (is.null(E)) {
    E <- sum(theta[lower.tri(theta, diag = countDiagonal)] != 0)
  }
  p <- nrow(theta)
  EBIC <- -2 * ll + E * log(n) + 4 * E * gamma * log(p)
  return(EBIC)
}


#### Testing -------------------------------------------------------------------
if (sys.nframe() == 0) {
  vector_to_adj(start_state, ncol(sample_cov)) # works fine
  
  test_zeroes <- retrieve_zeros(adjacency)
}

