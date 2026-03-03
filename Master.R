#-------------------------------------------------------------------------------
#### Loading all functions and packages needed to simulate
#-------------------------------------------------------------------------------
packages <- c(
  "bootnet",
  "qgraph",
  "mantar",
  "ggplot2",
  "future",
  "future.apply",
  "glasso",
  "fastmap",
  "progressr"
)
# install.packages(packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)

# Simulated annealing functions
source("functions_cost_function.R") # contains estimation methods of theta, loglikelihood, (E)BIC etc.
source("cost_function.R") # contains the cost function which estimates BIC or EBIC based on binary vector
source("functions_simulated_annealing.R") # Temperature, neighbourhood step, acceptance rate, bitstrings
source("simulated_annealing.R") # core algorithm

#-------------------------------------------------------------------------------
#### Setting up condition parameters
#-------------------------------------------------------------------------------
sample_sizes <- c(150, 400, 600, 1000, 2500, 5000)
repetitions <- 500 # set to 500 in real deal, 1 to 10 to try out

#-------------------------------------------------------------------------------
#### Setting up the true network / data
#-------------------------------------------------------------------------------
true_net <- read.csv("data/MAGNA_ptsd.csv")
true_net <- true_net[1:length(true_net)-1]
true_net <- as.matrix(true_net)

#-------------------------------------------------------------------------------
#### Setting up data generation scheme from Isvoranu & Epskamp (2023), extended to correlation matrices for direct use with methods
#-------------------------------------------------------------------------------
dataGenerator <- function(
    trueNet,
    sampleSize,
    type = c("normal"),
    nLevels = 4
){
  type <- match.arg(type)
  nNode <- ncol(trueNet)
  
  if (type == "normal"){
    # Generator function:
    gen <- ggmGenerator()
    
    # Generate data:
    Data <- gen(sampleSize,trueNet)
    
  }
  sample_correlation <- cor(Data)
  #return(data = Data)
  return(sample_correlation)
}

#-------------------------------------------------------------------------------
#### Setting up a replication structure
#-------------------------------------------------------------------------------
replication <- function(
    sample_correlation,
    sample_size
){
  
  # BIC Estimation Methods
  BIC_SA <- simulated_annealing(S = sample_correlation, n = sample_size, EBIC = FALSE)
  
  t0 <- proc.time()
  BIC_nonconvex <- mantar::regularization_net(mat = sample_correlation, ns = sample_size, likelihood = "mat_based")
  BIC_nonconvex_runtime <- proc.time() - t0
  
  t0 <- proc.time()
  BIC_ggmModSelect <- qgraph::ggmModSelect(S = sample_correlation, n = sample_size)
  BIC_ggmModSelect_runtime <- proc.time() - t0
  
  # EBIC Estimation Methods
  EBIC_SA <- simulated_annealing(S = sample_correlation, n = sample_size, EBIC = TRUE)
  
  t0 <- proc.time()
  EBIC_nonconvex <- mantar::regularization_net(mat = sample_correlation, ns = sample_size, likelihood = "mat_based", extended = TRUE)
  EBIC_nonconvex_runtime <- proc.time() - t0
  
  t0 <- proc.time()
  EBIC_ggmModSelect <- qgraph::ggmModSelect(S = sample_correlation, n = sample_size, gamma = 0.5)
  EBIC_ggmModSelect_runtime <- proc.time() - t0
  
  return(list(
    BIC_SA = list(
      graph = BIC_SA$best_vector,
      runtime = BIC_SA$runtime[1],
      convergence = BIC_SA$convergence,
      cache_uses = BIC_SA$cache_uses,
      initial_temperature = BIC_SA$temperature,
      first_accept_rate = BIC_SA$first_accept,
      start = BIC_SA$start_vector,
      start_ic = BIC_SA$start_value,
      best_ic = BIC_SA$best_value
    ),
    BIC_nonconvex = list(
      graph = BIC_nonconvex$pcor,
      runtime = BIC_nonconvex_runtime[1]
    ),
    BIC_ggmModSelect = list(
      graph = BIC_ggmModSelect$graph,
      runtime = BIC_ggmModSelect_runtime[1]
    ),
    EBIC_SA = list(
      graph = EBIC_SA$best_vector,
      runtime = EBIC_SA$runtime[1],
      convergence = EBIC_SA$convergence,
      cache_uses = EBIC_SA$cache_uses,
      initial_temperature = EBIC_SA$temperature,
      first_accept_rate = EBIC_SA$first_accept,
      start = EBIC_SA$start_vector,
      start_ic = EBIC_SA$start_value,
      best_ic = EBIC_SA$best_value
    ),
    EBIC_nonconvex = list(
      graph = EBIC_nonconvex$pcor,
      runtime = EBIC_nonconvex_runtime[1]
    ),
    EBIC_ggmModSelect = list(
      graph = EBIC_ggmModSelect$graph,
      runtime = EBIC_ggmModSelect_runtime[1]
    )
  ))
}

#-------------------------------------------------------------------------------
#### Parallel simulation
#-------------------------------------------------------------------------------

plan(multisession) # initial: change to multicore for Linux server, now to multisession for safety with c-compiled packages

#handlers(global = TRUE)
#handlers("txtprogressbar") # progress bar

set.seed(123)

for (i_n in 1:length(sample_sizes)) {
  
  n <- sample_sizes[i_n]
  design <- expand.grid(
    n = n,
    r = 1:repetitions
  )
  
  results <- with_progress({
    p <- progressor(along = 1:nrow(design))
    
    future_lapply(
      1:nrow(design),
      function(i) {
        
        n <- design$n[i]
        r <- design$r[i]
        
        # Generate data once per replication
        sample_correlation <- dataGenerator(
          trueNet = true_net,
          sampleSize = n,
          type = "normal"
        )
        
        # Apply all methods to same data
        est <- replication(
          sample_correlation = sample_correlation,
          sample_size = n
        )
        
        p()
        
        # Collect in list
        list(
          n = n,
          sample_correlation = sample_correlation,
          rep = r,
          estimates = est
        )
      },
      future.seed = TRUE
    )
  })
  
  saveRDS(results, file = paste0("output/condition_n_", sample_sizes[i_n], ".RDS"))
  
}

