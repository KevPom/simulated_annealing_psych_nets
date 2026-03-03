# Set Working Directory to analysis folder (analysis folder with this script on same level as data and output folder, all folders within simulation scripts folder)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(purrr)
library(tidyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

#-------------------------------------------------------------------------------
#### Creating true net and top 50% / 25% true net for recovery metrics
#-------------------------------------------------------------------------------
# load true net, convert to matrix
true_net <- read.csv("../data/MAGNA_ptsd.csv")
true_net <- true_net[1:length(true_net)-1]
true_net <- as.matrix(true_net)

# convert to edge-vector
true_edges <- true_net[lower.tri(true_net)]

# top 50%
true_top50 <- which(true_edges > median(abs(true_edges[true_edges!=0])))
true_top50_edges <- rep(0, length(true_edges))
true_top50_edges[true_top50] <- 1

# top 25%
true_top25 <- which(true_edges > quantile(abs(true_edges[true_edges!=0]), 0.75))
true_top25_edges <- rep(0, length(true_edges))
true_top25_edges[true_top25] <- 1

# top 10%
true_top10 <- which(true_edges > quantile(abs(true_edges[true_edges!=0]), 0.90))
true_top10_edges <- rep(0, length(true_edges))
true_top10_edges[true_top10] <- 1

#-------------------------------------------------------------------------------
#### Loading in the lists of the conditions, turning them into datasets
#-------------------------------------------------------------------------------
datapath <- "../output"

files <- list.files(datapath, full.names = T)

for (f in files) {
  obj_name <- tools::file_path_sans_ext(basename(f))
  assign(obj_name, readRDS(f))
}

#-------------------------------------------------------------------------------
#### Combining datasets using tidyverse packages
#-------------------------------------------------------------------------------
# put into shared list
results <- list(
  `150` = condition_n_150,
  `400` = condition_n_400,
  `600` = condition_n_600,
  `1000` = condition_n_1000,
  `2500` = condition_n_2500,
  `5000` = condition_n_5000
)

# turn large list into dataframe that reduces to just graph and runtime
# for within-algorithm analysis rely on different df later!
comparison_df <- map_dfr(results, ~ {
  map_dfr(.x, function(x) {
    map_dfr(names(x$estimates), function(method) {
      tibble(
        n = x$n,
        sample_correlation = list(x$sample_correlation),
        rep = x$rep,
        method = method,
        graph = list(x$estimates[[method]]$graph),
        runtime = x$estimates[[method]]$runtime
      )
    })
  })
})

# saveRDS(comparison_df, file = "comparison.RDS") # save base version of comparison_df

# Use adjacency vector to estimate final networks of SA-methods
# Using same reestimate/vector_to_adj/retrieve_zeros functions from algorithm for parity
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
reestimate <- function(covMat, adjMat){
  #find indexes for zero entries
  zero_indices <- retrieve_zeros(matrix = adjMat)
  if (length(zero_indices) == 0) zero_indices <- NULL
  
  #refit network with glasso and edges forced to zero
  inverse <- Matrix::forceSymmetric(glasso::glasso(covMat, rho = 0, zero = zero_indices, penalize.diagonal = FALSE)$wi) # penalize.diagonal not needed with 1 on diagonal of adjacency but doesnt hurt as safety
  inverse <- as.matrix(inverse)
  return(inverse)
}
#calculate partial from precision matrix
precision_to_partial <- function(inverse) {
  D <- diag(1 / sqrt(diag(inverse)))
  P <- -D %*% inverse %*% D
  diag(P) <- 0
  return(P)
}

comparison_df <- comparison_df %>%
  mutate(
    graph = if_else(
      method %in% c("BIC_SA", "EBIC_SA"),
      map2(
        graph,
        sample_correlation,
        ~ {
          adj  <- vector_to_adj(.x, 17)
          prec <- reestimate(.y, adj)
          precision_to_partial(prec)
        }
      ),
      graph
    )
  )


# turn graphs into lower triangular vectors for easier comparison metrics
comparison_df <- comparison_df %>%
  mutate(lower_tri = map(graph, ~ {
    if (is.matrix(.x)) {
      .x[lower.tri(.x)]
    } else {
      .x
    }
  }))

# check that everything got correctly converted into a vector of length 136 (17*(17-1)/2)
lengths(comparison_df$graph) %>% table()

#-------------------------------------------------------------------------------
#### Calculate the metrics of interest (means and deviations)
#-------------------------------------------------------------------------------
# Sensitivity
# Top 50% / 25% / 10% Sensitivity
sensitivity <- function(pred, truth) {
  TP <- sum(pred == 1 & truth == 1)
  FN <- sum(pred == 0 & truth == 1)
  return(TP / (TP + FN))
}
# Specificity
specificity <- function(pred, truth) {
  TN <- sum(pred == 0 & truth == 0)
  FP <- sum(pred == 1 & truth == 0)
  return(TN / (TN + FP))
}

# binary version of true_edges
true_edges_bin <- as.integer(true_edges != 0)

# add entries for specificity, sensitivity
comparison_df <- comparison_df %>%
  mutate(
    lower_tri_bin = map(lower_tri, ~ as.integer(.x != 0)),  # ensure 0/1
    sens_true = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_edges_bin)),
    spec_true = map_dbl(lower_tri_bin, ~ specificity(.x, true_edges_bin)),
    sens_top10 = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_top10_edges)),
    sens_top25 = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_top25_edges)),
    sens_top50 = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_top50_edges))
  )

## Strength Correlation
strength <- function(mat) {
  rowSums(abs(mat)) - abs(diag(mat))
}
# generate true net strength
true_strength <- strength(true_net)

# add entries for strength correlation and z-scaled strength correlation
comparison_df <- comparison_df %>%
  mutate(
    # store estimated strength per graph
    est_strength = map(graph, strength),
    # compute correlation with true_strength
    strength_cor = map_dbl(est_strength, ~ cor(.x, true_strength)),
    # Fisher-Z-transform it
    z_strength_cor = atanh(strength_cor)
  )
# report this in the summary df listing rather than the prior 
mean_r <- tanh(mean(comparison_df$z_strength_cor))
sd_r <- tanh(sd(comparison_df$z_strength_cor))

# plot one to check for linearity - normality should be fine given sample size
est_vec <- comparison_df$est_strength[[1]]

plot(est_vec, true_strength,
     xlab = "Estimated Strength",
     ylab = "True Strength")
abline(lm(true_strength ~ est_vec), col = "red")
cor_val <- cor(est_vec, true_strength)
text(0.5, 0.5, paste("r =", round(cor_val, 2)))


# Convergence Rate is a 100%; no missing data - no need to record here
# instead its numerical values (meaning iterations) will be used later-on for caching-percentage analysis


# Estimation Time + prior metrics in summarized form
summary_df <- comparison_df %>%
  group_by(method, n) %>%
  summarise(
    runtime_mean = mean(runtime, na.rm = TRUE),
    runtime_sd   = sd(runtime, na.rm = TRUE),
    
    strength_cor_mean = tanh(mean(z_strength_cor, na.rm = TRUE)),
    strength_cor_sd = tanh(sd(z_strength_cor, na.rm = TRUE)),
    
    # old
    # strength_cor_mean = mean(strength_cor, na.rm = TRUE),
    # strength_cor_sd   = sd(strength_cor, na.rm = TRUE),
    
    sens_true_mean = mean(sens_true, na.rm = TRUE) * 100,
    sens_true_sd   = sd(sens_true, na.rm = TRUE) * 100,
    
    spec_true_mean = mean(spec_true, na.rm = TRUE) * 100,
    spec_true_sd   = sd(spec_true, na.rm = TRUE) * 100,
    
    sens_top10_mean = mean(sens_top10, na.rm = TRUE) * 100,
    sens_top10_sd   = sd(sens_top10, na.rm = TRUE) * 100,
    
    sens_top25_mean = mean(sens_top25, na.rm = TRUE) * 100,
    sens_top25_sd   = sd(sens_top25, na.rm = TRUE) * 100,
    
    sens_top50_mean = mean(sens_top50, na.rm = TRUE) * 100,
    sens_top50_sd   = sd(sens_top50, na.rm = TRUE) * 100,
    
    .groups = "drop"
  ) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))

summary_long <- summary_df %>%
  pivot_longer(
    cols = -c(method, n),
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  ) %>%
  arrange(metric, n, method)


#-------------------------------------------------------------------------------
#### Detailed SA dataframe
#-------------------------------------------------------------------------------
desired_methods <- c("EBIC_SA", "BIC_SA")

SA_df <- map_dfr(results, ~ {
  map_dfr(.x, function(x) {
    map_dfr(intersect(names(x$estimates), desired_methods), function(method) {
      tibble(
        n = x$n,
        sample_correlation = list(x$sample_correlation),
        rep = x$rep,
        method = method,
        graph = list(x$estimates[[method]]$graph),
        runtime = x$estimates[[method]]$runtime,
        convergence = x$estimates[[method]]$convergence,
        cache_uses = x$estimates[[method]]$cache_uses,
        initial_temperatures = x$estimates[[method]]$initial_temperatures,
        first_accept_rate = x$estimates[[method]]$first_accept_rate,
        start = list(x$estimates[[method]]$start),
        start_ic = x$estimates[[method]]$start_ic,
        best_ic = x$estimates[[method]]$best_ic
      )
    })
  })
})

# saveRDS(SA_df, file = "Simulated_Annealing.RDS") # save base version of SA_df

#-------------------------------------------------------------------------------
#### Additional analysis of SA algorithm going past the preregistration
#-------------------------------------------------------------------------------
# Percentage cached
mean((1-(SA_df$convergence-SA_df$cache_uses)/SA_df$convergence))
sd((1-(SA_df$convergence-SA_df$cache_uses)/SA_df$convergence))

# Percentage of total state space covered
edge_amount <- (17*16)/2
total_state_space <- 2**(edge_amount)
mean(SA_df$convergence/total_state_space * 100)

# Accuracy of temperature estimation method
mean(SA_df$first_accept_rate)
sd(SA_df$first_accept_rate)

# Influence of start state accuracy on end state accuracy
# Start and end matrices
SA_df <- SA_df %>%
  mutate(
    graph = if_else(
      method %in% c("BIC_SA", "EBIC_SA"),
      map2(
        graph,
        sample_correlation,
        ~ {
          adj  <- vector_to_adj(.x, 17)
          prec <- reestimate(.y, adj)
          precision_to_partial(prec)
        }
      ),
      graph
    )
  )
SA_df <- SA_df %>%
  mutate(
    start = if_else(
      method %in% c("BIC_SA", "EBIC_SA"),
      map2(
        start,
        sample_correlation,
        ~ {
          adj  <- vector_to_adj(.x, 17)
          prec <- reestimate(.y, adj)
          precision_to_partial(prec)
        }
      ),
      start
    )
  )
SA_df <- SA_df %>%
  mutate(lower_tri = map(graph, ~ {
    if (is.matrix(.x)) {
      .x[lower.tri(.x)]
    } else {
      .x
    }
  }))
SA_df <- SA_df %>%
  mutate(start_lower_tri = map(start, ~ {
    if (is.matrix(.x)) {
      .x[lower.tri(.x)]
    } else {
      .x
    }
  }))

SA_df <- SA_df %>%
  mutate(
    lower_tri_bin = map(lower_tri, ~ as.integer(.x != 0)),  # ensure 0/1
    sens_true = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_edges_bin)),
    spec_true = map_dbl(lower_tri_bin, ~ specificity(.x, true_edges_bin)),
    sens_top10 = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_top10_edges)),
    sens_top25 = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_top25_edges)),
    sens_top50 = map_dbl(lower_tri_bin, ~ sensitivity(.x, true_top50_edges)),
    start_lower_tri_bin = map(start_lower_tri, ~ as.integer(.x != 0)),  # ensure 0/1
    start_sens_true = map_dbl(start_lower_tri_bin, ~ sensitivity(.x, true_edges_bin)),
    start_spec_true = map_dbl(start_lower_tri_bin, ~ specificity(.x, true_edges_bin)),
    start_sens_top10 = map_dbl(start_lower_tri_bin, ~ sensitivity(.x, true_top10_edges)),
    start_sens_top25 = map_dbl(start_lower_tri_bin, ~ sensitivity(.x, true_top25_edges)),
    start_sens_top50 = map_dbl(start_lower_tri_bin, ~ sensitivity(.x, true_top50_edges))
  )
## Correlation of start and end_state metrics
# Metrics to check
metrics <- c("sens_true", "spec_true", "ic")#, "sens_top10", "sens_top25", "sens_top50")

# Compute delta metrics (improvement)
SA_df <- SA_df %>%
  mutate(
    delta_sens_true = sens_true - start_sens_true,
    delta_spec_true = spec_true - start_spec_true,
    delta_sens_top10 = sens_top10 - start_sens_top10,
    delta_sens_top25 = sens_top25 - start_sens_top25,
    delta_sens_top50 = sens_top50 - start_sens_top50,
    delta_ic = best_ic - start_ic
  )

# Function to compute both correlations
SA_correlations_df <- SA_df %>%
  group_by(method, n) %>%
  group_modify(~ {
    
    map_dfr(metrics, function(metric) {
      start_metric <- paste0("start_", metric)
      delta_metric <- paste0("delta_", metric)
      
      tibble(
        metric = metric,
        pearson_start_end = cor(.x[[start_metric]], .x[[metric]], method = "pearson"),
        spearman_start_end = cor(.x[[start_metric]], .x[[metric]], method = "spearman"),
        pearson_start_delta = cor(.x[[start_metric]], .x[[delta_metric]], method = "pearson"),
        spearman_start_delta = cor(.x[[start_metric]], .x[[delta_metric]], method = "spearman")
      )
    })
    
  }) %>%
  ungroup()

# There is no NAs, yet there are NAs in correlations_df -> NAs are caused due to constant values!

# SA_NA_check <- SA_df %>%
#   group_by(method, n) %>%
#   summarise(across(everything(), ~ sum(is.na(.))), .groups = "drop")

# Yes, NAs are caused by lack of variance in entries

# SA_variance_table <- SA_df %>%
#   group_by(method, n) %>%
#   summarise(across(where(is.numeric), ~ var(., na.rm = TRUE)), .groups = "drop")


# Quick visualization: Start vs End and Start vs Delta
for (metric in metrics) {
  
  start_metric <- paste0("start_", metric)
  delta_metric <- paste0("delta_", metric)
  
  p1 <- ggplot(SA_df, aes_string(x = start_metric, y = metric)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE) +
    facet_grid(method ~ n) +
    labs(title = paste0("Start vs End: ", metric),
         x = "Start", y = "End")
  
  p2 <- ggplot(SA_df, aes_string(x = start_metric, y = delta_metric)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    facet_grid(method ~ n) +
    labs(title = paste0("Start vs Improvement: ", metric),
         x = "Start", y = "Delta")
  print(p1)
  print(p2)
}


# Start vs Best IC
p1 <- ggplot(SA_df, aes(x = start_ic, y = best_ic)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE) +
  facet_grid(method ~ n, scales = "free") +
  labs(title = "Start vs Best IC",
       x = "Start IC",
       y = "Best IC")

print(p1)

p2 <- ggplot(SA_df, aes(x = start_ic, y = delta_ic)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  facet_grid(method ~ n, scales = "free") +
  labs(title = "Start vs Improvement IC",
       x = "Start IC",
       y = "Delta IC")

print(p2)

# A look at the SD / Mean of the values between start and end, only to look at since we now have correlation already
# mostly relic from before correlation
SA_summary_df <- SA_df %>%
  group_by(method, n) %>%
  summarise(
    # Analysis of mean of sensitivity
    mean_start_sens = mean(start_sens_true),
    mean_end_sens = mean(sens_true),
    # Analysis of standard deviations of sensitivity
    sd_start_sens = sd(start_sens_true),
    sd_end_sens = sd(sens_true),
    # Analysis of mean of specificity
    mean_start_spec = mean(start_spec_true),
    mean_end_spec = mean(spec_true),
    # Analysis of standard deviations of specificity
    sd_start_spec = sd(start_spec_true),
    sd_end_spec = sd(spec_true),
    # Mean ICs between start and end
    mean_start_ic = mean(start_ic),
    mean_end_ic = mean(best_ic),
    # Now finally look at IC deviation
    sd_start_ic = sd(start_ic),
    sd_end_ic = sd(best_ic)
  )

#-------------------------------------------------------------------------------
#### Generate initial Table using rempsyc package
#-------------------------------------------------------------------------------
summary_df_ordered <- summary_df %>%
  arrange(n, method)

table_2_base <- rempsyc::nice_table(summary_df_ordered)

flextable::save_as_docx(table_2_base, path = "table_2_template.docx")

#-------------------------------------------------------------------------------
#### ggplot2 for main part
#-------------------------------------------------------------------------------
# Filter to two main metrics
summary_long_filtered <- summary_long %>%
  filter(!metric %in% c("runtime"))

summary_long_filtered <- summary_long_filtered %>%
  mutate(metric = recode(metric,
                         "sens_top10" = "Sensitivity (Top 10%)",
                         "sens_top25" = "Sensitivity (Top 25%)",
                         "sens_top50" = "Sensitivity (Top 50%)",
                         "sens_true" = "Sensitivity",
                         "spec_true" = "Specificity",
                         "strength_cor" = "Strength Correlation"))


# APA-ish theme
theme_apa <- jtools::theme_apa()

# Line plot comparing methods
ggplot(summary_long_filtered,
       aes(x = n, y = mean, color = method, shape = method, group = method)) +
  geom_line(linewidth = 1, linetype = "dashed") +
  geom_point(size = 3) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Sample Size (n)",
    y = "Mean Value",
    color = "Method",
    shape = "Method"
  ) +
  theme_apa + scale_color_hue(h = c(180, 300))



#-------------------------------------------------------------------------------
#### Linear mixed model to compare methods <- ONLY for initial experimentation / exploration, not relevant to analysis!
#-------------------------------------------------------------------------------
# define names of metrics that are analyzed (all but strength_cor)
metrics <- c(
  "sens_true",
  "spec_true",
  "sens_top10",
  "sens_top25",
  "sens_top50"
)

# structure of single model, since we want to compare methods we cluster by unique datasets on intercept
# metric ~ method * n + (1 | n:rep)

# Turn n and method into factor for easier comparison
comparison_df$n <- factor(comparison_df$n)
comparison_df$method <- factor(comparison_df$method)

# add a dataset factor for simpler representation of n and rep
comparison_df <- comparison_df %>%
  mutate(dataset = interaction(n, rep, drop = TRUE))

# apply models for our metrics of interest
models <- lapply(metrics, function(metric) {
  lmer(as.formula(paste(metric, "~ method * n + (1 | dataset)")),
       data = comparison_df)
})

# name models after metrics
models <- setNames(models, metrics)

# Use emmeans to look at model comparisons
emm_options(pbkrtest.limit = 18000)

# Sensitivity
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "BIC_SA")
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "EBIC_SA")

# Specificity
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "BIC_SA")
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "EBIC_SA")

# Sensitivity - Top 10%
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "BIC_SA")
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "EBIC_SA")

# Sensitivity - Top 25%
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "BIC_SA")
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "EBIC_SA")

# Sensitivity - Top 50%
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "BIC_SA")
emmeans(models$sens_true, trt.vs.ctrl ~ method, ref = "EBIC_SA")


