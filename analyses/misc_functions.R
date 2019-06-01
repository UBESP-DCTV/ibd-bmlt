# Miscellaneous functions ----------------------------------------------
# This script contains several functions that will be used to reproduce
# the analysis

# Bayesian network estimation ------------------------------------------
bn_estimation <- function(df, algorithm, outcome, R, k, cl, seed){

  # Measure arc strenghts with bootstrap
  set.seed(seed)
  arcs <- boot.strength(
    data = df,
    algorithm = algorithm,
    R = R,
    cpdag = F,
    cluster = cl
  )

  # Learn structure of the network
  set.seed(seed)
  bn <- cextend(averaged.network(arcs))

  #  Cross-validate the network
  set.seed(seed)
  bn_cv <- bn.cv(
    data = df,
    bn = bn,
    fit = 'bayes',
    loss = 'pred',
    loss.args = list(target = outcome),
    method = 'k-fold',
    k = k,
    runs = 1,
    cluster = cl
  )

  # Get network loss
  loss_v <- NULL

  for (i in 1:length(bn_cv)) {
    loss_v[i] <- bn_cv[[i]]$loss
  }

  message(glue::glue('{algorithm} finished'))

  value_list <- list(
    'arcs' = arcs,
    'net_struc' = bn,
    'net_cv' = bn_cv,
    'algorithm' = algorithm,
    'expected_loss' = mean(loss_v)
  )
}

# Choose the network with the lower MCR --------------------------------
exp_loss_fun <- function(bn_list_object){
  exp_loss_list <- NULL

  for (i in 1:length(bn_list_object)) {
    exp_loss_list[[i]] <- list(
      'index' = i,
      'exp_loss' = bn_list_object[[i]]$expected_loss
    )
  }

  df <- bind_rows(exp_loss_list)

  which.min(df$exp_loss)
}

# Choose the cv network with the lower MCR -----------------------------
loss_fun <- function(cv_object){

  loss_list <- NULL

  for (i in 1:length(cv_object)) {
    loss_list[[i]] <- list(
      'index' = i,
      'loss' = cv_object[[i]]$loss
    )
  }

  df <- bind_rows(loss_list)

  which.min(df$loss)

}
