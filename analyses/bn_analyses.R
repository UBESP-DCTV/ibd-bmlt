# Bayesian Networks analyses -------------------------------------------
# If "bnlearn"l, doParallel", "caret", "glue" and "graph" are not
# installed, please run the following command before executing the
# script:
# install.packages(
#   c("bnlearn", "doParallel", "caret", "glue", "graph"),
#   dependencies = TRUE
# )
library(tidyverse)
library(bnlearn)
library(caret)
library(graph)
library(glue)
library(doParallel)

# The aim of this script is to reproduce the analyses with the Bayesian
# Network models on the simulated data.
seed <- 554477

# Load the simulated dataset in the global environment
load(here::here("data", "simulated_db.rda"))
source(here::here("analyses", "misc_functions.R"))

# 1) Model with only clinical covariates -------------------------------
db_clinical <- db %>%
  dplyr::select(contains("clinical"), y)

outcome <- "y"

# Define list with all possible algorithms of 'bnlearn'
alg_list <- c(
  'gs', 'iamb', 'fast.iamb', 'inter.iamb', 'hc', 'tabu',
  'mmhc', 'rsmax2', 'mmpc', 'si.hiton.pc', 'chow.liu', 'aracne'
)

# Learn network structures for each algorithm
cl <- makeCluster(detectCores() - 1)

bn_list <- map(
  .x = alg_list,
  ~ {

    set.seed(seed)

    bn_estimation(
      df = as.data.frame(db_clinical), algorithm = .x,
      outcome = outcome, R = 1000, k = 10, cl = cl, seed = seed
    )
  }
)

stopCluster(cl)

# Choose the network with the lowest expected classification error
expected_loss_df <- exp_loss_fun(bn_list)

bn_cv <- bn_list[[expected_loss_df]]

# Take the network with lowest classification error
loss_df <- loss_fun(bn_cv$net_cv)

bn_final <- bn_cv$net_cv[[loss_df]]

# Plot of the learned network
graphviz.plot(bn_cv$net_struc, highlight = list(nodes = 'y'))

# Performances training
bn_train_prob <- predict(
  object = bn_final$fitted,
  data   = as.data.frame(db_clinical),
  node   = outcome,
  prob   = T
)

bn_train_pred <- predict(
  object = bn_final$fitted,
  data   = as.data.frame(db_clinical),
  node   = outcome
)

cf_df <- tibble(
  pred = bn_train_pred,
  obs = db_clinical[["y"]]
)

bn_pred_prob <- attr(bn_train_prob, which = 'prob') %>%
  as_tibble(.) %>%
  slice(2) %>%
  as_vector(.)

# Test the model on bootstrap samples (1000)
n_boot <- 1000L

# Assign id to each patient
db_clinical <- db_clinical %>%
  mutate(id = seq(from = 1, to = nrow(.), by = 1))

# Creation of bootstrap samples
boot_seeds <- seq(from = 1, to = n_boot, by = 1)

boot_s <- map(
  .x = boot_seeds, ~ {

    set.seed(.x)

    bs <- sample_n(
      db_clinical, size = nrow(db_clinical), replace = TRUE
    )

    message(glue("Bootstrap sample {.x} is over"))

    bs

  }
)

# Creation of lists with performances on each bootstrap sample
bn_test_list <- NULL

bn_perf_metr_list <- NULL

for (i in 1:length(boot_s)) {

  bn_test_prob <- predict(
    object = bn_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id)),
    node = outcome,
    prob = T
  )

  bn_test_pred <- predict(
    object = bn_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id)),
    node = outcome
  )

  cf_df <- tibble(
    pred = bn_test_pred,
    obs = boot_s[[i]][["y"]]
  )

  # Get predicted probabilties
  bn_pred_prob <- attr(bn_test_prob, which = 'prob') %>%
    as_tibble(.) %>%
    slice(2) %>%
    as_vector(.)

  # Performance metrics
  bn_cf <- confusionMatrix(table(cf_df))

  bn_acc <- bn_cf$overall['Accuracy']

  bn_mcr <- 1 - bn_acc

  bn_sens <- bn_cf$byClass['Sensitivity']

  bn_spec <- bn_cf$byClass['Specificity']

  bn_ppv <- bn_cf$byClass['Pos Pred Value']

  bn_npv <- bn_cf$byClass['Neg Pred Value']

  bn_auc <- pROC::auc(
    response = boot_s[[i]][["y"]],
    predictor = bn_pred_prob
  )

  bn_somers <- 2*(bn_auc - 0.5)

  # Attach predicted probabilities and predicted class to bootstrap
  # sample
  boot_s[[i]] <- boot_s[[i]] %>%
    mutate(
      pred_prob = bn_pred_prob,
      pred_class = bn_test_pred
    )

  bn_test_list[[i]] <- list(
    'boot_pred_df' = boot_s[[i]],
    'individual_prob' = bn_test_prob,
    'predicted_values' = bn_test_pred,
    'confusion_matrix' = bn_cf
  )

  bn_perf_metr_list[[i]] <- list(
    'Accuracy' = bn_acc,
    'Overall_misc_rates' = bn_mcr,
    "Sensitivity" = bn_sens,
    "Specificity" = bn_spec,
    'Pos_pred_v' = bn_ppv,
    'Neg_pred_v' = bn_npv,
    'AUC' = bn_auc,
    'SomersD' = bn_somers
  )

  message(glue('{i} clinical bayesian network prediction finished'))
}

# Performances results averaged over boostrap samples
bn_clinical_perf <- bind_rows(bn_perf_metr_list) %>%
  summarise_all(.tbl = ., ~ mean(.))

# Individual predicted probability averaged over boostrap samples
bn_clinical_ind_pred_prob <- map_df(
  .x = bn_test_list, ~ .x[["boot_pred_df"]]
) %>%
  group_by(id) %>%
  summarise(average_pred_prob = mean(pred_prob))

# 2) Model with clinical covariates and genes --------------------------
# Learn network structure for each algorithm
cl <- makeCluster(detectCores() - 1)

bn_list <- map(
  .x = alg_list,
  ~ {

    set.seed(seed)

    bn_estimation(
      df = as.data.frame(db), algorithm = .x,
      outcome = outcome, R = 1000, k = 10, cl = cl, seed = seed
    )
  }
)

stopCluster(cl)

# Choose the network with the lowest expected classification error
expected_loss_df <- exp_loss_fun(bn_list)

bn_cv <- bn_list[[expected_loss_df]]

# Take the network with lowest classification error
loss_df <- loss_fun(bn_cv$net_cv)

bn_final <- bn_cv$net_cv[[loss_df]]

# Plot of the learned network
graphviz.plot(bn_cv$net_struc, highlight = list(nodes = 'y'))

# Performances training
bn_train_prob <- predict(
  object = bn_final$fitted,
  data   = as.data.frame(db),
  node   = outcome,
  prob   = T
)

bn_train_pred <- predict(
  object = bn_final$fitted,
  data   = as.data.frame(db),
  node   = outcome
)

cf_df <- tibble(
  pred = bn_train_pred,
  obs = db_clinical[["y"]]
)

bn_pred_prob <- attr(bn_train_prob, which = 'prob') %>%
  as_tibble(.) %>%
  slice(2) %>%
  as_vector(.)

# Test the model on bootstrap samples (1000)
n_boot <- 1000L

# Assign id to each patient
db <- db %>%
  mutate(id = seq(from = 1, to = nrow(.), by = 1))

# Creation of bootstrap samples
boot_seeds <- seq(from = 1, to = n_boot, by = 1)

boot_s <- map(
  .x = boot_seeds, ~ {

    set.seed(.x)

    bs <- sample_n(
      db, size = nrow(db), replace = TRUE
    )

    message(glue("Bootstrap sample {.x} is over"))

    bs

  }
)

# Creation of lists with performances on each bootstrap sample
bn_test_list <- NULL

bn_perf_metr_list <- NULL

for (i in 1:length(boot_s)) {

  bn_test_prob <- predict(
    object = bn_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id)),
    node = outcome,
    prob = T
  )

  bn_test_pred <- predict(
    object = bn_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id)),
    node = outcome
  )

  cf_df <- tibble(
    pred = bn_test_pred,
    obs = boot_s[[i]][["y"]]
  )

  # Get predicted probabilties
  bn_pred_prob <- attr(bn_test_prob, which = 'prob') %>%
    as_tibble(.) %>%
    slice(2) %>%
    as_vector(.)

  # Performance metrics
  bn_cf <- confusionMatrix(table(cf_df))

  bn_acc <- bn_cf$overall['Accuracy']

  bn_mcr <- 1 - bn_acc

  bn_sens <- bn_cf$byClass['Sensitivity']

  bn_spec <- bn_cf$byClass['Specificity']

  bn_ppv <- bn_cf$byClass['Pos Pred Value']

  bn_npv <- bn_cf$byClass['Neg Pred Value']

  bn_auc <- pROC::auc(
    response = boot_s[[i]][["y"]],
    predictor = bn_pred_prob
  )

  bn_somers <- 2*(bn_auc - 0.5)

  # Attach predicted probabilities and predicted class to bootstrap
  # sample
  boot_s[[i]] <- boot_s[[i]] %>%
    mutate(
      pred_prob = bn_pred_prob,
      pred_class = bn_test_pred
    )

  bn_test_list[[i]] <- list(
    'boot_pred_df' = boot_s[[i]],
    'individual_prob' = bn_test_prob,
    'predicted_values' = bn_test_pred,
    'confusion_matrix' = bn_cf
  )

  bn_perf_metr_list[[i]] <- list(
    'Accuracy' = bn_acc,
    'Overall_misc_rates' = bn_mcr,
    "Sensitivity" = bn_sens,
    "Specificity" = bn_spec,
    'Pos_pred_v' = bn_ppv,
    'Neg_pred_v' = bn_npv,
    'AUC' = bn_auc,
    'SomersD' = bn_somers
  )

  message(glue('{i} total bayesian network prediction finished'))
}

# Performances results averaged over boostrap samples
bn_total_perf <- bind_rows(bn_perf_metr_list) %>%
  summarise_all(.tbl = ., ~ mean(.))

# Individual predicted probability averaged over boostrap samples
bn_total_ind_pred_prob <- map_df(
  .x = bn_test_list, ~ .x[["boot_pred_df"]]
) %>%
  group_by(id) %>%
  summarise(average_pred_prob = mean(pred_prob))

