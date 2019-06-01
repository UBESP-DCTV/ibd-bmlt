# Naive-Bayes analyses -------------------------------------------------
# If "bnlearn" and "doParallel" are not installed, please run the
# following command before executing the script:
# install.packages(
#   c("bnlearn", "doParallel", "caret", "glue"), dependencies = TRUE
# )
library(tidyverse)
library(bnlearn)
library(caret)
library(glue)
library(doParallel)

# The aim of this script is to reproduce the analyses with the Naive
# Bayes models on the simulated data.
seed <- 554477

# Load the simulated dataset in the global environment
load(here::here("data", "simulated_db.rda"))
source(here::here("analyses", "misc_functions.R"))

# 1) Model with only clinical covariates -------------------------------
db_clinical <- db %>%
  dplyr::select(contains("clinical"), y)

# Store the covariates and the outcome names into different objects
covs <- db_clinical %>%
  dplyr::select(contains("clinical")) %>%
  names()

outcome <- "y"

# Learn network structure
nb_struc <- naive.bayes(
  x = as.data.frame(db_clinical), training = outcome, explanatory = covs
)

# Training of the network
cl <- makeCluster(detectCores() - 1)

set.seed(seed)
nb_cv <- bn.cv(
  data = as.data.frame(db_clinical), bn = nb_struc, fit = 'bayes',
  loss = 'pred', loss.args = list(target = 'y'), method = 'k-fold',
  k = 10, cluster = cl
)

stopCluster(cl)

# Takes the one with lowest classification error
loss_df <- loss_fun(nb_cv)

nb_final <- nb_cv[[loss_df]]

# Performances training
nb_train_prob <- predict(
  object = nb_final$fitted,
  data   = as.data.frame(db_clinical),
  prob   = T
)

nb_train_pred <- predict(
  object = nb_final$fitted,
  data   = as.data.frame(db_clinical)
)

cf_df <- tibble(
  pred = nb_train_pred,
  obs = db_clinical[["y"]]
)

nb_pred_prob <- attr(nb_train_prob, which = 'prob') %>%
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
nb_test_list <- NULL

nb_perf_metr_list <- NULL

for (i in 1:length(boot_s)) {

  nb_test_prob <- predict(
    object = nb_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id)),
    prob = T
  )

  nb_test_pred <- predict(
    object = nb_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id))
  )

  cf_df <- tibble(
    pred = nb_test_pred,
    obs = boot_s[[i]][["y"]]
  )

  # Get predicted probabilties
  nb_pred_prob <- attr(nb_test_prob, which = 'prob') %>%
    as_tibble(.) %>%
    slice(2) %>%
    as_vector(.)

  # Performance metrics
  nb_cf <- confusionMatrix(table(cf_df))

  nb_acc <- nb_cf$overall['Accuracy']

  nb_mcr <- 1 - nb_acc

  nb_sens <- nb_cf$byClass['Sensitivity']

  nb_spec <- nb_cf$byClass['Specificity']

  nb_ppv <- nb_cf$byClass['Pos Pred Value']

  nb_npv <- nb_cf$byClass['Neg Pred Value']

  nb_auc <- pROC::auc(
    response = boot_s[[i]][["y"]],
    predictor = nb_pred_prob
  )

  nb_somers <- 2*(nb_auc - 0.5)

  # Attach predicted probabilities and predicted class to bootstrap
  # sample
  boot_s[[i]] <- boot_s[[i]] %>%
    mutate(
      pred_prob = nb_pred_prob,
      pred_class = nb_test_pred
    )

  nb_test_list[[i]] <- list(
    'boot_pred_df' = boot_s[[i]],
    'individual_prob' = nb_test_prob,
    'predicted_values' = nb_test_pred,
    'confusion_matrix' = nb_cf
  )

  nb_perf_metr_list[[i]] <- list(
    'Accuracy' = nb_acc,
    'Overall_misc_rates' = nb_mcr,
    "Sensitivity" = nb_sens,
    "Specificity" = nb_spec,
    'Pos_pred_v' = nb_ppv,
    'Neg_pred_v' = nb_npv,
    'AUC' = nb_auc,
    'SomersD' = nb_somers
  )

  message(glue('{i} clinical naive bayes prediction finished'))
}

# Performances results averaged over boostrap samples
nb_clinical_perf <- bind_rows(nb_perf_metr_list) %>%
  summarise_all(.tbl = ., ~ mean(.))

# Individual predicted probability averaged over boostrap samples
nb_clinical_ind_pred_prob <- map_df(
  .x = nb_test_list, ~ .x[["boot_pred_df"]]
) %>%
  group_by(id) %>%
  summarise(average_pred_prob = mean(pred_prob))

# 2) Model with clinical covariates and genes --------------------------
# Store the covariates and the outcome names into different objects
covs <- db %>%
  dplyr::select(-y) %>%
  names()

outcome <- "y"

# Learn network structure
nb_struc <- naive.bayes(
  x = as.data.frame(db), training = outcome, explanatory = covs
)

# Training of the network
cl <- makeCluster(detectCores() - 1)

set.seed(seed)
nb_cv <- bn.cv(
  data = as.data.frame(db), bn = nb_struc, fit = 'bayes',
  loss = 'pred', loss.args = list(target = 'y'), method = 'k-fold',
  k = 10, cluster = cl
)

stopCluster(cl)

# Takes the one with lowest classification error
loss_df <- loss_fun(nb_cv)

nb_final <- nb_cv[[loss_df]]

# Performances training
nb_train_prob <- predict(
  object = nb_final$fitted,
  data   = as.data.frame(db),
  prob   = T
)

nb_train_pred <- predict(
  object = nb_final$fitted,
  data   = as.data.frame(db)
)

cf_df <- tibble(
  pred = nb_train_pred,
  obs = db[["y"]]
)

nb_pred_prob <- attr(nb_train_prob, which = 'prob') %>%
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
nb_test_list <- NULL

nb_perf_metr_list <- NULL

for (i in 1:length(boot_s)) {

  nb_test_prob <- predict(
    object = nb_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id)),
    prob = T
  )

  nb_test_pred <- predict(
    object = nb_final$fitted,
    data = as.data.frame(boot_s[[i]] %>% dplyr::select(-id))
  )

  cf_df <- tibble(
    pred = nb_test_pred,
    obs = boot_s[[i]][["y"]]
  )

  # Get predicted probabilties
  nb_pred_prob <- attr(nb_test_prob, which = 'prob') %>%
    as_tibble(.) %>%
    slice(2) %>%
    as_vector(.)

  # Performance metrics
  nb_cf <- confusionMatrix(table(cf_df))

  nb_acc <- nb_cf$overall['Accuracy']

  nb_mcr <- 1 - nb_acc

  nb_sens <- nb_cf$byClass['Sensitivity']

  nb_spec <- nb_cf$byClass['Specificity']

  nb_ppv <- nb_cf$byClass['Pos Pred Value']

  nb_npv <- nb_cf$byClass['Neg Pred Value']

  nb_auc <- pROC::auc(
    response = boot_s[[i]][["y"]],
    predictor = nb_pred_prob
  )

  nb_somers <- 2*(nb_auc - 0.5)

  # Attach predicted probabilities and predicted class to bootstrap
  # sample
  boot_s[[i]] <- boot_s[[i]] %>%
    mutate(
      pred_prob = nb_pred_prob,
      pred_class = nb_test_pred
    )

  nb_test_list[[i]] <- list(
    'boot_pred_df' = boot_s[[i]],
    'individual_prob' = nb_test_prob,
    'predicted_values' = nb_test_pred,
    'confusion_matrix' = nb_cf
  )

  nb_perf_metr_list[[i]] <- list(
    'Accuracy' = nb_acc,
    'Overall_misc_rates' = nb_mcr,
    "Sensitivity" = nb_sens,
    "Specificity" = nb_spec,
    'Pos_pred_v' = nb_ppv,
    'Neg_pred_v' = nb_npv,
    'AUC' = nb_auc,
    'SomersD' = nb_somers
  )

  message(glue('{i} total naive bayes prediction finished'))
}

# Performances results averaged over boostrap samples
nb_total_perf <- bind_rows(nb_perf_metr_list) %>%
  summarise_all(.tbl = ., ~ mean(.))

# Individual predicted probability averaged over boostrap samples
nb_total_ind_pred_prob <- map_df(
  .x = nb_test_list, ~ .x[["boot_pred_df"]]
) %>%
  group_by(id) %>%
  summarise(average_pred_prob = mean(pred_prob))
