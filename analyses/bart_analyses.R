# BART analyses --------------------------------------------------------
# If doParallel", "caret", "glue" and "bartMachine" are not
# installed, please run the following command before executing the
# script:
# install.packages(
#   c("bartMachine", "doParallel", "caret", "glue"),
#   dependencies = TRUE
# )
library(tidyverse)
library(caret)
library(glue)
library(doParallel)

# The aim of this script is to reproduce the analyses with the BART
# models on the simulated data.
seed <- 554477

# BART computational settings ------------------------------------------
options(java.parameters = "-Xmx8g") # must be set initially
library(bartMachine)
cl <- detectCores() - 1
set_bart_machine_num_cores(cl)

# Load the simulated dataset in the global environment
load(here::here("data", "simulated_db.rda"))

# 1) Model with only clinical covariates -------------------------------
db_clinical <- db %>%
  dplyr::select(contains("clinical"), y)

method <- 'bartMachine'
k_folds <- 10
repetitions <- 10

# Parameters grid
bart_grid <-  expand.grid(
  num_trees = c(100, 200, 500),
  alpha = 0.95,
  beta = 2,
  k = c(-3, -0.5, 0.5 ,3),
  nu = c(3, 10)
)

# Set the seeds for resampling
seq_fold <- k_folds * repetitions
seeds <- c(
  map(
    seq_len(seq_fold),
    ~ seq.int(
      from = 1 + (. - 1)*nrow(bart_grid),
      to = .*nrow(bart_grid),
      by = 1
    )
  )
)

seeds[[seq_fold + 1]] <- sample.int(10001)

bart_cv <- train(
  y ~ .,
  data = db_clinical,
  method = method,
  trControl = trainControl(
    method = 'repeatedcv',
    number = k_folds,
    repeats = repetitions,
    seeds = seeds,
    classProbs = TRUE,
    savePredictions = 'final',
    allowParallel = TRUE,
    search = 'grid',
    selectionFunction = "best"
  ),
  tuneGrid = bart_grid,
  serialize = T,      # serialize the model for future sessions
  num_burn_in = 500,  # MCMC samples as burn-in
  num_iterations_after_burn_in = 2000 # MCMC draws from posterior
)

bart_final <- bart_cv

# Performances training
bart_train_prob <- predict(
  object = bart_final,
  new_data = db_clinical,
  type = 'prob'
)

bart_train_pred <- predict(
  object = bart_final,
  data   = db_clinical
)

cf_df <- tibble(
  pred = bart_train_pred,
  obs = db_clinical[["y"]]
)

bart_pred_prob <- bart_train_prob %>%
  select(yes) %>%
  as_vector()

bart_cf <- confusionMatrix(table(cf_df))

bart_train_perf_no_genetic <- list(
  'Accuracy' = bart_cf$overall['Accuracy'],
  'Overall_misc_rates' = 1 - bart_cf$overall['Accuracy'],
  "Sensitivity" = bart_cf$byClass['Sensitivity'],
  "Specificity" = bart_cf$byClass['Specificity'],
  'Pos_pred_v' = bart_cf$byClass['Pos Pred Value'],
  'Neg_pred_v' = bart_cf$byClass['Neg Pred Value'],
  'AUC' = pROC::auc(
    response = db_clinical[["y"]],
    predictor = bart_pred_prob
  ),
  'SomersD' = 2*(
    pROC::auc(
      response = db_clinical[["y"]],
      predictor = bart_pred_prob
    ) - 0.5
  )
)

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
bart_test_list <- NULL

bart_perf_metr_list <- NULL

for (i in 1:length(boot_s)) {
  bart_test_prob <- predict(
    object = bart_final,
    new_data = boot_s[[i]] %>% dplyr::select(-id),
    type = 'prob'
  )

  bart_test_pred <- predict(
    object = bart_final,
    new_data = boot_s[[i]] %>% dplyr::select(-id)
  )

  cf_df <- data_frame(
    pred = bart_test_pred,
    obs = boot_s[[i]][["y"]]
  )

  # Get predicted probabilties
  bart_test_prob <- bart_test_prob %>%
    select(yes) %>%
    as_vector()

  # Performance metrics
  bart_cf <- confusionMatrix(table(cf_df))

  bart_acc <- bart_cf$overall['Accuracy']

  bart_mcr <- 1 - bart_acc

  bart_sens <- bart_cf$byClass['Sensitivity']

  bart_spec <- bart_cf$byClass['Specificity']

  bart_ppv <- bart_cf$byClass['Pos Pred Value']

  bart_npv <- bart_cf$byClass['Neg Pred Value']

  bart_auc <- pROC::auc(
    response = boot_s[[i]][["y"]],
    predictor = bart_test_prob
  )

  bart_somers <- 2*(bart_auc - 0.5)

  # Attach predicted probabilities and predicted class to bootstrap
  # sample
  boot_s[[i]] <- boot_s[[i]] %>%
    mutate(
      pred_prob = bart_test_prob,
      pred_class = bart_test_pred
    )

  bart_test_list[[i]] <- list(
    'boot_pred_df' = boot_s[[i]],
    'individual_prob' = bart_test_prob,
    'predicted_values' = bart_test_pred,
    'confusion_matrix' = bart_cf
  )

  bart_perf_metr_list[[i]] <- list(
    'Accuracy' = bart_acc,
    'Overall_misc_rates' = bart_mcr,
    'Pos_pred_v' = bart_ppv,
    'Neg_pred_v' = bart_npv,
    'AUC' = bart_auc,
    'SomersD' = bart_somers
  )

  message(glue::glue('{i} clinical bart prediction finished'))
}

# Performances results averaged over boostrap samples
bart_clinical_perf <- bind_rows(bart_perf_metr_list) %>%
  summarise_all(.tbl = ., ~ mean(.))

# Individual predicted probability averaged over boostrap samples
bart_clinical_ind_pred_prob <- map_df(
  .x = bart_test_list, ~ .x[["boot_pred_df"]]
) %>%
  group_by(id) %>%
  summarise(average_pred_prob = mean(pred_prob))


# 2) Model with clinical covariates and genes --------------------------
method <- 'bartMachine'
k_folds <- 10
repetitions <- 10

# Parameters grid
bart_grid <-  expand.grid(
  num_trees = c(100, 200, 500),
  alpha = 0.95,
  beta = 2,
  k = c(-3, -0.5, 0.5 ,3),
  nu = c(3, 10)
)

# Set the seeds for resampling
seq_fold <- k_folds * repetitions
seeds <- c(
  map(
    seq_len(seq_fold),
    ~ seq.int(
      from = 1 + (. - 1)*nrow(bart_grid),
      to = .*nrow(bart_grid),
      by = 1
    )
  )
)

seeds[[seq_fold + 1]] <- sample.int(10001)

bart_cv <- train(
  y ~ .,
  data = db,
  method = method,
  trControl = trainControl(
    method = 'repeatedcv',
    number = k_folds,
    repeats = repetitions,
    seeds = seeds,
    classProbs = TRUE,
    savePredictions = 'final',
    allowParallel = TRUE,
    search = 'grid',
    selectionFunction = "best"
  ),
  tuneGrid = bart_grid,
  serialize = T,      # serialize the model for future sessions
  num_burn_in = 500,  # MCMC samples as burn-in
  num_iterations_after_burn_in = 2000 # MCMC draws from posterior
)

bart_final <- bart_cv

# Performances training
bart_train_prob <- predict(
  object = bart_final,
  new_data = db,
  type = 'prob'
)

bart_train_pred <- predict(
  object = bart_final,
  data   = db
)

cf_df <- tibble(
  pred = bart_train_pred,
  obs = db[["y"]]
)

bart_pred_prob <- bart_train_prob %>%
  select(yes) %>%
  as_vector()

bart_cf <- confusionMatrix(table(cf_df))

bart_train_perf_no_genetic <- list(
  'Accuracy' = bart_cf$overall['Accuracy'],
  'Overall_misc_rates' = 1 - bart_cf$overall['Accuracy'],
  "Sensitivity" = bart_cf$byClass['Sensitivity'],
  "Specificity" = bart_cf$byClass['Specificity'],
  'Pos_pred_v' = bart_cf$byClass['Pos Pred Value'],
  'Neg_pred_v' = bart_cf$byClass['Neg Pred Value'],
  'AUC' = pROC::auc(
    response = db[["y"]],
    predictor = bart_pred_prob
  ),
  'SomersD' = 2*(
    pROC::auc(
      response = db[["y"]],
      predictor = bart_pred_prob
    ) - 0.5
  )
)

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
bart_test_list <- NULL

bart_perf_metr_list <- NULL

for (i in 1:length(boot_s)) {
  bart_test_prob <- predict(
    object = bart_final,
    new_data = boot_s[[i]] %>% dplyr::select(-id),
    type = 'prob'
  )

  bart_test_pred <- predict(
    object = bart_final,
    new_data = boot_s[[i]] %>% dplyr::select(-id)
  )

  cf_df <- data_frame(
    pred = bart_test_pred,
    obs = boot_s[[i]][["y"]]
  )

  # Get predicted probabilties
  bart_test_prob <- bart_test_prob %>%
    select(yes) %>%
    as_vector()

  # Performance metrics
  bart_cf <- confusionMatrix(table(cf_df))

  bart_acc <- bart_cf$overall['Accuracy']

  bart_mcr <- 1 - bart_acc

  bart_sens <- bart_cf$byClass['Sensitivity']

  bart_spec <- bart_cf$byClass['Specificity']

  bart_ppv <- bart_cf$byClass['Pos Pred Value']

  bart_npv <- bart_cf$byClass['Neg Pred Value']

  bart_auc <- pROC::auc(
    response = boot_s[[i]][["y"]],
    predictor = bart_test_prob
  )

  bart_somers <- 2*(bart_auc - 0.5)

  # Attach predicted probabilities and predicted class to bootstrap
  # sample
  boot_s[[i]] <- boot_s[[i]] %>%
    mutate(
      pred_prob = bart_test_prob,
      pred_class = bart_test_pred
    )

  bart_test_list[[i]] <- list(
    'boot_pred_df' = boot_s[[i]],
    'individual_prob' = bart_test_prob,
    'predicted_values' = bart_test_pred,
    'confusion_matrix' = bart_cf
  )

  bart_perf_metr_list[[i]] <- list(
    'Accuracy' = bart_acc,
    'Overall_misc_rates' = bart_mcr,
    'Pos_pred_v' = bart_ppv,
    'Neg_pred_v' = bart_npv,
    'AUC' = bart_auc,
    'SomersD' = bart_somers
  )

  message(glue::glue('{i} total bart prediction finished'))
}

# Performances results averaged over boostrap samples
bart_total_perf <- bind_rows(bart_perf_metr_list) %>%
  summarise_all(.tbl = ., ~ mean(.))

# Individual predicted probability averaged over boostrap samples
bart_total_ind_pred_prob <- map_df(
  .x = bart_test_list, ~ .x[["boot_pred_df"]]
) %>%
  group_by(id) %>%
  summarise(average_pred_prob = mean(pred_prob))
