# Simulated data -------------------------------------------------------
# If "tidyverse" and "here" packages are not installed, please run the
# following command before executing the script:
# install.packages(c("tidyverse", "here"), dependencies = TRUE)
library(tidyverse)

# The aim of this script is to simulate a dataset similar to those
# used in the analysis. The dataset will be used to reproduce the
# analysis of the paper.

seed <- 554477

# Set a seed for reproducibility of the data
set.seed(seed)

# All the variables are factors. Simulate 7 clinical variables and 8
# genes. All the variables are simulated from a binomial distribution.
n <- 200L

set.seed(seed)
df_covs <- tibble(
  clinical_1 = factor(rbinom(n = n, size = 1, prob = 0.3)),
  clinical_2 = factor(rbinom(n = n, size = 2, prob = 0.3)),
  clinical_3 = factor(rbinom(n = n, size = 3, prob = 0.4)),
  clinical_4 = factor(rbinom(n = n, size = 1, prob = 0.2)),
  clinical_5 = factor(rbinom(n = n, size = 1, prob = 0.3)),
  clinical_6 = factor(rbinom(n = n, size = 2, prob = 0.4)),
  clinical_7 = factor(rbinom(n = n, size = 1, prob = 0.3)),
  gene_1 = factor(rbinom(n = n, size = 2, prob = 0.6)),
  gene_2 = factor(rbinom(n = n, size = 2, prob = 0.3)),
  gene_3 = factor(rbinom(n = n, size = 1, prob = 0.3)),
  gene_4 = factor(rbinom(n = n, size = 2, prob = 0.3)),
  gene_5 = factor(rbinom(n = n, size = 2, prob = 0.5)),
  gene_6 = factor(rbinom(n = n, size = 2, prob = 0.2)),
  gene_7 = factor(rbinom(n = n, size = 1, prob = 0.3)),
  gene_8 = factor(rbinom(n = n, size = 4, prob = 0.4))
)

# Simulate the outcome from a binomial distribution. The probability
# will be computed following a logit model with main effects and two-way
# interactions
design_matrix <- model.matrix(~ . * ., data = df_covs)

# Simulate coefficients from a uniform distribution bounded between
# -0.8 and 0.8
set.seed(seed)
coefs <- runif(n = ncol(design_matrix) - 1, min = -0.8, max = 0.8)

# Intercept such that the proportion of the outcome is approximately 0.5
int <- 0.8

# Simulate the linear predictor, transform it and simulate the outcome
inv_logit <- function(x) {exp(x)/(1 + exp(x))}

eta <- design_matrix %*% c(int, coefs)
mu <- inv_logit(eta)
set.seed(seed)
y <- rbinom(n = n, size = 1, prob = mu)

# Add the outcome to the tibble as factor (absence and presence of
# an hypothetical disease) and save the data
db <- df_covs %>%
  mutate(
    y = factor(if_else(y == 0, "no", "yes"), levels = c("no", "yes"))
  )

save(db, file = here::here("data", "simulated_db.rda"))
