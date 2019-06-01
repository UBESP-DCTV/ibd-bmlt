# ibd-bmlt
This project contains the R scripts and the data to reproduce the
statistical analysis of the study "__Role of Genetic Factors in
Characterizing Extra-Intestinal Manifestations in Crohn’s Disease 
Patients: Are Bayesian Machine-Learning Methods Improving Outcome
Predictions?__".

Before running the analysis, please install the following R packages:

``` r
install.packages(
  c("tidyverse", "doParallel", "caret", "glue", "bnlearn",
    "bartMachine", "graph"), dependencies = TRUE
)
```

The project is organized as follows:

- The file _ibd_bmlt.Rproj_ is the file of the Rstudio project.

- The folder _data_ contains the following files:

    * _sim_data.R_, which contains the R script used to simulate the
      data
      
    * _simulated_db.rda_, which is the file that contains the simulated
      dataset used to reproduce the analyses
      
- The folder _analyses_ contains the following files:

    * _nb_analyses.R_, which contains the R script used to run the
      analyses with Naive-Bayes models
      
    * _bn_analyses.R_, which contains the R script used to run the
      analyses with Bayesian Network models
      
    * _bart_analyses.R_, which contains the R script used to run the
      analyses with BART models
      
    * _misc_functions.R_, which contains the R script with 
      various functions used during the analyses
  
