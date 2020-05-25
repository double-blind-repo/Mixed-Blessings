# Loading Packages ----
list.of.packages <- c("broom","curl","dotwhisker","dplyr","essurvey","fastDummies","GPArotation",
                     "gridExtra","lavaan","mice","parallel","psych","stringr","xtable") 

lapply(list.of.packages, library, character.only = TRUE)

# Preparing the data.
source("setup.R")

# Main Analysis using Imputation
source("imputation_mice.R")

# List-Wise Version of Analysis
source("listwise.R")

# Supplementary Analysis using Imputation
source("supplementary_mice.R")