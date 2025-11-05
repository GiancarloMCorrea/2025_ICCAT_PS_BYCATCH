# Read libraries:
library(ggplot2)
library(dplyr)  
library(sf)
library(tibble)
library(readr)
library(purrr)
library(gridExtra)
library(viridis)
library(sdmTMB)
library(DescTools)
library(boot)
library(stars)
library(lubridate)
library(blockCV)
library(fmesher)
library(corrplot)
library(car)
library(future)
library(tidyr)
library(reshape2)
require(here)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# SELECT YOUR SET TYPE HERE!:
this_type = 'FOB' # FOB or FSC

# Define plot and data folder:
# Data folder:
data_folder = here(file.path("data", "3_sampsize", this_type))
dir.create(data_folder, showWarnings = FALSE, recursive = TRUE)
# Plot folder:
plot_folder = here(file.path("figures", "3_sampsize", this_type))
dir.create(plot_folder, showWarnings = FALSE, recursive = TRUE)
# Model folder:
model_folder = here(file.path("model", "3_sampsize", this_type))
dir.create(model_folder, showWarnings = FALSE, recursive = TRUE)
