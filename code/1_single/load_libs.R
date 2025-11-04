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
library(ggh4x)
library(foreach)
library(snowfall)
library(scales)
library(here)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# SELECT YOUR SET TYPE HERE!:
this_type = 'FOB' # FOB or FSC

# First year to make analyses:
str_yr = 2015

# Define plot and data folder:
# Data folder:
data_folder = here(file.path("data", "1_single", this_type))
dir.create(data_folder, showWarnings = FALSE, recursive = TRUE)
# Plot folder:
plot_folder = here(file.path("figures", "1_single", this_type))
dir.create(plot_folder, showWarnings = FALSE, recursive = TRUE)
# Model folder:
model_folder = here(file.path("model", "1_single", this_type))
dir.create(model_folder, showWarnings = FALSE, recursive = TRUE)
