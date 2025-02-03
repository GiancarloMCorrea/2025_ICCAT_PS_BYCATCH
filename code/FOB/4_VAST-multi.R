rm(list = ls())
# -------------------------------------------------------------------------
library(VAST)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(boot)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
n_sp = 15 # first N species (only applies to weight or numbers datasets)
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'VAST-multi')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'VAST-multi')
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
user_region = readRDS(file = file.path('data', this_type, 'user_region.rds'))
load(file.path('data', this_type, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Filter first N species:
mod_data = weight_data
cumsp_data = mod_data %>% group_by(sp_code) %>% summarise(value = sum(value))
cumsp_data = arrange(cumsp_data, desc(value))
cumsp_data = cumsp_data[1:n_sp, ]
mod_data = mod_data %>% dplyr::filter(sp_code %in% cumsp_data$sp_code) %>%
              mutate(sp_number = factor(sp_code, levels = cumsp_data$sp_code, labels = 1:n_sp))

# Select variables to be used in the model:
mod_data = mod_data %>% select(ID, id_set, year, month, vessel_code, flag_country,
                               latitude, longitude, sst, sunrise_diference, tuna_catch,
                               yft_catch, bet_catch, skj_catch, sp_code, sp_number, value)
# Remove NA's:
mod_data = mod_data[complete.cases(mod_data), ]
# Rename using VAST format:
mod_data = mod_data %>% dplyr::rename(Year = year, Lat = latitude, Lon = longitude,
                                      Catch = value, tSunrise = sunrise_diference)
# Define variable type:
mod_data = mod_data %>% mutate(sp_number = as.numeric(as.character(sp_number)),
                               tSunrise = as.numeric(tSunrise),
                               Year = as.integer(Year), AreaSwept_km2 = 1)
glimpse(mod_data)

# Make settings:
settings <- make_settings(n_x = 100, Region='User',
                          purpose = "ordination", bias.correct = FALSE,
                          use_anisotropy = FALSE, 
                          fine_scale = TRUE,
                          knot_method = 'grid',
                          n_categories = 2)

settings$ObsModel = c(10,2) # Tweedie
settings$FieldConfig[c('Omega','Epsilon','Beta'),'Component_1'] = c(0,0,'IID')
settings$RhoConfig[c('Beta1','Epsilon1')] = c(3,0)
# Specify full rank spatial term and eliminate spatio-temporal term
# this makes the demo run faster
settings$FieldConfig[c('Omega','Epsilon','Beta'),'Component_2'] = c(2,0,'IID')
settings$RhoConfig[c('Beta2','Epsilon2')] = c(0,0)

fit <- fit_model(settings=settings,
                 Lat_i=mod_data$Lat, 
                 Lon_i=mod_data$Lon,
                 t_i=mod_data$Year, 
                 b_i=mod_data$Catch,
                 a_i=mod_data$AreaSwept_km2,
                 c_i=mod_data$sp_number-1,
                 covariate_data = as.data.frame(mod_data[,c('Lat', 'Lon', 'Year', 'tuna_catch', 'sst', 'tSunrise')]),
                 X2_formula = ~ tuna_catch + sst + tSunrise,
                 newtonsteps = 0,
                 getsd = FALSE,
                 input_grid=user_region)
save(fit, file = file.path(model_folder, 'fit.RData'))

# Plot results
results = plot( fit,
                plot_set = c(3,17),
                category_names = cumsp_data$sp_code,
                working_dir = file.path(getwd(), plot_folder))
