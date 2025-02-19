rm(list = ls())
require(ggplot2)
require(sf)
require(dplyr)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Define type of school to be analyzed:
this_type = 'FOB'

# -------------------------------------------------------------------------
# Data and plot folders:
save_data_folder = file.path('data', this_type)
dir.create(save_data_folder, recursive = TRUE, showWarnings = FALSE)
plot_dir = file.path('figures', this_type)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load data
load(file.path(save_data_folder, 'obsPoints.RData'))
st_geometry(obsPoints) = NULL
my_data = obsPoints

# -------------------------------------------------------------------------
# Filter data based on some criteria:
summary(my_data)
my_data = my_data %>% dplyr::filter(year >= 2013) # same as effort data

# -------------------------------------------------------------------------
# # define sp to be deleted from dataset (usually target species)
# del_sp = c('Katsuwonus pelamis', 'Thunnus albacares', 'Thunnus obesus',
#            'Auxis rochei', 'Auxis thazard', 'Auxis thazard, A. rochei',
#            'Euthynnus affinis', 'Euthynnus alletteratus',
#            'Thunnus alalunga')
# # define sp to be used as covariates in the model (BET, YFT, SKJ)
# cov_sp = c('Thunnus albacares', 'Thunnus obesus', 'Katsuwonus pelamis')

# -------------------------------------------------------------------------
# 'number' data frame:
del_cols = c('weight_')
sel_col = 'number_'
tmp_data = my_data
# tmp_data = tmp_data %>% mutate(tuna_catch = rowSums(select(tmp_data, paste0('weight_', cov_sp))),
#                                yft_catch = rowSums(select(tmp_data, 'weight_Thunnus albacares')),
#                                bet_catch = rowSums(select(tmp_data, 'weight_Thunnus obesus')),
#                                skj_catch = rowSums(select(tmp_data, 'weight_Katsuwonus pelamis')))
tmp_data = tmp_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
numbers_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                   names_to = 'sp_name', values_to = 'value') %>% 
                  mutate(sp_name = gsub(pattern = sel_col, replacement = '', x = sp_name)) 
                  # dplyr::filter(!(sp_name %in% del_sp))
# Final edits:
numbers_data = numbers_data %>% dplyr::rename(Year = year, Lat = latitude, Lon = longitude,
                                            Catch = value, tunaCatch = tons_target_tuna)
numbers_data = numbers_data %>% mutate(Year = as.integer(Year), AreaSwept_km2 = 1)

# 'weight' data frame:
del_cols = c('number_')
sel_col = 'weight_'
tmp_data = my_data
tmp_data = tmp_data %>% select(-grep(pattern = paste(del_cols, collapse = "|"), x = colnames(my_data)))
weight_data = tidyr::pivot_longer(data = tmp_data, cols = grep(pattern = sel_col, x = colnames(tmp_data)),
                                  names_to = 'sp_name', values_to = 'value') %>% 
                  mutate(sp_name = gsub(pattern = sel_col, replacement = '', x = sp_name)) 
# Final edits:
weight_data = weight_data %>% dplyr::rename(Year = year, Lat = latitude, Lon = longitude,
                                    Catch = value, tunaCatch = tons_target_tuna)
weight_data = weight_data %>% mutate(Year = as.integer(Year), AreaSwept_km2 = 1)

# Save created data:
saveRDS(object = numbers_data, file = file.path(save_data_folder, 'numbers_data.rds'))
saveRDS(object = weight_data, file = file.path(save_data_folder, 'weight_data.rds'))
