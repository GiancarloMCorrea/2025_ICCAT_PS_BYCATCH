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
sel_sp = 'Carcharhinus falciformis'
this_type = 'FOB'
min_year = 2015

# Define folder to save results
plot_folder = file.path('figures', this_type, 'VAST', sel_sp)
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'VAST', sel_sp)
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
user_region = readRDS(file = file.path('data', this_type, 'user_region.rds'))
load(file.path('data', this_type, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Run model:
mod_data = weight_data %>% dplyr::filter(year >= min_year, sp_code %in% sel_sp) 
mod_data = mod_data %>% mutate(year = as.numeric(year), AreaSwept_km2 = 1)

settings <- make_settings(n_x = 689, Region='User',
                          purpose = "index2", bias.correct = FALSE,
                          use_anisotropy = FALSE, 
                          fine_scale = TRUE,
                          knot_method = 'grid')
settings$ObsModel = c(4,0) # lognormal = 4 and logit-link = 0
settings$FieldConfig[1:2,1] = c(0, 0) # no spatial or spatiotemporal effect comp 1
settings$FieldConfig
settings$RhoConfig
fit <- fit_model(settings=settings,
                 Lat_i=mod_data$latitude, Lon_i=mod_data$longitude,
                 t_i=mod_data$year, b_i=mod_data$value,
                 a_i=mod_data$AreaSwept_km2,
                 input_grid=user_region)
save(fit, file = file.path(model_folder, 'fit.RData'))

# Plot results:
plot_results(fit, plot_set=3, working_dir = file.path(getwd(), plot_folder))

# Plot Omega 2nd component:
plot_dat_omega = data.frame(lon = fit$extrapolation_list$Data_Extrap$Lon, 
                            lat = fit$extrapolation_list$Data_Extrap$Lat, 
                            omega = fit$Report$Omega2_gc[,1])
p1 = ggplot(data = plot_dat_omega, aes(x = lon, y = lat)) +
  geom_point(aes(color = omega), size = 1) +
  scale_colour_gradient2() + 
  coord_fixed() + theme(legend.position = 'bottom', legend.direction="horizontal")
ggsave(filename = paste0('spatial_effect', img_type), path = plot_folder, plot = p1, 
       width = img_width*0.5, height = 75, units = 'mm', dpi = img_res)

# Plot Epsilon 2nd component:
epsilon_df = list()
for(i in seq_along(fit$year_labels)) {
  epsilon_df[[i]] = data.frame(year = fit$year_labels[i],
                               lon = fit$extrapolation_list$Data_Extrap$Lon, 
                               lat = fit$extrapolation_list$Data_Extrap$Lat, 
                               epsilon = fit$Report$Epsilon2_gct[,1,i])
}
plot_dat_epsilon = bind_rows(epsilon_df)
p1 = ggplot(data = plot_dat_epsilon, aes(x = lon, y = lat)) +
  geom_point(aes(color = epsilon), size = 1) +
  scale_colour_gradient2() +
  coord_fixed() + theme(legend.position = 'bottom', legend.direction="horizontal") +
  facet_wrap(~year)
ggsave(filename = paste0('sptemp_effect', img_type), path = plot_folder, plot = p1, 
       width = img_width*0.75, height = 140, units = 'mm', dpi = img_res)

# Plot predicted values: invlogit(comp1)*exp(comp2):
tmp_df = list()
for(i in seq_along(fit$year_labels)) {
  tmp_df[[i]] = data.frame(ID = fit$extrapolation_list$Data_Extrap$ID,
                           year = fit$year_labels[i],
                           lon = fit$extrapolation_list$Data_Extrap$Lon, 
                           lat = fit$extrapolation_list$Data_Extrap$Lat, 
                           cpue_pred = as.numeric(fit$Report$D_gct[,1,i]))
}
pred_df = bind_rows(tmp_df)
PredGrid = left_join(MyGrid, pred_df, by = c('ID')) %>% na.omit
p1 = plot_predictions(PredGrid, legend_position = 'bottom', nCol = 3)
ggsave(filename = paste0('map_predictions', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 140, units = 'mm', dpi = img_res)
