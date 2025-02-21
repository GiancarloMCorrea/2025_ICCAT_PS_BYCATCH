rm(list = ls())
# -------------------------------------------------------------------------
library(sdmTMB)
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
plot_folder = file.path('figures', this_type, 'sdmTMB', sel_sp)
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'sdmTMB', sel_sp)
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type,'weight_data.rds'))
load(file.path('data', this_type, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Filter:
mod_data = weight_data %>% dplyr::filter(year >= min_year, sp_code %in% sel_sp) 
hist(mod_data$value[mod_data$value > 0])

# -------------------------------------------------------------------------
# Make mesh:
# Make some changes before running model:
# mod_data = sdmTMB::add_utm_columns(dat = mod_data, ll_names = c('longitude', 'latitude'), units = 'km')
# FishStatsUtils::project_coordinates(X = mod_data$longitude, Y = mod_data$latitude)
# oce::lonlat2utm(longitude = mod_data$longitude, latitude = mod_data$latitude, km = TRUE)

# Make mesh:
# mesh = make_mesh(mod_data, c("X", "Y"), cutoff = 75) # larger cutoff = less complexity
mesh = make_mesh(mod_data, c("longitude", "latitude"), cutoff = 1) # larger cutoff = less complexity
png(filename = file.path(plot_folder, 'mesh.png'), width = 95, height = 100, units = 'mm', res = 500)
par(mar = c(1, 1, 1, 1))
plot(mesh)
dev.off()


# -------------------------------------------------------------------------
# Spatiotemporal model:
mod_data = mod_data %>% mutate(time_factor = factor(year))
my_st_model <- sdmTMB(
  data = mod_data,
  formula = list(
    value ~ 0 + time_factor, # presence/absence
    value ~ 0 + time_factor # positive 
  ),
  mesh = mesh, 
  time = 'year',
  family = delta_lognormal(),
  spatial = list("off", "on"),
  spatiotemporal = list("off", 'iid')
)
save(my_st_model, file = file.path(model_folder, 'my_model.RData'))

# Get residuals:
mod_data$resids1 = residuals(my_st_model, type = "mle-mvn", model = 1) # randomized quantile residuals
mod_data$resids2 = residuals(my_st_model, type = "mle-mvn", model = 2) # randomized quantile residuals

png(filename = file.path(plot_folder, 'rand_residuals.png'), width = 170, height = 170, units = 'mm', res = 500)
par(mfrow = c(2,2))
hist(mod_data$resids1, main = 'Component 1')
qqnorm(mod_data$resids1, main = 'Component 1')
abline(a = 0, b = 1)
hist(mod_data$resids2, main = 'Component 2')
qqnorm(mod_data$resids2, main = 'Component 2')
abline(a = 0, b = 1)
dev.off()

# Check by time:
p1 = ggplot(mod_data, aes(longitude, latitude, col = resids1)) + 
  scale_colour_gradient2(name = 'Resid (comp 1)') +
  geom_point(size = 0.5) + facet_wrap(~year, ncol = 3) + 
  coord_fixed() + theme(legend.position = 'bottom', legend.direction="horizontal") 
ggsave(filename = paste0('residuals_spacetime_1', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 140, units = 'mm', dpi = img_res)
p2 = ggplot(mod_data %>% dplyr::filter(value > 0), aes(longitude, latitude, col = resids2)) + 
  scale_colour_gradient2(name = 'Resid (comp 2)') +
  geom_point(size = 0.5) + facet_wrap(~year, ncol = 3) + 
  coord_fixed() + theme(legend.position = 'bottom', legend.direction="horizontal") 
ggsave(filename = paste0('residuals_spacetime_2', img_type), path = plot_folder, plot = p2, 
       width = 170, height = 140, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Predictions:
varData = expand.grid(time_factor = as.factor(levels(unique(my_st_model$data$time_factor))),
                      ID = unique(my_st_model$data$ID))
varData = varData %>% mutate(year = as.character(time_factor))
MyGrid2 = st_centroid(MyGrid) %>% dplyr::mutate(longitude = sf::st_coordinates(.)[,1], latitude = sf::st_coordinates(.)[,2])
st_geometry(MyGrid2) = NULL
predData = left_join(varData, MyGrid2[,c('ID', 'longitude', 'latitude')], by = 'ID')
predictions = predict(my_st_model, newdata = predData, return_tmb_object = TRUE)
PredGrid = left_join(MyGrid, predictions$data, by = c('ID')) %>% na.omit
PredGrid$cpue_pred = boot::inv.logit(PredGrid$est1) * exp(PredGrid$est2)
saveRDS(PredGrid, file = file.path(model_folder, 'predictions.rds'))

# Plot spatial random effects:
p1 = ggplot(predictions$data, aes(longitude, latitude, col = omega_s2)) + 
  scale_colour_gradient2() +
  geom_point(size = 0.5) + facet_wrap(~year, ncol = 3) + 
  coord_fixed() + theme(legend.position = 'bottom', legend.direction="horizontal")
ggsave(filename = paste0('spatial_effect', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 140, units = 'mm', dpi = img_res)

# Plot spatiotemporal random effects:
p1 = ggplot(predictions$data, aes(longitude, latitude, col = epsilon_st2)) + 
  scale_colour_gradient2() +
  geom_point(size = 0.5) + facet_wrap(~year, ncol = 3) + 
  coord_fixed() + theme(legend.position = 'bottom', legend.direction="horizontal")
ggsave(filename = paste0('sptemp_effect', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 140, units = 'mm', dpi = img_res)

# Plot CPUE by year:
plot_data = PredGrid %>% group_by(year, ID) %>% summarise(cpue_pred = mean(cpue_pred)) # aggregate by year and grid, just to display in a map
p1 = plot_predictions(plot_data, legend_position = 'bottom', nCol = 3)
ggsave(filename = paste0('map_predictions_hurdle', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 140, units = 'mm', dpi = img_res)
