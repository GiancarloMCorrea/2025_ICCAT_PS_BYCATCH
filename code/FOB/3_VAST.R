rm(list = ls())
# -------------------------------------------------------------------------
library(VAST)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(boot)
library(ggeffects)
library(pdp)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'VAST')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'VAST')
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
user_region = readRDS(file = file.path('data', this_type, 'user_region.rds'))
MyGridSets = readRDS(file = file.path('data', this_type, 'MyGridSets.rds'))
load(file.path('data', this_type, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
cumsp_data = mod_data %>% group_by(sp_code) %>% summarise(value = sum(value))
cumsp_data = arrange(cumsp_data, desc(value))
this_sp = cumsp_data$sp_code[10] # change this index if want to make a loop
sp_data = mod_data %>% dplyr::filter(sp_code %in% this_sp) 

# Remove NA's:
sp_data = sp_data[complete.cases(sp_data), ]
# Rename using VAST format:
sp_data = sp_data %>% dplyr::rename(Year = year, Lat = latitude, Lon = longitude,
                                      Catch = value, tunaCatch = tons_target_tuna)
# Define variable type:
sp_data = sp_data %>% mutate(Year = as.integer(Year), AreaSwept_km2 = 1)
glimpse(sp_data)

# Make settings:
settings <- make_settings(n_x = 100, Region='User',
                          purpose = "index2", bias.correct = TRUE,
                          use_anisotropy = FALSE, 
                          fine_scale = TRUE,
                          knot_method = 'grid'
                          # Options = c( "Calculate_Range"=TRUE, 
                          #              "Calculate_effective_area"=TRUE, 
                          #              "treat_nonencounter_as_zero"=TRUE,
                          #              "SD_site_density" = TRUE,
                          #              "SD_site_logdensity" = TRUE)
                          )
settings$ObsModel = c(4,4) # lognormal = 4 and logit-link = 4 (fix 1 and 0)
settings$FieldConfig[1:2,1] = c(0, 0) # no spatial or spatiotemporal effect comp 1
# settings$FieldConfig[3,] = 0 # fixed effects betas
settings$FieldConfig
settings$RhoConfig
fit <- fit_model(settings=settings,
                 Lat_i=sp_data$Lat, Lon_i=sp_data$Lon,
                 t_i=sp_data$Year, b_i=sp_data$Catch,
                 a_i=sp_data$AreaSwept_km2,
                 input_grid=user_region)
dir.create(file.path(model_folder, this_sp), showWarnings = FALSE)
save(fit, file = file.path(model_folder, this_sp, 'fit.RData'))

# Plot results:
dir.create(file.path(plot_folder, this_sp), showWarnings = FALSE)
plot_results(fit, plot_set=c(3), working_dir = file.path(getwd(), plot_folder, this_sp))

# Calculate annual bycatch estimates (by hand, TODO: calculate SE):
nset_matrix = MyGridSets %>% select(ID, year, n_sets) %>% 
  tidyr::pivot_wider(names_from = year, values_from = n_sets)
nset_matrix = as.matrix(nset_matrix)
rownames(nset_matrix) = nset_matrix[,1]
nset_matrix = nset_matrix[,-1]
# Order rows:
nset_matrix = nset_matrix[match(user_region$ID, rownames(nset_matrix)), ]
# Order columns:
nset_matrix = nset_matrix[ , match(fit$year_labels, colnames(nset_matrix))]
# Multiply by D_gct:
bycatch_est = colSums(nset_matrix * fit$Report$D_gct[,1,])
bycatch_est = data.frame(year = names(bycatch_est), est = as.vector(bycatch_est))
plot(bycatch_est$year, bycatch_est$est, type = 'b')

# Ratio estimator (confirm?)
sp_data %>% group_by(Year, ID) %>% summarise(Catch = mean(Catch))

# -------------------------------------------------------------------------
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

# -------------------------------------------------------------------------
# Plot covariates effects (ggeffects):

# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = NULL

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( mod = fit,
                         focal.predictors = c("sst"),
                         which_formula = "X2",
                         transformation = list(link=identity, inverse=identity) )
plot(pred)

# -------------------------------------------------------------------------
# Plot covariates effects (pdp):

# Make function to interface with pdp
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = newdata,
           do_checks = FALSE )
}

# Run partial
Partial = partial( object = fit,
                   pred.var = "sst",
                   pred.fun = pred.fun,
                   train = fit$covariate_data )

# Make plot using ggplot2
autoplot(Partial)
