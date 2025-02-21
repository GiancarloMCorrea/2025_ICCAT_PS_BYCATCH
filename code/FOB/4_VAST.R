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
extraRegion_VAST = readRDS(file = file.path('data', this_type, 'extraRegion_VAST.rds'))
extraRegion_tinyVAST = readRDS(file = file.path('data', this_type, 'extraRegion_tinyVAST.rds'))

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
cumsp_data = mod_data %>% group_by(sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = arrange(cumsp_data, desc(Catch))

n_mod_sp = 5 # model the first N species
for(i in 1:n_mod_sp) {
  
  # Create dir to save plots and model outputs:
  this_sp = cumsp_data$sp_name[i]
  dir.create(file.path(model_folder, this_sp), showWarnings = FALSE)
  dir.create(file.path(plot_folder, this_sp), showWarnings = FALSE)
  
  # Filter data:
  sp_data = mod_data %>% dplyr::filter(sp_name %in% this_sp) 
  # Prepare data for model:
  sp_data = as.data.frame(sp_data)
  
  # Make settings:
  settings <- make_settings(n_x = 100, Region='User',
                            purpose = "index2", bias.correct = TRUE,
                            use_anisotropy = FALSE, 
                            fine_scale = TRUE,
                            knot_method = 'grid'
                            )
  settings$ObsModel = c(4,1) # lognormal = 4 and poisson-link delta = 1 
  settings$FieldConfig[1:2,1] = c(0, 0) # no spatial or spatiotemporal effect comp 1
  # settings$FieldConfig[3,] = 0 # fixed effects betas
  # settings$RhoConfig["Epsilon2"] = 4 # AR1 on epsilon2
  VModel <- fit_model(settings=settings,
                   Lat_i=sp_data$Lat, Lon_i=sp_data$Lon,
                   t_i=sp_data$Year, b_i=sp_data$Catch,
                   a_i=sp_data$AreaSwept_km2,
                   input_grid=extraRegion_VAST)
  save(VModel, file = file.path(model_folder, this_sp, 'fit.RData'))
  
  # Plot results:
  plot_results(VModel, plot_set=c(3), working_dir = file.path(getwd(), plot_folder, this_sp))
  
  # -------------------------------------------------------------------------
  # Calculate annual bycatch estimates by hand 
  # VAST cannot do it automatically since index is weighted by area (time-invariant)
  # SE cannot be estimated. Do it using tinyVAST
  nset_matrix = extraRegion_tinyVAST %>% select(ID, Year, n_sets) %>% 
    tidyr::pivot_wider(names_from = Year, values_from = n_sets)
  nset_matrix = as.matrix(nset_matrix)
  rownames(nset_matrix) = nset_matrix[,1]
  nset_matrix = nset_matrix[,-1]
  # Order rows and columns:
  nset_matrix = nset_matrix[match(extraRegion_VAST$ID, rownames(nset_matrix)), ]
  nset_matrix = nset_matrix[ , match(VModel$year_labels, colnames(nset_matrix))]
  # Multiply by D_gct:
  bycatch_est = colSums(nset_matrix * VModel$Report$D_gct[,1,])
  Index = data.frame(Year = as.integer(names(bycatch_est)), est = as.vector(bycatch_est))
  Index$model = 'VAST'
  Index$type = 'single'
  Index$species = this_sp
  write.csv(Index, file = file.path(plot_folder, this_sp, 'Bycatch_est_ts.csv'), row.names = FALSE)
  
  # Plot annual estimates with observed estimate:
  obs_df = read.csv(file = file.path('data', this_type, 'bycatch_est_obs', paste0(this_sp, '.csv')))
  p1 = plot_time_predictions(Index, obs_df, Year, est, yLab = 'Estimated bycatch (tons)')
  ggsave(filename = paste0('Bycatch_est_ts', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width*0.5, height = 70, units = 'mm', dpi = img_res)
  
  # -------------------------------------------------------------------------
  # Plot Omega 2nd component:
  plot_dat = data.frame(Lon = VModel$spatial_list$latlon_s[,2], 
                        Lat = VModel$spatial_list$latlon_s[,1], 
                        omega2 = VModel$Report$Omega2_sc[,1])
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = omega2), size = 2) + scale_colour_gradient2() + labs(color = 'Omega2')
  p1 = add_sf_map(p1)  
  ggsave(filename = paste0('omega2', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)


  # -------------------------------------------------------------------------
  # Plot Epsilon 2nd component:
  tmp_df = list()
  for(i in seq_along(VModel$year_labels)) {
    tmp_df[[i]] = data.frame(year = VModel$year_labels[i],
                             Lon = VModel$spatial_list$latlon_s[,2], 
                             Lat = VModel$spatial_list$latlon_s[,1], 
                             epsilon2 = VModel$Report$Epsilon2_sct[,1,i])
  }
  plot_dat = bind_rows(tmp_df)
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = epsilon2), size = 1) + scale_colour_gradient2() + labs(color = 'Epsilon2')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ year)
  ggsave(filename = paste0('epsilon2', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)

  # # -------------------------------------------------------------------------
  # # Plot covariates effects (ggeffects):
  # 
  # # Must add data-frames to global environment (hope to fix in future)
  # covariate_data_full = fit$effects$covariate_data_full
  # catchability_data_full = NULL
  # 
  # # Plot 1st linear predictor, but could use `transformation` to apply link function
  # pred = Effect.fit_model( mod = fit,
  #                          focal.predictors = c("sst"),
  #                          which_formula = "X2",
  #                          transformation = list(link=identity, inverse=identity) )
  # plot(pred)
  # 
  # # -------------------------------------------------------------------------
  # # Plot covariates effects (pdp):
  # # Make function to interface with pdp
  # pred.fun = function( object, newdata ){
  #   predict( x=object,
  #            Lat_i = object$data_frame$Lat_i,
  #            Lon_i = object$data_frame$Lon_i,
  #            t_i = object$data_frame$t_i,
  #            a_i = object$data_frame$a_i,
  #            what = "P1_iz",
  #            new_covariate_data = newdata,
  #            do_checks = FALSE )
  # }
  # # Run partial
  # Partial = partial( object = fit,
  #                    pred.var = "sst",
  #                    pred.fun = pred.fun,
  #                    train = fit$covariate_data )
  # 
  # # Make plot using ggplot2
  # autoplot(Partial)

}
