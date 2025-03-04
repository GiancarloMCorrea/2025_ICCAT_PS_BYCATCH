rm(list = ls())
# -------------------------------------------------------------------------
library(tinyVAST)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(boot)
library(ggeffects)
library(pdp)
library(fmesher)
require(tidyr)
library(sdmTMB)
require(sem)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
n_sp = 20 # first N species (only applies to weight or numbers datasets)
n_fac_vec = 2:6 # Number of factors to test (vector)
this_type = 'FOB'

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
extraRegion_tinyVAST = readRDS(file = file.path('data', this_type, 'extraRegion_tinyVAST.rds'))


# -------------------------------------------------------------------------
# Main folders:
# Define folder to save results
main_plot_folder = file.path('figures', this_type, 'tinyVAST-multi')
dir.create(main_plot_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
cumsp_data = mod_data %>% group_by(sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = arrange(cumsp_data, desc(Catch))
cumsp_data = cumsp_data[1:n_sp, ]
sp_data = mod_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>%
  mutate(sp_number = factor(sp_name, levels = cumsp_data$sp_name, labels = 1:n_sp))
# Define sp ID:
sp_data = sp_data %>% mutate(sp_number = as.integer(as.character(sp_number)))
sp_data = as.data.frame(sp_data)

# Make mesh
my_mesh = fm_mesh_2d( sp_data[,c('Lon', 'Lat')], cutoff = 3 )
p1 = ggplot() + geom_fm(data = my_mesh)
ggsave(paste0('map_mesh', img_type), path = main_plot_folder, 
       plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)

# Plot mesh with observations:
mesh_df = data.frame(Lon = my_mesh$loc[,1], Lat = my_mesh$loc[,2]) %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
plot_dat = sp_data %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
p1 = ggplot() + 
  geom_sf(data = plot_dat, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
  geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
p1 = add_sf_map(p1)
ggsave(paste0('map_mesh_obs', img_type), path = main_plot_folder, 
       plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)

# Loop over number of factors
for(m in seq_along(n_fac)) {
  
  n_fac = n_fac_vec[m]
  # Define folder to save results
  plot_folder = file.path('figures', this_type, 'tinyVAST-multi', n_fac)
  dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Define model folder:
  model_folder = file.path('models', this_type, 'tinyVAST-multi', n_fac)
  dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------------------------
  # Construc SEM file:
  sem_p1 = NULL
  for(k in 1:n_fac){
    sem_p1 = c(sem_p1, 
      paste(
      rep(paste0('f', k), times = n_sp-k+1),
      (1+k-1):n_sp,
      sep = ' -> '))
  }
  sem_p1 = paste(sem_p1, paste0('o', 1:(n_sp*n_fac - sum(0:(n_fac-1)))), sep = ', ')
  sem_p2 = paste(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' <-> '),
    paste0('sd_o_f', 1:n_fac),
    sep = ', '
  )
  sem_p3 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' <-> '),
    ', NA, 0'
  )
  sem_file = c(sem_p1, sem_p2, sem_p3)
  writeLines(sem_file, file.path(model_folder, paste0("sem_input_sp", n_sp, ".txt")))
  sem_mod = paste0("\n  ", paste(sem_file, collapse = '\n  '), "\n")
  
  
  # -------------------------------------------------------------------------
  # Construc DSEM file:
  dsem_p1 = NULL
  for(k in 1:n_fac){
    dsem_p1 = c(dsem_p1, 
               paste(
                 rep(paste0('f', k), times = n_sp-k+1),
                 (1+k-1):n_sp,
                 sep = ' -> '))
  }
  dsem_p1 = paste(dsem_p1, paste0('e', 1:(n_sp*n_fac - sum(0:(n_fac-1)))), sep = ', 0, ')
  dsem_p2 = paste0(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' -> '),
    ', 1, NA, 0'
  )
  dsem_p3 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' -> '),
    ', 1, NA, 0'
  )
  dsem_p4 = paste(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' <-> '),
    paste0('sd_e_f', 1:n_fac),
    sep = ', 0, '
  )
  dsem_p5 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' <-> '),
    ', 0, NA, 0'
  )
  dsem_file = c(dsem_p1, dsem_p2, dsem_p3, dsem_p4, dsem_p5)
  writeLines(dsem_file, file.path(model_folder, paste0("dsem_input_sp", n_sp, ".txt")))
  dsem_mod = paste0("\n  ", paste(dsem_file, collapse = '\n  '), "\n")
  
  
  # -------------------------------------------------------------------------
  # Fit tinyVAST models
  
  # delta lognormal Poisson link
  jtV_mod_ln = tryCatch(tinyVAST( data = sp_data,
                                     formula = Catch ~ 0 + factor(Year):factor(sp_number),
                                     delta_options = list(
                                       delta_formula = ~ 0 + factor(Year):factor(sp_number),
                                       delta_sem = sem_mod,
                                       delta_dsem = dsem_mod),
                                     variables = c(paste0('f', 1:n_fac), 1:n_sp ),
                                     variable_column = 'sp_number',
                                     space_columns = c("Lon", "Lat"),
                                     time_column = 'Year',
                                     family = delta_lognormal(),
                                     spatial_graph = my_mesh,
                                     control = tinyVASTcontrol(gmrf="proj")), error = function(e) conditionMessage(e))
  # delta lognormal Poisson link
  jtV_mod_ln_pl = tryCatch(tinyVAST( data = sp_data,
                      formula = Catch ~ 0 + factor(Year):factor(sp_number),
                      delta_options = list(
                           delta_formula = ~ 0 + factor(Year):factor(sp_number),
                           delta_sem = sem_mod,
                           delta_dsem = dsem_mod),
                      variables = c(paste0('f', 1:n_fac), 1:n_sp ),
                      variable_column = 'sp_number',
                      space_columns = c("Lon", "Lat"),
                      time_column = 'Year',
                      family = delta_lognormal(type="poisson-link"),
                      spatial_graph = my_mesh,
                      control = tinyVASTcontrol(gmrf="proj")), error = function(e) conditionMessage(e))
  # delta gamma
  jtV_mod_gm = tryCatch(tinyVAST( data = sp_data,
                                     formula = Catch ~ 0 + factor(Year):factor(sp_number),
                                     delta_options = list(
                                       delta_formula = ~ 0 + factor(Year):factor(sp_number),
                                       delta_sem = sem_mod,
                                       delta_dsem = dsem_mod),
                                     variables = c(paste0('f', 1:n_fac), 1:n_sp ),
                                     variable_column = 'sp_number',
                                     space_columns = c("Lon", "Lat"),
                                     time_column = 'Year',
                                     family = delta_gamma(),
                                     spatial_graph = my_mesh,
                                     control = tinyVASTcontrol(gmrf="proj")), error = function(e) conditionMessage(e))
  # delta gamma Poisson link
  jtV_mod_gm_pl = tryCatch(tinyVAST( data = sp_data,
                                     formula = Catch ~ 0 + factor(Year):factor(sp_number),
                                     delta_options = list(
                                       delta_formula = ~ 0 + factor(Year):factor(sp_number),
                                       delta_sem = sem_mod,
                                       delta_dsem = dsem_mod),
                                     variables = c(paste0('f', 1:n_fac), 1:n_sp ),
                                     variable_column = 'sp_number',
                                     space_columns = c("Lon", "Lat"),
                                     time_column = 'Year',
                                     family = delta_gamma(type="poisson-link"),
                                     spatial_graph = my_mesh,
                                     control = tinyVASTcontrol(gmrf="proj")), error = function(e) conditionMessage(e))
  # Tweedie
  jtV_mod_tw = tryCatch(tinyVAST( data = sp_data,
                                     formula = Catch ~ 0 + factor(Year):factor(sp_number),
                                     sem = sem_mod,
                                     dsem = dsem_mod,
                                     variables = c(paste0('f', 1:n_fac), 1:n_sp ),
                                     variable_column = 'sp_number',
                                     space_columns = c("Lon", "Lat"),
                                     time_column = 'Year',
                                     family = tweedie(),
                                     spatial_graph = my_mesh,
                                     control = tinyVASTcontrol(gmrf="proj")), error = function(e) conditionMessage(e))
  
  # Model order:
  mod_order = c('jtV_mod_ln', 'jtV_mod_ln_pl', 'jtV_mod_gm', 'jtV_mod_gm_pl', 'jtV_mod_tw')
  
  # Check AIC:
  aic_vec = numeric(5) # because 5 models
  aic_vec[1] = ifelse(!is.character(jtV_mod_ln), AIC(jtV_mod_ln), NA)
  aic_vec[2] = ifelse(!is.character(jtV_mod_ln_pl), AIC(jtV_mod_ln_pl), NA)
  aic_vec[3] = ifelse(!is.character(jtV_mod_gm), AIC(jtV_mod_gm), NA)
  aic_vec[4] = ifelse(!is.character(jtV_mod_gm_pl), AIC(jtV_mod_gm_pl), NA)
  aic_vec[5] = ifelse(!is.character(jtV_mod_tw), AIC(jtV_mod_tw), NA)
  
  # save aic df:
  aic_df = data.frame(model = mod_order, aic = aic_vec)
  write.csv(aic_df, file = file.path(model_folder, 'AIC.csv'), row.names = FALSE)
  
  # -------------------------------------------------------------------------
  # Select best model based on AIC:
  best_mod = NA
  
  if(length(which.min(aic_vec)) > 0) {

  best_mod = mod_order[which.min(aic_vec)[1]]
  jtVModel = get(best_mod)  
  # Save best model only
  save(jtVModel, file = file.path(model_folder, 'jtVModel.RData'))
  
  # Check residuals
  y_ir = replicate( n = 100, expr = jtVModel$obj$simulate()$y_i )
  res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                              observedResponse = sp_data$Catch, 
                              fittedPredictedResponse = fitted(jtVModel) )
  png(file = file.path(plot_folder, 'resid_check.png'), width = img_width, 
      height = 100, res = img_res, units = 'mm')
  plot(res)
  dev.off()

  # -------------------------------------------------------------------------
  # Plot Omega 2nd component (Factors):
  n_comps = as.vector(jtVModel$tmb_inputs$tmb_data$components_e)
  omega_slot = ifelse(n_comps == 2, 'omega2_sc', 'omega_sc')
  epsilon_slot = ifelse(n_comps == 2, 'epsilon2_stc', 'epsilon_stc')  
  plot_dat = NULL
  for(i in 1:n_fac) {
    tmp = data.frame(Lon = jtVModel$spatial_graph$loc[,1], 
                     Lat = jtVModel$spatial_graph$loc[,2], 
                     omega2 = jtVModel$internal$parlist[[omega_slot]][,i],
                     type_fac = paste0('factor_', i))
    plot_dat = rbind(plot_dat, tmp)
  }
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = omega2), size = 2) + scale_colour_gradient2() + labs(color = 'Omega2')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ type_fac)
  ggsave(filename = paste0('omega2_fac', img_type), path = plot_folder, plot = p1, 
         width = img_width, height = 90, units = 'mm', dpi = img_res)
  
  # -------------------------------------------------------------------------
  # Plot Epsilon 2nd component (Factors):
  all_years = sort(unique(sp_data$Year))
  for(j in 1:n_fac) {
    plot_dat = NULL
    for(i in seq_along(all_years)) {
      tmp_df = data.frame(year = all_years[i],
                          Lon = jtVModel$spatial_graph$loc[,1], 
                          Lat = jtVModel$spatial_graph$loc[,2], 
                          epsilon2 = jtVModel$internal$parlist[[epsilon_slot]][,i,j])
      plot_dat = rbind(plot_dat, tmp_df)
    }
    plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = epsilon2), size = 1) + scale_colour_gradient2() + labs(color = 'Epsilon2')
    p1 = add_sf_map(p1)
    p1 = p1 + facet_wrap(~ year)
    ggsave(filename = paste0('epsilon2_Factor', j, img_type), path = plot_folder, plot = p1, 
           width = img_width, height = 130, units = 'mm', dpi = img_res)
  }

  # -------------------------------------------------------------------------
  # Plot predicted values 
  pred_df = as.data.frame(extraRegion_tinyVAST)
  for(i in 1:n_sp) {
    tmp_pred = pred_df
    this_sp = cumsp_data$sp_name[i]
    tmp_pred$sp_number = i
    tmp_pred$mu_g = predict(jtVModel, newdata = tmp_pred, what = "mu_g")
    
    plot_dat = tmp_pred %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = mu_g), size = 0.5) + scale_colour_viridis() + labs(color = 'mu_g')
    p1 = add_sf_map(p1)
    p1 = p1 + facet_wrap(~ Year)
    ggsave(filename = paste0('mu_g_', this_sp, img_type), path = plot_folder, plot = p1, 
           width = img_width, height = 130, units = 'mm', dpi = img_res)
    
  }
  
  # -------------------------------------------------------------------------
  # Association sp to factors (theta):
  theta_slot = ifelse(n_comps == 2, 'theta2_z', 'theta_z') 
  Lhat_cf = matrix(0, nrow = n_sp, ncol = n_fac)
  theta2_vec = as.list(jtVModel$sdrep, what="Estimate")$theta2_z[1:(n_fac*n_sp - sum(0:(n_fac-1)))]
  Lhat_cf[lower.tri(Lhat_cf, diag=TRUE)] = theta2_vec
  Lhat_cf = rotate_pca( L_tf = Lhat_cf, order = "decreasing" )$L_tf
  dimnames(Lhat_cf) = list( cumsp_data$sp_name[1:n_sp],
                            paste0("Factor ", 1:ncol(Lhat_cf)) )
  load_df = as.data.frame(Lhat_cf)
  load_df$species = rownames(load_df)
  load_df = pivot_longer(load_df, cols = starts_with("Factor"), names_to = 'Factor', values_to = 'loading')
  
  p1 = ggplot(data = load_df) + geom_col(aes(x = species, y = loading), position = 'dodge') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_rect(fill="white")) +
    labs(y = 'Factor loading', x = NULL) +
    facet_wrap(~ Factor, scales = 'free_y')
  ggsave(filename = paste0('Loading_Omega2', img_type), path = plot_folder, plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
  
  # -------------------------------------------------------------------------
  # Association sp to factors (Epsilon2):
  beta_slot = ifelse(n_comps == 2, 'beta2_z', 'beta_z') 
  Lhat_cf = matrix(0, nrow = n_sp, ncol = n_fac)
  beta2_vec = as.list(jtVModel$sdrep, what="Estimate")$beta2_z[1:(n_fac*n_sp - sum(0:(n_fac-1)))]
  Lhat_cf[lower.tri(Lhat_cf, diag=TRUE)] = beta2_vec
  Lhat_cf = rotate_pca( L_tf = Lhat_cf, order = "decreasing" )$L_tf
  dimnames(Lhat_cf) = list( cumsp_data$sp_name[1:n_sp],
                            paste0("Factor ", 1:ncol(Lhat_cf)) )
  load_df = as.data.frame(Lhat_cf)
  load_df$species = rownames(load_df)
  load_df = pivot_longer(load_df, cols = starts_with("Factor"), names_to = 'Factor', values_to = 'loading')
  
  p1 = ggplot(data = load_df) + geom_col(aes(x = species, y = loading), position = 'dodge') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_rect(fill="white")) +
    labs(y = 'Factor loading', x = NULL) +
    facet_wrap(~ Factor, scales = 'free_y')
  ggsave(filename = paste0('Loading_Epsilon2', img_type), path = plot_folder, plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
  
  cat("Nfac finished:", n_fac_vec[m], "-Best model:", best_mod, "\n")
  
  } # conditional
  
} # n_fac loop  
  

# -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
  # Get index:
  Index = NULL
  # To save annual estimates by species:
  plot_dir_sp = file.path(plot_folder, 'Bycatch_est_ts_sp')
  dir.create(plot_dir_sp, showWarnings = FALSE)
  for(i in 1:n_sp) {
    this_sp = cumsp_data$sp_name[i]
    tmp_pred = extraRegion_tinyVAST
    tmp_pred$sp_number = i
    Est = sapply( sort(unique(sp_data$Year)), FUN=\(t) {
      pred_data = subset(tmp_pred, Year == t)
      integrate_output(jtVModel, newdata = pred_data, area = pred_data$n_sets) 
    })
    tmp = data.frame( Year = sort(unique(sp_data$Year)), t(Est))
    colnames(tmp) = c('Year', 'est', 'se', 'est_corr', 'se_corr')
    # Calculate CI:
    tmp$est = tmp$est_corr # replace uncorrected values by corrected?
    tmp$lwr = tmp[,'est'] - 1.96*tmp[,'se']
    tmp$upr = tmp[,'est'] + 1.96*tmp[,'se']
    tmp$species = this_sp
    Index = rbind(Index, tmp)
    
    # Plot annual estimates with observed estimate:
    obs_df = read.csv(file = file.path('data', this_type, 'bycatch_est_obs', paste0(this_sp, '.csv')))
    p1 = plot_time_predictions(tmp, obs_df, Year, est, lwr, upr, yLab = 'Estimated bycatch (tons)')
    ggsave(filename = paste0(this_sp, img_type), path = plot_dir_sp, plot = p1, 
           width = img_width*0.5, height = 70, units = 'mm', dpi = img_res)
    
  }
  Index$model = best_mod
  Index$type = 'joint'
  write.csv(Index, file = file.path(plot_folder, 'Bycatch_est_ts.csv'), row.names = FALSE)

