# remotes::install_github(repo = "vast-lib/tinyVAST", ref='dev')
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
library(sdmTMB)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
this_type = 'FSC'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'tinyVAST')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'tinyVAST')
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# N sp to build models
n_mod_sp = 35 # model the first N species

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
extraRegion_tinyVAST = readRDS(file = file.path('data', this_type, 'extraRegion_tinyVAST.rds'))

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
cumsp_data = read.csv(file.path('data', this_type, 'cumsp_data_weight.csv'))
if(n_mod_sp > nrow(cumsp_data)) stop('You selected more N species than present in data.')
# Tmp data, to make mesh. the selection of sp does not matter since all sp have the same obs locations
tmp_data = mod_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name[1])

# -------------------------------------------------------------------------
# Make mesh (will be the same for all sp):
my_mesh = fm_mesh_2d( tmp_data[,c('Lon', 'Lat')], cutoff = 1 )
p1 = ggplot() + geom_fm(data = my_mesh)
ggsave(paste0('map_mesh', img_type), path = file.path(plot_folder), 
       plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)

# Plot mesh with observations:
mesh_df = data.frame(Lon = my_mesh$loc[,1], Lat = my_mesh$loc[,2]) %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
plot_dat = tmp_data %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
p1 = ggplot() + 
  geom_sf(data = plot_dat, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
  geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
p1 = add_sf_map(p1)
ggsave(paste0('map_mesh_obs', img_type), path = file.path(plot_folder), 
       plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Start loop over species:
for(i in 1:n_mod_sp) {

  # Create dir to save plots and model outputs:
  this_sp = cumsp_data$sp_name[i]
  dir.create(file.path(model_folder, this_sp), showWarnings = FALSE)
  dir.create(file.path(plot_folder, this_sp), showWarnings = FALSE)
  
  # Filter data:
  sp_data = mod_data %>% dplyr::filter(sp_name %in% this_sp) 
  # For plotting later:
  pres_data = sp_data %>% mutate(presence = as.factor(if_else(Catch > 0, 1, 0)))
  
  # Check sum catch by year:
  plot_dat = pres_data %>% group_by(Year) %>% summarise(Catch = sum(Catch), prop = (length(which(presence == 1))/n())*100)
  tmp_zero = plot_dat %>% arrange(Year)
  zero_years = which(tmp_zero$prop == 0) # to fix alphas for years with zero catch
  png(file = file.path(plot_folder, this_sp, 'catch_prop_ts.png'), width = img_width*0.75, 
      height = 150, res = img_res, units = 'mm')
  par(mfrow = c(2,1))
  par(mar = c(3, 4, 0.5, 0.5))
  plot(plot_dat$Year, plot_dat$Catch, type = 'l', xlab = '', ylab = 'Catch (t)', ylim = c(0, max(plot_dat$Catch)*1.05))
  par(mar = c(3, 4, 0.5, 0.5))
  plot(plot_dat$Year, plot_dat$prop, type = 'l', xlab = '', ylab = 'Positive sets (%)', ylim = c(0, max(plot_dat$prop)*1.05))
  abline(h = 0, lty = 2)
  dev.off()
  
  # Make maps by sp group (weight):
  pres_data = pres_data %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  
  p1 = ggplot(pres_data) + 
          geom_sf(data = pres_data %>% dplyr::filter(presence == 0), color = 'gray70', fill = 'gray70', 
                  shape = 21, size = 0.5, alpha = 0.5) +
          geom_sf(data = pres_data %>% dplyr::filter(presence == 1), color = 'red', fill = 'red', 
                  shape = 21, size = 0.5) 
  p1 = add_sf_map(p1)
  p1 = p1 + ggtitle(label = this_sp) + theme(legend.position = 'none') + 
    facet_wrap(~ Year)
  ggsave(paste0('obs_map_presence', img_type), path = file.path(plot_folder, this_sp), 
         plot = p1, width = img_width, height = 140, units = 'mm', dpi = img_res)

  # Prepare data for model:
  sp_data = sp_data %>% mutate(sp_id = 'sp')
  sp_data = sp_data %>% mutate(fyear = factor(Year))
  sp_data = as.data.frame(sp_data)
  
  # Fit tinyVAST model
  tV_mod_ln = tryCatch(tinyVAST( data = sp_data,
                           formula = Catch ~ 0 + fyear,
                           space_term = "sp <-> sp, sd_p1",
                           spacetime_term = "sp -> sp, 1, NA, 0",
                           delta_options = list(
                             formula = ~ 0 + fyear,
                             space_term = "sp <-> sp, sd_p2",
                             spacetime_term = "sp -> sp, 1, NA, 0"), # fix rho parameter at zero, so epsilon IID
                           variable_column = 'sp_id',
                           space_columns = c("Lon", "Lat"),
                           time_column = 'Year',
                           family = delta_lognormal(),
                           spatial_domain = my_mesh), error = function(e) conditionMessage(e))
  tV_mod_ln_pl = tryCatch(tinyVAST( data = sp_data,
                      formula = Catch ~ 0 + fyear,
                      space_term = "sp <-> sp, sd_p1",
                      spacetime_term = "sp -> sp, 1, NA, 0",
                      delta_options = list(
                        formula = ~ 0 + fyear,
                        space_term = "sp <-> sp, sd_p2",
                        spacetime_term = "sp -> sp, 1, NA, 0"), # fix rho parameter at zero, so epsilon IID
                      variable_column = 'sp_id',
                      space_columns = c("Lon", "Lat"),
                      time_column = 'Year',
                      family = delta_lognormal(type="poisson-link"),
                      spatial_domain = my_mesh), error = function(e) conditionMessage(e))
  tV_mod_gm = tryCatch(tinyVAST( data = sp_data,
                       formula = Catch ~ 0 + fyear,
                       space_term = "sp <-> sp, sd_p1",
                       spacetime_term = "sp -> sp, 1, NA, 0",
                       delta_options = list(
                         formula = ~ 0 + fyear,
                         space_term = "sp <-> sp, sd_p2",
                         spacetime_term = "sp -> sp, 1, NA, 0"), # fix rho parameter at zero, so epsilon IID
                       variable_column = 'sp_id',
                       space_columns = c("Lon", "Lat"),
                       time_column = 'Year',
                       family = delta_gamma(),
                       spatial_domain = my_mesh), error = function(e) conditionMessage(e))
  tV_mod_gm_pl = tryCatch(tinyVAST( data = sp_data,
                        formula = Catch ~ 0 + fyear,
                        space_term = "sp <-> sp, sd_p1",
                        spacetime_term = "sp -> sp, 1, NA, 0",
                        delta_options = list(
                          formula = ~ 0 + fyear,
                          space_term = "sp <-> sp, sd_p2",
                          spacetime_term = "sp -> sp, 1, NA, 0"), # fix rho parameter at zero, so epsilon IID
                        variable_column = 'sp_id',
                        space_columns = c("Lon", "Lat"),
                        time_column = 'Year',
                        family = delta_gamma(type = "poisson-link"),
                        spatial_domain = my_mesh), error = function(e) conditionMessage(e))
  tV_mod_tw = tryCatch(tinyVAST(data = sp_data,
                      formula = Catch ~ 0 + fyear,
                      space_term = "sp <-> sp, sd_p",
                      spacetime_term = "sp -> sp, 1, NA, 0", # fix rho parameter at zero, so epsilon IID
                      variable_column = 'sp_id',
                      space_columns = c("Lon", "Lat"),
                      time_column = 'Year',
                      family = tweedie(link = 'log'),
                      spatial_domain = my_mesh), error = function(e) conditionMessage(e))

  # Model order:
  mod_order = c('tV_mod_ln', 'tV_mod_ln_pl', 'tV_mod_gm', 'tV_mod_gm_pl', 'tV_mod_tw')
  
  # Check sd reports:
  # tV_mod_ln$sdrep
  # tV_mod_ln_pl$sdrep
  # tV_mod_gm$sdrep
  # tV_mod_gm_pl$sdrep
  # tV_mod_tw$sdrep
  
  # Check AIC:
  aic_vec = numeric(5) # because 5 models
  aic_vec[1] = ifelse(!is.character(tV_mod_ln), AIC(tV_mod_ln), NA)
  aic_vec[2] = ifelse(!is.character(tV_mod_ln_pl), AIC(tV_mod_ln_pl), NA)
  aic_vec[3] = ifelse(!is.character(tV_mod_gm), AIC(tV_mod_gm), NA)
  aic_vec[4] = ifelse(!is.character(tV_mod_gm_pl), AIC(tV_mod_gm_pl), NA)
  aic_vec[5] = ifelse(!is.character(tV_mod_tw), AIC(tV_mod_tw), NA)

  # Add another way to compare models here
  
  # save aic df:
  aic_df = data.frame(model = mod_order, aic = aic_vec)
  write.csv(aic_df, file = file.path(model_folder, this_sp, 'AIC.csv'), row.names = FALSE)

  # -------------------------------------------------------------------------
  # Select best model based on AIC:
  best_mod = NA
  
  if(length(which.min(aic_vec)) > 0) { # if minimum exists
  
  best_mod = mod_order[which.min(aic_vec)[1]]
  tVModel = get(best_mod)
  # Save best model only to save storage space
  save(tVModel, file = file.path(model_folder, this_sp, 'tVModel.RData'))
  
  # Check residuals
  y_ir = replicate( n = 100, expr = tVModel$obj$simulate()$y_i )
  res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                              observedResponse = sp_data$Catch, 
                              fittedPredictedResponse = fitted(tVModel) )
  png(file = file.path(plot_folder, this_sp, 'resid_check.png'), width = img_width, 
      height = 100, res = img_res, units = 'mm')
  plot(res)
  dev.off()
  

  # -------------------------------------------------------------------------
  # Produce tables for parameters and convergence 
  n_comps = as.vector(tVModel$tmb_inputs$tmb_data$components_e)
  ParHat = tVModel$obj$env$parList()
  SE = as.list( tVModel$sdrep, report=FALSE, what="Std. Error")
  if(n_comps == 1) { parNames = c('alpha_j', 'beta_z', 'theta_z', 
                                  'log_sigma', 'log_kappa') }
  if(n_comps == 2) { parNames = c('alpha_j', 'beta_z', 'theta_z', 
                               'alpha2_j', 'beta2_z', 'theta2_z', 
                               'log_sigma', 'log_kappa') }
  # Get estimates :
  save_coeff = list()
  for(k in seq_along(parNames)) {
    tmp_df = data.frame(parameter = names(ParHat[parNames[k]]), 
                        estimate = ParHat[parNames[k]][[1]],
                        se = SE[parNames[k]][[1]])
    tmp_df = tmp_df %>% mutate(zvalue = estimate/se,
                               pvalue = pnorm(-abs(zvalue))*2)
    save_coeff[[k]] = tmp_df
  }
  par_df = bind_rows(save_coeff)
  par_df$max_grad = max(abs(tVModel$obj$gr()))
  par_df$convergence = tVModel$opt$convergence
  write.csv(par_df, file = file.path(model_folder, this_sp, 'parameter_estimates.csv'), row.names = FALSE)
  
  # -------------------------------------------------------------------------
  # Plot Omega by component:
  if(n_comps == 1) omega_slot_vec = 'omega_sc'
  if(n_comps == 2) omega_slot_vec = c('omega_sc', 'omega2_sc')
  for(cp in seq_along(omega_slot_vec)) {
    omega_slot = omega_slot_vec[cp]
    plot_dat = data.frame(Lon = tVModel$spatial_domain$loc[,1], 
                          Lat = tVModel$spatial_domain$loc[,2], 
                          omega = tVModel$internal$parlist[[omega_slot]][,1])
    plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 2) + scale_colour_gradient2() + labs(color = 'Omega')
    p1 = add_sf_map(p1)  
    ggsave(filename = paste0('omega', cp, img_type), path = file.path(plot_folder, this_sp), plot = p1, 
           width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)
  }
  
  # -------------------------------------------------------------------------
  # Plot Epsilon by component:
  if(n_comps == 1) epsilon_slot_vec = 'epsilon_stc'
  if(n_comps == 2) epsilon_slot_vec = c('epsilon_stc', 'epsilon2_stc')
  all_years = sort(unique(tVModel$data$Year))
  for(cp in seq_along(epsilon_slot_vec)) {
    epsilon_slot = epsilon_slot_vec[cp]
    tmp_df = list()
    for(j in seq_along(all_years)) {
      tmp_df[[j]] = data.frame(year = all_years[j],
                               Lon = tVModel$spatial_domain$loc[,1], 
                               Lat = tVModel$spatial_domain$loc[,2], 
                               epsilon = tVModel$internal$parlist[[epsilon_slot]][,j,1])
    }
    plot_dat = bind_rows(tmp_df)
    plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = epsilon), size = 1) + scale_colour_gradient2() + labs(color = 'Epsilon')
    p1 = add_sf_map(p1)
    p1 = p1 + facet_wrap(~ year)
    ggsave(filename = paste0('epsilon', cp, img_type), path = file.path(plot_folder, this_sp), plot = p1, 
           width = img_width, height = 130, units = 'mm', dpi = img_res)
  }
  
  # -------------------------------------------------------------------------
  # Plot predicted values
  pred_df = as.data.frame(extraRegion_tinyVAST)
  pred_df$sp_id = 'sp'
  pred_df$fyear = factor(pred_df$Year)
  pred_df$mu_g = predict(tVModel, newdata = pred_df, what = "mu_g")
  plot_dat = pred_df %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  
  p1 = ggplot(plot_dat) + geom_sf(aes(color = mu_g), size = 0.5) + scale_colour_viridis() + labs(color = 'mu_g')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ Year)
  ggsave(filename = paste0('mu_g', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
  
  } # conditional convergence
  
  cat("Model finished:", i, "-", this_sp, "-Best model:", best_mod, "\n")
  
} # modelling loop


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Make predictions and get index:
for(i in 1:n_mod_sp) {
  
  # Select sp:
  this_sp = cumsp_data$sp_name[i]
  
  # Remove objects just in case
  remove(list = c('tVModel')) # remove objects just in case
  
  # Check if model converged
  mod_exist = file.exists(file.path(model_folder, this_sp, 'tVModel.RData'))
  
  if(mod_exist) {
  # Load best model:
  load(file.path(model_folder, this_sp, 'tVModel.RData'))
  
  # Get index:
  all_years = sort(unique(tVModel$data$Year))
  Est = sapply( all_years, FUN=\(t) {
    pred_data = subset(as.data.frame(extraRegion_tinyVAST), Year == t)
    pred_data$sp_id = 'sp'
    pred_data$fyear = factor(pred_data$Year)
    integrate_output(tVModel, newdata = pred_data, area = pred_data$n_sets) 
  })
  Index = data.frame( Year = all_years, t(Est))
  colnames(Index) = c('Year', 'est', 'se', 'est_corr', 'se_corr')
  Index$est = Index$est_corr # replace uncorrected values by corrected?
  Index$lwr = Index[,'est'] - 1.96*Index[,'se']
  Index$upr = Index[,'est'] + 1.96*Index[,'se']
  Index$type = 'single'
  Index$species = this_sp
  write.csv(Index, file = file.path(plot_folder, this_sp, 'Bycatch_est_ts.csv'), row.names = FALSE)
  
  # Plot annual estimates with observed estimate:
  obs_df = read.csv(file = file.path('data', this_type, 'bycatch_est_obs', paste0(this_sp, '.csv')))
  p1 = plot_time_predictions(Index, obs_df, Year, est, lwr, upr, yLab = 'Estimated bycatch (tons)')
  ggsave(filename = paste0('Bycatch_est_ts', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width*0.5, height = 70, units = 'mm', dpi = img_res)
  
  cat("Index done for species:", i, "-", this_sp, "\n")
  
  } else {
    cat("No model found for species:", i, "-", this_sp, "\n")
  }
  
} # sp loop
