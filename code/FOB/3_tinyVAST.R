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
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type, 'tinyVAST')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Define model folder:
model_folder = file.path('models', this_type, 'tinyVAST')
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path('data', this_type, 'weight_data.rds'))
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
  # Check sum catch by year:
  plot_dat = sp_data %>% group_by(Year) %>% summarise(Catch = sum(Catch), prop = (length(which(Catch > 0))/n())*100)
  png(file = file.path(plot_folder, this_sp, 'catch_prop_ts.png'), width = img_width*0.75, 
      height = 150, res = img_res, units = 'mm')
  par(mfrow = c(2,1))
  par(mar = c(3, 4, 0.5, 0.5))
  plot(plot_dat$Year, plot_dat$Catch, type = 'l', xlab = '', ylab = 'Catch (t)', ylim = c(0, max(plot_dat$Catch)*1.05))
  par(mar = c(3, 4, 0.5, 0.5))
  plot(plot_dat$Year, plot_dat$prop, type = 'l', xlab = '', ylab = 'Positive sets (%)', ylim = c(0, max(plot_dat$prop)*1.05))
  dev.off()
  
  # Make maps by sp group (weight):
  pres_data = sp_data %>% mutate(presence = as.factor(if_else(Catch > 0, 1, 0)))
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
  sp_data = as.data.frame(sp_data)
  
  # Make mesh
  my_mesh = fm_mesh_2d( sp_data[,c('Lon', 'Lat')], cutoff = 2 )
  p1 = ggplot() + geom_fm(data = my_mesh)
  ggsave(paste0('map_mesh', img_type), path = file.path(plot_folder, this_sp), 
         plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)
  
  # Plot mesh with observations:
  mesh_df = data.frame(Lon = my_mesh$loc[,1], Lat = my_mesh$loc[,2]) %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot() + 
    geom_sf(data = pres_data, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
    geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
  p1 = add_sf_map(p1)
  ggsave(paste0('map_mesh_obs', img_type), path = file.path(plot_folder, this_sp), 
         plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)

  # Fit tinyVAST model
  tVModel = tinyVAST( data = sp_data,
                      formula = Catch ~ 0 + factor(Year),
                      delta_options = list(
                           delta_formula = ~ 0 + factor(Year),
                           delta_sem = "sp <-> sp, spatial_sd",
                           delta_dsem = "sp -> sp, 1, NA, 0"), # fix rho parameter at zero, so epsilon IID
                      variable_column = 'sp_id',
                      space_columns = c("Lon", "Lat"),
                      time_column = 'Year',
                      family = delta_lognormal(type="poisson-link"),
                      spatial_graph = my_mesh )
  tVModel$opt$convergence
  save(tVModel, file = file.path(model_folder, this_sp, 'fit.RData'))
  
  # Check residuals
  y_ir = replicate( n = 100, expr = tVModel$obj$simulate()$y_i )
  res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                              observedResponse = sp_data$Catch, 
                              fittedPredictedResponse = fitted(tVModel) )
  png(file = file.path(plot_folder, this_sp, 'resid_check.png'), width = img_width, 
      height = 100, res = img_res, units = 'mm')
  plot(res)
  dev.off()
  
  # Get index:
  all_years = sort(unique(sp_data$Year))
  Est = sapply( all_years, FUN=\(t) {
    pred_data = subset(as.data.frame(extraRegion_tinyVAST), Year == t)
    pred_data$sp_id = 'sp'
    integrate_output(tVModel, newdata = pred_data, area = pred_data$n_sets) 
    })
  Index = data.frame( Year = all_years, t(Est))
  colnames(Index) = c('Year', 'est', 'se', 'est_corr', 'se_corr')
  Index$est = Index$est_corr # replace uncorrected values by corrected?
  Index$lwr = Index[,'est'] - 1.96*Index[,'se']
  Index$upr = Index[,'est'] + 1.96*Index[,'se']
  Index$model = 'tinyVAST'
  Index$type = 'single'
  Index$species = this_sp
  write.csv(Index, file = file.path(plot_folder, this_sp, 'Bycatch_est_ts.csv'), row.names = FALSE)
  
  # Plot annual estimates with observed estimate:
  obs_df = read.csv(file = file.path('data', this_type, 'bycatch_est_obs', paste0(this_sp, '.csv')))
  p1 = plot_time_predictions(Index, obs_df, Year, est, lwr, upr, yLab = 'Estimated bycatch (tons)')
  ggsave(filename = paste0('Bycatch_est_ts', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width*0.5, height = 70, units = 'mm', dpi = img_res)

  # -------------------------------------------------------------------------
  # Plot Omega 2nd component:
  plot_dat = data.frame(Lon = tVModel$spatial_graph$loc[,1], 
                        Lat = tVModel$spatial_graph$loc[,2], 
                        omega2 = tVModel$internal$parlist$omega2_sc[,1])
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = omega2), size = 2) + scale_colour_gradient2() + labs(color = 'Omega2')
  p1 = add_sf_map(p1)  
  ggsave(filename = paste0('omega2', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)
  
  
  # -------------------------------------------------------------------------
  # Plot Epsilon 2nd component:
  tmp_df = list()
  for(i in seq_along(all_years)) {
    tmp_df[[i]] = data.frame(year = all_years[i],
                             Lon = tVModel$spatial_graph$loc[,1], 
                             Lat = tVModel$spatial_graph$loc[,2], 
                             epsilon2 = tVModel$internal$parlist$epsilon2_stc[,i,1])
  }
  plot_dat = bind_rows(tmp_df)
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = epsilon2), size = 1) + scale_colour_gradient2() + labs(color = 'Epsilon2')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ year)
  ggsave(filename = paste0('epsilon2', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
  

  # -------------------------------------------------------------------------
  # Plot predicted values
  pred_df = as.data.frame(extraRegion_tinyVAST)
  pred_df$sp_id = 'sp'
  pred_df$mu_g = predict(tVModel, newdata = pred_df, what = "mu_g")
  plot_dat = pred_df %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  
  p1 = ggplot(plot_dat) + geom_sf(aes(color = mu_g), size = 0.5) + scale_colour_viridis() + labs(color = 'mu_g')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~ Year)
  ggsave(filename = paste0('mu_g', img_type), path = file.path(plot_folder, this_sp), plot = p1, 
         width = img_width, height = 130, units = 'mm', dpi = img_res)
  
} # modelling loop
