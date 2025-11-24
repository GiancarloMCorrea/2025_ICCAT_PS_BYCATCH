rm(list = ls())

# Define type of school to be analyzed:
source('code/1_single/load_libs.R')
if(this_type == 'FOB') mesh_cutoff = 1.5
if(this_type == 'FSC') mesh_cutoff = 1
n_cores = 4 # to run in parallel

# -------------------------------------------------------------------------
# Specify full formula:
# Always include two components regardless of the statistical family
full_formula = list(as.formula(bycatch ~ 0 + fyear + quarter + trop_catch),
                    as.formula(bycatch ~ 0 + fyear + quarter + trop_catch)) # for sampling coverage

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path(data_folder, 'weight_data.rds'))
effPoints = readRDS(file = file.path(data_folder, 'effPoints.rds')) # for prediction
load(file.path(data_folder, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Filter species:
mod_data = weight_data
selsp_data = readRDS(file.path(data_folder, 'model_cat_sp.rds'))
# Tmp data, to make mesh. the selection of sp does not matter since all sp have the same obs locations
tmp_data = mod_data %>% dplyr::filter(sp_name %in% selsp_data$sp_name[1])
n_mod_sp = nrow(selsp_data) # IMPORANT!: select based on discussion with Jon!


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Make mesh (will be the same for all sp):
prev_mesh = fm_mesh_2d( tmp_data[,c('lon', 'lat')], cutoff = mesh_cutoff )
p1 = ggplot() + geom_fm(data = prev_mesh)
ggsave(paste0('map_mesh', img_type), path = file.path(plot_folder), 
       plot = p1, width = img_width*0.5, height = 80, units = 'mm', dpi = img_res)

# Plot mesh with observations:
mesh_df = data.frame(lon = prev_mesh$loc[,1], lat = prev_mesh$loc[,2]) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
plot_dat = tmp_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
p1 = ggplot() + 
  geom_sf(data = plot_dat, color = 'gray80', fill = 'gray80', shape = 21, size = 0.5) +
  geom_sf(data = mesh_df, color = 'black', fill = 'black', shape = 21, size = 0.5) 
p1 = add_sf_map(p1)
ggsave(paste0('map_mesh_obs', img_type), path = plot_folder, 
       plot = p1, width = img_width*0.75, height = 120, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Run models in parallel:

source('code/1_single/sdmTMB-config.R')
snowfall::sfInit(parallel=TRUE, cpus=n_cores)
snowfall::sfExportAll()
trash = snowfall::sfLapply(1:n_mod_sp, function(i) {
  
  # Load required libraries:
  require(dplyr)
  require(ggplot2)
  require(sdmTMB)
  require(sf)
  require(fmesher)
  
  # Create dir to save plots and model outputs:
  this_sp = selsp_data$sp_name[i]
  this_model_folder = file.path(model_folder, this_sp)
  this_plot_folder = file.path(plot_folder, "model_sp", this_sp)
  dir.create(this_model_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(this_plot_folder, showWarnings = FALSE, recursive = TRUE)
  this_formula = full_formula
  
  # Filter data:
  sp_data = mod_data %>% dplyr::filter(sp_name %in% this_sp) 
  # For plotting later:
  pres_data = sp_data %>% mutate(presence = as.factor(if_else(bycatch > 0, 1, 0)))
  
  # Check sum catch by year:
  plot_dat = pres_data %>% group_by(year) %>% summarise(bycatch = sum(bycatch), prop = (length(which(presence == 1))/n())*100)
  png(file = file.path(this_plot_folder, 'catch_prop_ts.png'), width = img_width*0.75, 
      height = 150, res = img_res, units = 'mm')
  par(mfrow = c(2,1))
  par(mar = c(3, 4, 0.5, 0.5))
  plot(plot_dat$year, plot_dat$bycatch, type = 'l', xlab = '', ylab = 'Bycatch (t)', ylim = c(0, max(plot_dat$bycatch)*1.05))
  par(mar = c(3, 4, 0.5, 0.5))
  plot(plot_dat$year, plot_dat$prop, type = 'l', xlab = '', ylab = 'Positive sets (%)', ylim = c(0, max(plot_dat$prop)*1.05))
  abline(h = 0, lty = 2)
  dev.off()
  
  # Make maps by sp group (weight):
  pres_data = pres_data %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326, remove = FALSE)
  p1 = ggplot(pres_data) + 
    geom_sf(data = pres_data %>% dplyr::filter(presence == 0), color = 'gray70', fill = 'gray70', 
            shape = 21, size = 0.5, alpha = 0.5) +
    geom_sf(data = pres_data %>% dplyr::filter(presence == 1), color = 'red', fill = 'red', 
            shape = 21, size = 0.5) 
  p1 = add_sf_map(p1)
  p1 = p1 + ggtitle(label = this_sp) + theme(legend.position = 'none') + 
    facet_wrap(~ year)
  ggsave(paste0('obs_map_presence', img_type), path = this_plot_folder, 
         plot = p1, width = img_width, height = 140, units = 'mm', dpi = img_res)
  
  # Remove years with no presence (if any):
  yr_pa = sp_data %>% group_by(year) %>% summarise(bycatch = sum(bycatch))
  yr_keep = yr_pa$year[which(yr_pa$bycatch > 0)]
  sp_data = sp_data %>% filter(year %in% yr_keep)
  if(length(yr_keep) == 1) {
    this_formula = list(update(this_formula[[1]], ~ . - fyear),
                        update(this_formula[[2]], ~ . - fyear))
  }
  
  # Remove quarters with no presence (if any):
  qt_pa = sp_data %>% group_by(quarter) %>% summarise(bycatch = sum(bycatch))
  qt_keep = qt_pa$quarter[which(qt_pa$bycatch > 0)]
  sp_data = sp_data %>% filter(quarter %in% qt_keep)
  if(length(qt_keep) == 1) {
    this_formula = list(update(this_formula[[1]], ~ . - quarter),
                        update(this_formula[[2]], ~ . - quarter))
  }
  
  # Prepare data for model:
  sp_data = sp_data %>% mutate(fyear = factor(year, levels = sort(unique(sp_data$year))))
  
  # Plot distribution of values
  p1 = ggplot(data = sp_data, aes(x = fyear, y = bycatch)) +
    geom_boxplot() + 
    theme_bw() +
    labs(x = 'Year', y = 'Observed bycatch (t)') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  ggsave(paste0('bycatch_obs_val', img_type), path = this_plot_folder, plot = p1,
         width = img_width*0.75, height = 100, units = 'mm', dpi = img_res)
  
  # Make mesh again (with data that will be used in model)
  sp_mesh = sdmTMB::make_mesh(data = sp_data, xy_cols = c('lon', 'lat'),
                              mesh = fm_mesh_2d( sp_data[,c('lon', 'lat')], cutoff = mesh_cutoff ))
  
  # Run model ---------------------------------------------------------------
  model_est = run_sdmTMB_model(sp_data, this_formula, sp_mesh, this_goal, 
                               this_model_folder, this_plot_folder, this_sp,
                               effPoints, yr_keep, qt_keep, MyGrid)

} )

# Stop cluster
snowfall::sfStop()
