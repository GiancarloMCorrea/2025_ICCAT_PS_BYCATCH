rm(list = ls())

# Define type of school to be analyzed:
source('code/1_single/load_libs.R')
mesh_cutoff = 1.5
n_cores = 3 # to run in parallel

# -------------------------------------------------------------------------
# Specify full formula:
# full_formula = list(as.formula(bycatch ~ 0 + fyear + quarter + trop_catch + sst),
#                     as.formula(bycatch ~ 0 + fyear + quarter + trop_catch + sst)) # for bycatch estimation
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
n_mod_sp = 26 # IMPORANT!: select based on discussion with Jon!


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
# Run in parallel:

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
    this_formula = list(update(full_formula[[1]], ~ . - fyear),
                        update(full_formula[[2]], ~ . - fyear))
  }
  
  # Remove quarters with no presence (if any):
  qt_pa = sp_data %>% group_by(quarter) %>% summarise(bycatch = sum(bycatch))
  qt_keep = qt_pa$quarter[which(qt_pa$bycatch > 0)]
  sp_data = sp_data %>% filter(quarter %in% qt_keep)
  if(length(qt_keep) == 1) {
    this_formula = list(update(full_formula[[1]], ~ . - quarter),
                        update(full_formula[[2]], ~ . - quarter))
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
  check_df = list(all_ok = FALSE) # initial state
  this_cat = 1 # initial category
  n_comps = 1 # for category 1 and 2
  this_family = tweedie() # for category 1 and 2
  
  # Category 1:
  mod_init <- tryCatch(sdmTMB(
    data = sp_data,
    formula = this_formula[[1]],
    mesh = sp_mesh, 
    time = 'year',
    family = this_family,
    spatial = "on",
    spatiotemporal = "iid"
  ), error = function(e) conditionMessage(e))
  if(!is.character(mod_init)) check_df = sanity(mod_init, silent = TRUE)
  
  # Check if model passed
  
  if( !check_df$all_ok ) {
    # Category 1 did not pass, run Category 2:
    this_cat = 2
    mod_init <- tryCatch(sdmTMB(
      data = sp_data,
      formula = this_formula[[1]],
      mesh = sp_mesh, 
      time = 'year',
      family = this_family,
      spatial = "on",
      spatiotemporal = "off"
    ), error = function(e) conditionMessage(e))
    if(!is.character(mod_init)) check_df = sanity(mod_init, silent = TRUE)
    
  } # category 2 conditional
  
  # Run model with updated formula if it passed sanity checks
  if(check_df$all_ok) {
    mod_upd = mod_init
    check_df = sanity(mod_upd, silent = TRUE)
    check_df2 = data.frame(name = names(unlist(check_df)), check = unlist(check_df))
    check_df2$category = this_cat
    write.csv(check_df2, file = file.path(this_model_folder, 'mod_check.csv'), row.names = FALSE)
    summ_table = get_summary_sdmTMB(model = mod_upd, n_comp = n_comps, model_label = this_sp)
    summ_table$check = check_df$all_ok
    summ_table$category = this_cat
    saveRDS(summ_table, file = file.path(this_model_folder, 'mod_summ.rds'))
    save(mod_upd, file = file.path(this_model_folder, 'mod.RData'))
  }
  
  # Now plot outputs and prediction if model passed:
  if(check_df$all_ok) {
    
    # Plot residuals ----------------------------------------------------------
    this_model = mod_upd
    plot_residuals(this_model, plot_dir = this_plot_folder)
    
    # Plot omega --------------------------------------------------------------
    
    p1 = plot_omega(this_model, n_comps = n_comps)
    ggsave(filename = paste0('omega', img_type), path = this_plot_folder, plot = p1, 
           width = img_width*0.7, height = 90, units = 'mm', dpi = img_res)
    
    
    # Plot epsilon ------------------------------------------------------------
    
    p1 = plot_epsilon(this_model, n_comps = n_comps)
    for(p in 1:n_comps) {
      ggsave(filename = paste0('epsilon', p, img_type), path = this_plot_folder, plot = p1[[p]], 
             width = img_width, height = 130, units = 'mm', dpi = img_res)
    }
    
    
    # Predictions -------------------------------------------------------------
    
    # Filter prediction data for years and quarters when observations are available
    this_pred_df = effPoints %>% filter(year %in% yr_keep, quarter %in% qt_keep)
    this_pred_df = this_pred_df %>% mutate(fyear = factor(year, levels = sort(unique(this_pred_df$year))))
    predictions = predict(this_model, newdata = this_pred_df, return_tmb_object = TRUE)
    pred_time = get_index(predictions, area = 1, bias_correct = TRUE)
    pred_time = pred_time %>% dplyr::rename(log_se = se) %>% 
      mutate(log_lwr = log_est - 1.96*log_se,
             log_upr = log_est + 1.96*log_se,
             cv = sqrt(exp(log_se^2) - 1)) 
    pred_time$sp_name = this_sp
    pred_time$check = check_df$all_ok
    pred_time$category = this_cat
    write.csv(pred_time, file = file.path(this_model_folder, 'pred_est_time.csv'), row.names = FALSE)
    
    # Plot ratio and sdmTMB estimates:
    ratioEst = readRDS(file.path(data_folder, 'ratio_estimates.rds'))
    ratioEstSp = ratioEst %>% filter(sp_name %in% this_sp) %>% group_by(year) %>%
      summarise(est_prod = sum(est_prod), est_sets = sum(est_sets), .groups = 'drop')
    
    p1 = plot_time_predictions(model_est = pred_time, ratio_est = ratioEstSp)
    ggsave(paste0('compare_estimates', img_type), path = this_plot_folder, plot = p1,
           width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)
    
    # -------------------------------------------------------------------------
    # Plot predictions spatially:
    PredGrid = predictions$data
    PredGrid$bycatch_est = this_model$family$linkinv(PredGrid$est)
    
    p1 = plot_predictions(plot_data = PredGrid)
    ggsave(filename = paste0('pred_spt_time', img_type), path = this_plot_folder, plot = p1, 
           width = img_width, height = 140, units = 'mm', dpi = img_res)
  
  } # if model passed

} )

# Stop cluster
snowfall::sfStop()
