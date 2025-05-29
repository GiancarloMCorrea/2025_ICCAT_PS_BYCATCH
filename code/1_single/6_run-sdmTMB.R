# -------------------------------------------------------------------------
# Start loop over species:
for(i in 1:n_mod_sp) {
  
  # Remove models:
  rm(mod_init, mod_upd)
  
  # Create dir to save plots and model outputs:
  this_sp = selsp_data$sp_name[i]
  this_cat = selsp_data$category[i]
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
    labs(x = 'Year', y = 'Observed bycatch (t)') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  ggsave(paste0('bycatch_obs_val', img_type), path = this_plot_folder, plot = p1,
        width = img_width*0.75, height = 100, units = 'mm', dpi = img_res)
  
  # Make mesh again (with data that will be used in model)
  sp_mesh = sdmTMB::make_mesh(data = sp_data, xy_cols = c('lon', 'lat'),
                              mesh = fm_mesh_2d( sp_data[,c('lon', 'lat')], cutoff = mesh_cutoff ))
  
  # Run models:
  
  # Category 1
  if(this_cat == 1) { 
    source("code/1_single/6-cat1-sdmTMB.R")
  }
  
  # Category 2
  if(this_cat == 2) { 
    source("code/1_single/6-cat2-sdmTMB.R")
  }
  
  # Category 3
  # TODO: remove year or quarter if only one level
  if(this_cat == 3) { 
    # Decide what to do here
  }
  
  cat("Model for species", i, this_sp, "completed.\n")
  
} # sp loop
