# -------------------------------------------------------------------------
# Start loop over species:
for(i in 1:n_sp) {
  
  # Remove model:
  rm(tVModel)
  
  # Create dir to save plots and model outputs:
  this_sp = selsp_data$sp_name[i]
  this_model_folder = file.path(model_folder, 'single', this_sp)
  this_plot_folder = file.path(plot_folder, "single", this_sp)
  dir.create(this_model_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(this_plot_folder, showWarnings = FALSE, recursive = TRUE)

  # Filter data:
  sp_data = mod_data %>% dplyr::filter(sp_name %in% this_sp) %>% mutate(sp_id = 'sp')
  sp_data = droplevels(sp_data)
  
  # Make mesh again (with data that will be used in model)
  sp_mesh = fm_mesh_2d( sp_data[,c('lon', 'lat')], cutoff = mesh_cutoff )
  
  # Fit tinyVAST model
  tVModel = tryCatch(tinyVAST( data = sp_data,
                              formula = bycatch ~ 1,
                              space_term = "sp <-> sp, sd_p1",
                              delta_options = list(
                                  formula = ~ 1,
                                  space_term = "sp <-> sp, sd_p2"),
                              variable_column = 'sp_id',
                              space_columns = c("lon", "lat"),
                              time_column = 'year',
                              family = delta_lognormal(type="poisson-link"),
                              spatial_domain = sp_mesh), 
                    error = function(e) conditionMessage(e))


  if(!is.character(tVModel)) { # passed convergence
    
    # Save model
    save(tVModel, file = file.path(this_model_folder, 'tVModel.RData'))
    
    # Check residuals
    y_ir = simulate(tVModel, nsim=500, type="mle-mvn")
    res = DHARMa::createDHARMa( simulatedResponse = y_ir,
                                observedResponse = sp_data$bycatch,
                                fittedPredictedResponse = fitted(tVModel) )
    png(file = file.path(this_plot_folder, 'resid_check.png'), width = img_width, 
        height = 100, res = img_res, units = 'mm')
    plot(res)
    dev.off()
    
    # -------------------------------------------------------------------------
    # Produce tables for parameters and convergence 
    n_comps = as.vector(tVModel$tmb_inputs$tmb_data$components_e)
    ParHat = tVModel$obj$env$parList()
    SE = as.list( tVModel$sdrep, report=FALSE, what="Std. Error")
    if(n_comps == 1) { parNames = c('alpha_j', 'theta_z', 
                                    'log_sigma', 'log_kappa') }
    if(n_comps == 2) { parNames = c('alpha_j', 'theta_z', 
                                    'alpha2_j', 'theta2_z', 
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
    write.csv(par_df, file = file.path(this_model_folder, 'parameter_estimates.csv'), row.names = FALSE)
    
    # -------------------------------------------------------------------------
    # Plot Omega by component:
    if(n_comps == 1) omega_slot_vec = 'omega_sc'
    if(n_comps == 2) omega_slot_vec = c('omega_sc', 'omega2_sc')
    tmp_df = list()
    for(cp in seq_along(omega_slot_vec)) {
      omega_slot = omega_slot_vec[cp]
      tmp_df[[cp]] = data.frame(Lon = tVModel$spatial_domain$loc[,1], 
                            Lat = tVModel$spatial_domain$loc[,2], 
                            omega = tVModel$internal$parlist[[omega_slot]][,1],
                            comp = paste0('component_', cp))
    }
    plot_dat = bind_rows(tmp_df)
    plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 2) + scale_colour_gradient2() + labs(color = 'Omega')
    p1 = add_sf_map(p1)  
    p1 = p1 + facet_wrap(~comp)
    ggsave(filename = paste0('omega', cp, img_type), path = this_plot_folder, plot = p1, 
           width = img_width*n_comps*0.5, height = 70, units = 'mm', dpi = img_res)
    
    # -------------------------------------------------------------------------
    # Make predictions:
    pred_df = effPoints
    pred_df$sp_id = 'sp'
    pred_df$mu_g = predict(tVModel, newdata = pred_df, what = "mu_g")
    plot_dat = pred_df %>% group_by(ID) %>% summarise(pred = sum(mu_g), .groups = 'drop')
    plot_dat$sp_name = this_sp
    # Save predictions:
    saveRDS(plot_dat, file = file.path(this_model_folder, 'predictions.rds'))
    
    # Plot:
    grid_dat = left_join(MyGrid, plot_dat, by = 'ID')
    p1 = ggplot(grid_dat) + geom_sf(aes(color = pred, fill = pred)) + 
      scale_colour_viridis() + scale_fill_viridis() + 
      guides(color = 'none') + labs(fill = 'Prediction', title = this_sp) 
    p1 = add_sf_map(p1)
    ggsave(filename = paste0('mu_g', img_type), path = this_plot_folder, plot = p1, 
           width = img_width*0.6, height = 90, units = 'mm', dpi = img_res)
    
    cat("Model converged and done for species", i, "-", this_sp, "\n")
    
  } # conditional convergence
  
} # species loop
