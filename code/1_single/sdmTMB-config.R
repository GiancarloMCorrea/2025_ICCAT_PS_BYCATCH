
# -------------------------------------------------------------------------
run_sdmTMB_model = function(sp_data, this_formula, sp_mesh, this_goal = 'sample-coverage', 
                     this_model_folder = NULL, this_plot_folder = NULL, this_sp,
                     effPoints, yr_keep, qt_keep, MyGrid = NULL,
                     this_cat = 1, save_model = TRUE, make_plots = TRUE,
                     save_results = TRUE, check_residuals = TRUE) {
  # this_cat is initial model category
  
  check_df = list(all_ok = FALSE) # initial state
  if(this_goal == 'sample-coverage') {
    n_comps = 1 # always 1 for tweedie
    this_family = tweedie()
  }
  if(this_goal == 'estimation') {
    n_comps = 2 # for category 1 and 2
    this_family = delta_lognormal()
  }
  
  # -------------------------------------------------------------------------
  if(this_goal == 'sample-coverage') {
    
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
    
    # Do not update model with significant covariates:
    mod_upd = mod_init
    
  } # Goal: sample coverage
  
  # -------------------------------------------------------------------------
  if(this_goal == 'estimation') {
    
    # Category 1:
    mod_init <- tryCatch(sdmTMB(
      data = sp_data,
      formula = this_formula,
      mesh = sp_mesh, 
      time = 'year',
      family = this_family,
      spatial = list("on", "on"),
      spatiotemporal = list("off", "iid")
    ), error = function(e) conditionMessage(e))
    if(!is.character(mod_init)) check_df = sanity(mod_init, silent = TRUE)
    
    # Check if model passed
    
    if( !check_df$all_ok ) {
      # Category 1 did not pass, run Category 2:
      this_cat = 2
      mod_init <- tryCatch(sdmTMB(
        data = sp_data,
        formula = this_formula,
        mesh = sp_mesh, 
        time = 'year',
        family = this_family,
        spatial = list("on", "on"),
        spatiotemporal = list("off", "off")
      ), error = function(e) conditionMessage(e))
      if(!is.character(mod_init)) check_df = sanity(mod_init, silent = TRUE)
      
      if( !check_df$all_ok ) {
        this_cat = 3
        # Category 2 did not pass, run Category 3:
        mod_init <- tryCatch(sdmTMB(
          data = sp_data,
          formula = this_formula[[1]],
          mesh = sp_mesh, 
          time = 'year',
          family = tweedie(),
          spatial = "on",
          spatiotemporal = "off"
        ), error = function(e) conditionMessage(e))
        if(!is.character(mod_init)) check_df = sanity(mod_init, silent = TRUE)
        if(check_df$all_ok) {
          n_comps = 1 # for tweedie
          this_family = tweedie() # for tweedie
        }
      }# category 3 conditional
      
    } # category 2 conditional
    
    # Run model with updated formula if it passed sanity checks
    # This part will select significant variables and rerun model
    if(check_df$all_ok) {
      # Remove non-significant terms:
      update_formula = remove_terms_sdmTMB(mod_init, n_comps, this_formula)
      mod_upd <- tryCatch(sdmTMB(
        data = sp_data,
        formula = if(n_comps==1) update_formula[[1]] else update_formula,
        mesh = sp_mesh, 
        time = 'year',
        family = this_family,
        spatial = as.list(mod_init$spatial),
        spatiotemporal = as.list(mod_init$spatiotemporal)
      ), error = function(e) conditionMessage(e))
    }
    
  } # Goal: estimation
  
  pred_time = NULL # time predictions
  PredGrid = NULL # predictions eff
  if(check_df$all_ok) {
    # Plot residuals
    check_df = sanity(mod_upd, silent = TRUE)
    check_df2 = data.frame(name = names(unlist(check_df)), check = unlist(check_df))
    check_df2$category = this_cat
    check_df2$n_comps = n_comps
    if(check_residuals) {
      res_check = plot_residuals(mod_upd, plot_dir = this_plot_folder)
      check_df2$unif_test = round(res_check$unif_test, digits = 3)
      check_df2$outl_test = round(res_check$outl_test, digits = 3)
      check_df2$disp_test = round(res_check$disp_test, digits = 3)
    }
    if(save_results) write.csv(check_df2, file = file.path(this_model_folder, 'mod_check.csv'), row.names = FALSE)
    summ_table = get_summary_sdmTMB(model = mod_upd, n_comp = n_comps, model_label = this_sp)
    summ_table$check = check_df$all_ok
    summ_table$category = this_cat
    if(save_results) saveRDS(summ_table, file = file.path(this_model_folder, 'mod_summ.rds'))
    if(save_model) save(mod_upd, file = file.path(this_model_folder, 'mod.RData'))
    
    if(make_plots) {
      # Plot omega --------------------------------------------------------------
      p1 = plot_omega(mod_upd, n_comps = n_comps)
      ggsave(filename = paste0('omega', img_type), path = this_plot_folder, plot = p1, 
             width = img_width*n_comps*0.5, height = 90, units = 'mm', dpi = img_res)
      
      
      # Plot epsilon ------------------------------------------------------------
      p1 = plot_epsilon(mod_upd, n_comps = n_comps)
      for(p in 1:n_comps) {
        ggsave(filename = paste0('epsilon', p, img_type), path = this_plot_folder, plot = p1[[p]], 
               width = img_width, height = 130, units = 'mm', dpi = img_res)
      }
    }
    
    # Predictions -------------------------------------------------------------
    
    # Filter prediction data for years and quarters when observations are available
    this_pred_df = effPoints %>% filter(year %in% yr_keep, quarter %in% qt_keep)
    this_pred_df = this_pred_df %>% mutate(fyear = factor(year, levels = sort(unique(this_pred_df$year))))
    predictions = predict(mod_upd, newdata = this_pred_df, return_tmb_object = TRUE)
    pred_time = get_index(predictions, area = 1, bias_correct = TRUE)
    pred_time = pred_time %>% dplyr::rename(log_se = se) %>% 
      mutate(log_lwr = log_est - 1.96*log_se,
             log_upr = log_est + 1.96*log_se,
             cv = sqrt(exp(log_se^2) - 1)) 
    pred_time$sp_name = this_sp
    pred_time$check = check_df$all_ok
    pred_time$category = this_cat
    if(save_results) write.csv(pred_time, file = file.path(this_model_folder, 'pred_est_time.csv'), row.names = FALSE)
    
    # To plot predictions spatially:
    PredGrid = predictions$data
    if(n_comps == 1) PredGrid$bycatch_est = mod_upd$family$linkinv(PredGrid$est)
    if(n_comps == 2) {
      if(mod_upd$family$link[1] == 'logit' & mod_upd$family$link[2] == 'log') PredGrid$bycatch_est = inv.logit(PredGrid$est1)*exp(PredGrid$est2)
      if(mod_upd$family$link[1] == 'log' & mod_upd$family$link[2] == 'log') PredGrid$bycatch_est = exp(PredGrid$est1)*exp(PredGrid$est2)
    }
    
    if(make_plots) {
      # Plot ratio and sdmTMB estimates:
      ratioEst = calculate_ratio_bycatch(obs_df = sp_data, eff_df = effPoints, type = 'production')
      ratioEst = ratioEst %>% rename(est_prod = est)
      ratioEstSp = ratioEst %>% filter(sp_name %in% this_sp) %>% group_by(year) %>%
        summarise(est_prod = sum(est_prod), .groups = 'drop')
      
      p1 = plot_time_predictions(model_est = pred_time, ratio_est = ratioEstSp, this_sp = this_sp)
      ggsave(paste0('compare_estimates', img_type), path = this_plot_folder, plot = p1,
             width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)
      
      p1 = plot_predictions(plot_data = PredGrid, MyGrid = MyGrid)
      ggsave(filename = paste0('pred_spt_time', img_type), path = this_plot_folder, plot = p1, 
             width = img_width, height = 140, units = 'mm', dpi = img_res)
      
    } # if plot
    
  } # if model converged
  
  return(list(pred_time, PredGrid))
  
}

