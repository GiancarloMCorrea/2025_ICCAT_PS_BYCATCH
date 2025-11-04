
# Run model ---------------------------------------------------------------
check_df = list(all_ok = FALSE) # initial state
this_cat = 1 # initial category
n_comps = 2 # for category 1 and 2
this_family = delta_lognormal() # for category 1 and 2

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
  if(!is.character(mod_upd)) {
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
}

# Now plot outputs and prediction if model passed:
if(check_df$all_ok) {
  
  # Plot residuals ----------------------------------------------------------
  this_model = mod_upd
  plot_residuals(this_model)
  
  # Plot omega --------------------------------------------------------------
  
  p1 = plot_omega(this_model)
  ggsave(filename = paste0('omega', img_type), path = this_plot_folder, plot = p1, 
         width = img_width*n_comps*0.5, height = 70, units = 'mm', dpi = img_res)
  
  
  # Plot epsilon ------------------------------------------------------------
  
  p1 = plot_epsilon(this_model)
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
  if(n_comps == 1) PredGrid$bycatch_est = this_model$family$linkinv(PredGrid$est)
  if(n_comps == 2) {
    if(this_model$family$link[1] == 'logit' & this_model$family$link[2] == 'log') PredGrid$bycatch_est = inv.logit(PredGrid$est1)*exp(PredGrid$est2)
    if(this_model$family$link[1] == 'log' & this_model$family$link[2] == 'log') PredGrid$bycatch_est = exp(PredGrid$est1)*exp(PredGrid$est2)
  }
  
  p1 = plot_predictions(plot_data = PredGrid)
  ggsave(filename = paste0('pred_spt_time', img_type), path = this_plot_folder, plot = p1, 
         width = img_width, height = 140, units = 'mm', dpi = img_res)

} # if model passed
