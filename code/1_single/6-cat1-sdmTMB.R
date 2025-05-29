
# Run model ---------------------------------------------------------------
n_comps = 2

mod_init <- tryCatch(sdmTMB(
  data = sp_data,
  formula = this_formula,
  mesh = sp_mesh, 
  time = 'year',
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("off", "iid")
), error = function(e) conditionMessage(e))
if(!is.character(mod_init)) {
  # Remove non-significant terms:
  update_formula = remove_terms_sdmTMB(mod_init, n_comps, this_formula)
  mod_upd <- tryCatch(sdmTMB(
    data = sp_data,
    formula = update_formula,
    mesh = sp_mesh, 
    time = 'year',
    family = delta_lognormal(),
    spatial = list("on", "on"),
    spatiotemporal = list("off", "iid")
  ), error = function(e) conditionMessage(e))
}
if(!is.character(mod_upd)) {
  check_df = sanity(mod_upd)
  check_df = data.frame(name = names(unlist(check_df)), check = unlist(check_df))
  write.csv(check_df, file = file.path(this_model_folder, 'mod_check.csv'), row.names = FALSE)
  summ_table = get_summary_sdmTMB(model = mod_upd, n_comp = n_comps, model_label = this_sp)
  write.csv(summ_table, file = file.path(this_model_folder, 'mod_summ.csv'), row.names = FALSE)
  save(mod_upd, file = file.path(this_model_folder, 'mod.RData'))
}


# Plot residuals ----------------------------------------------------------

this_model = mod_upd
plot_residuals(this_model)

# Plot omega --------------------------------------------------------------

p1 = plot_omega(this_model)
ggsave(filename = paste0('omega', img_type), path = this_plot_folder, plot = p1, 
       width = img_width*n_comps*0.5, height = 70, units = 'mm', dpi = img_res)


# Plot epsilon ------------------------------------------------------------

p1 = plot_epsilon(this_model)
ggsave(filename = paste0('epsilon', img_type), path = this_plot_folder, plot = p1, 
       width = img_width, height = 130, units = 'mm', dpi = img_res)

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
write.csv(pred_time, file = file.path(this_model_folder, 'pred_est_time.csv'), row.names = FALSE)

# Plot ratio and sdmTMB estimates:
ratioEst = readRDS(file.path(data_folder, 'ratio_estimates.rds'))
ratioEstSp = ratioEst %>% filter(sp_name %in% this_sp) %>% group_by(year) %>%
  summarise(est_prod = sum(est_prod), est_sets = sum(est_sets))

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
