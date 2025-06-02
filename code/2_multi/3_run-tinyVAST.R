
# Loop over number of factors
for(m in seq_along(n_fac_vec)) {
  
  n_fac = n_fac_vec[m]
  # Define folder to save results
  this_plot_folder = file.path(plot_folder, n_fac)
  dir.create(this_plot_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Define model folder:
  this_model_folder = file.path(model_folder, n_fac)
  dir.create(this_model_folder, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------------------------
  # Construc space_term file:
  s_input = make_space_term(n_sp = n_sp, n_fac = n_fac, par_lab = 'o', save_folder = this_model_folder)
  
  # -------------------------------------------------------------------------
  # Fit tinyVAST models
  
  # delta lognormal Poisson link
  # this_mod = tryCatch(tinyVAST( data = mod_data,
  #                     formula = bycatch ~ 0 + fsp,
  #                     space_term = s_input,
  #                     delta_options = list(
  #                          formula = ~ 0 + fsp,
  #                          space_term = s_input
  #                          ),
  #                     variables = c(paste0('f', 1:n_fac), 1:n_sp ),
  #                     variable_column = 'sp_number',
  #                     space_columns = c("lon", "lat"),
  #                     time_column = 'year',
  #                     family = delta_lognormal(type="poisson-link"),
  #                     spatial_domain = my_mesh,
  #                     control = tinyVASTcontrol(gmrf="proj")), 
  #                     error = function(e) conditionMessage(e))
  
  # Tweedie:
  jtVModel = tryCatch(tinyVAST( data = mod_data,
                                formula = bycatch ~ 0 + fsp,
                                space_term = s_input,
                                variables = c(paste0('f', 1:n_fac), 1:n_sp ),
                                variable_column = 'sp_number',
                                space_columns = c("lon", "lat"),
                                time_column = 'year',
                                family = tweedie(),
                                spatial_domain = my_mesh,
                                control = tinyVASTcontrol(gmrf="proj")), 
                      error = function(e) conditionMessage(e))

  # Save model
  save(jtVModel, file = file.path(this_model_folder, 'jtVModel.RData'))
  
  if(!is.character(jtVModel)) {
  
  # Check residuals
  y_ir = replicate( n = 500, expr = jtVModel$obj$simulate()$y_i )
  res = DHARMa::createDHARMa( simulatedResponse = y_ir, 
                              observedResponse = mod_data$bycatch, 
                              fittedPredictedResponse = fitted(jtVModel) )
  png(file = file.path(this_plot_folder, 'resid_check.png'), width = img_width, 
      height = 100, res = img_res, units = 'mm')
  plot(res)
  dev.off()

  # -------------------------------------------------------------------------
  # Produce tables for parameters and convergence 
  ParHat = jtVModel$obj$env$parList()
  SE = as.list( jtVModel$sdrep, report=FALSE, what="Std. Error")
  parNames = c('alpha_j', 'theta_z', 'log_sigma', 'log_kappa') 
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
  par_df$max_grad = max(abs(jtVModel$obj$gr()))
  par_df$convergence = jtVModel$opt$convergence
  write.csv(par_df, file = file.path(this_model_folder, 'parameter_estimates.csv'), row.names = FALSE)
  
  # -------------------------------------------------------------------------
  # Plot Omega (Factors):
  if(n_comps == 1) omega_slot_vec = 'omega_sc'
  if(n_comps == 2) omega_slot_vec = c('omega_sc', 'omega2_sc')
  for(cp in seq_along(omega_slot_vec)) {
    plot_dat = NULL
    omega_slot = omega_slot_vec[cp]
    for(i in 1:n_fac) {
      tmp = data.frame(Lon = jtVModel$spatial_domain$loc[,1], 
                       Lat = jtVModel$spatial_domain$loc[,2], 
                       omega = jtVModel$internal$parlist[[omega_slot]][,i],
                       type_fac = paste0('factor_', i))
      plot_dat = rbind(plot_dat, tmp)
    }
    plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 2) + scale_colour_gradient2() + labs(color = 'Omega')
    p1 = add_sf_map(p1)
    p1 = p1 + facet_wrap(~ type_fac)
    ggsave(filename = paste0('omega', cp, img_type), path = this_plot_folder, plot = p1, 
           width = img_width, height = 90, units = 'mm', dpi = img_res)
  }

  # -------------------------------------------------------------------------
  # Association sp to factors (time_term, nu): COMBINED 
  # Vhat_total = save_Lnu_mat[[1]] %*% t(save_Lnu_mat[[1]]) + save_Lnu_mat[[2]] %*% t(save_Lnu_mat[[2]])
  # Vhat_total = t(save_Lnu_mat[[1]]) %*% save_Lnu_mat[[1]] + t(save_Lnu_mat[[2]]) %*% save_Lnu_mat[[2]]
  # 
  # # Plot cluster:
  # rownames(Vhat_total) = cumsp_data$sp_name[1:n_sp]
  # colnames(Vhat_total) = cumsp_data$sp_name[1:n_sp]
  # Dist = dist(Vhat_total, diag=TRUE, upper=TRUE)
  # png(file = file.path(plot_folder, paste0('Cluster_Nu', cp, img_type)), width = img_width, 
  #     height = 120, res = img_res, units = 'mm')
  # plot(hclust(Dist))
  # dev.off()
  # 
  # xa = Matrix::chol(Vhat_total, pivot = TRUE)
  # xa = Matrix::chol(Vhat_total)
  # matrixcalc::is.positive.semi.definite(Vhat_total)
  # expand(Matrix::lu(Vhat_total))
  
  # -------------------------------------------------------------------------
  # Association sp to factors (space_term, theta):
  if(n_comps == 1) theta_slot_vec = 'theta_z'
  if(n_comps == 2) theta_slot_vec = c('theta_z', 'theta2_z')
  for(cp in seq_along(theta_slot_vec)) {
    theta_slot = theta_slot_vec[cp]
    Lhat_cf = matrix(0, nrow = n_sp, ncol = n_fac)
    theta2_vec = as.list(jtVModel$sdrep, what="Estimate")[[theta_slot]][1:(n_fac*n_sp - sum(0:(n_fac-1)))]
    Lhat_cf[lower.tri(Lhat_cf, diag=TRUE)] = theta2_vec
    Lhat_cf_rot = rotate_pca( L_tf = Lhat_cf )$L_tf
    # Var explained df:
    varex_vec = numeric(n_fac)
    for(k in 1:n_fac) varex_vec[k] = paste0(round(100*sum(Lhat_cf_rot[,k]^2)/sum(Lhat_cf_rot^2),1), '%')
    # Continue..
    dimnames(Lhat_cf_rot) = list( selsp_data$sp_name[1:n_sp],
                              paste0("Factor ", 1:ncol(Lhat_cf_rot), " (", varex_vec, ")") )
    load_df = as.data.frame(Lhat_cf_rot)
    load_df$species = rownames(load_df)
    load_df = pivot_longer(load_df, cols = starts_with("Factor"), names_to = 'Factor', values_to = 'loading')

    p1 = ggplot(data = load_df) + geom_col(aes(x = species, y = loading), position = 'dodge') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            strip.background = element_rect(fill="white")) +
      labs(y = 'Factor loading', x = NULL) +
      facet_wrap(~ Factor, scales = 'free_y') 
    ggsave(filename = paste0('Loading_Theta', cp, img_type), path = this_plot_folder, plot = p1, 
           width = img_width, height = 130, units = 'mm', dpi = img_res)
    
    # Plot Correlation matrix:
    Cov_omega2 = Lhat_cf %*% t(Lhat_cf)
    cor_mat = cov2cor(Cov_omega2)
    dimnames(cor_mat) = list( selsp_data$sp_name[1:n_sp], selsp_data$sp_name[1:n_sp])
    png(file = file.path(this_plot_folder, paste0('Corrplot_Theta', cp, img_type)), width = img_width, 
        height = 150, res = img_res, units = 'mm')
    corrplot::corrplot(cor_mat , method="circle", type="lower", diag = FALSE)
    dev.off()
    
    # Plot cluster:
    rownames(Lhat_cf) = selsp_data$sp_name[1:n_sp]
    Dist = dist(Lhat_cf, diag=TRUE, upper=TRUE)
    png(file = file.path(this_plot_folder, paste0('Cluster_Theta', cp, img_type)), width = img_width, 
        height = 120, res = img_res, units = 'mm')
    plot(hclust(Dist))
    dev.off()
  }
  
  cat("Nfac finished:", n_fac, "\n")
  
  } # conditional
  
} # n_fac loop  
  