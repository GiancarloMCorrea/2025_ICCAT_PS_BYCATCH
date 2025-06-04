
# -------------------------------------------------------------------------
# Make grid plot:

add_sf_map = function(my_plot) {
  
  out_plot = my_plot + 
    geom_sf(data = worldmap, fill = "gray60", color = "gray60") +
    coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
    xlab(NULL) + ylab(NULL) +
    scale_y_continuous(breaks = yBreaks) + scale_x_continuous(breaks = xBreaks) 
  return(out_plot)
  
}


# -------------------------------------------------------------------------
calculate_area_on_land = function(dat) {
  
  library(rgeos)
  library(rnaturalearth)
  world_map = rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))
  wm = as(world_map, "Spatial")
  cs = gUnaryUnion(wm, id=as.character(wm$continent))
  cs_sf = st_as_sf(cs)
  inter_grid = st_intersection(cs_sf, dat)
  if(nrow(inter_grid) > 0) area_on_land = sum(as.numeric(st_area(inter_grid)))*1e-06 # in km2
  else area_on_land = 0
  
  return(area_on_land)
  
}


# -------------------------------------------------------------------------

plot_predictions = function(plot_data, nCol = 4) {
  
  # Plot CPUE by grid and year:
  plot_data = plot_data %>% group_by(year, ID) %>% summarise(bycatch_est = sum(bycatch_est), .groups = 'drop') 
  plot_data = left_join(MyGrid, plot_data, by = 'ID')
  plot_data = plot_data %>% na.omit
  plot_data = plot_data %>% mutate(est_disc = cut(bycatch_est, 
                                                  breaks = c(0, quantile(x = plot_data$bycatch_est, probs = c(0.5, 0.75, 0.95)), Inf), 
                                                  right=FALSE))
  
  p1 = ggplot() +  
    geom_sf(data = plot_data, aes(fill = est_disc, color = est_disc)) + 
    scale_fill_viridis_d() + scale_color_viridis_d() 
  p1 = add_sf_map(p1)
  p1 = p1 + labs(fill = "Estimates") + theme(legend.position = 'bottom') + 
    guides(color = 'none') +
    facet_wrap(~ year, ncol = nCol)
  
  return(p1)
  
}

# -------------------------------------------------------------------------
plot_time_predictions = function(model_est, ratio_est) {
  
  colvals = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
  # Legend for plot:
  annotations <- data.frame(
    xpos = c(Inf, Inf, Inf),
    ypos =  c(Inf, Inf, Inf),
    annotateText = c("Production", "Number of sets", "sdmTMB"),
    hjustvar = c(1, 1, 1) ,
    vjustvar = c(1.5, 3.5, 5.5))
  
  p1 = ggplot(data = ratio_est) +
    geom_point(aes(x = year, y = est_prod), color = colvals[1]) +
    geom_line(aes(x = year, y = est_prod), color = colvals[1]) +
    geom_point(aes(x = year, y = est_sets), color = colvals[2]) +
    geom_line(aes(x = year, y = est_sets), color = colvals[2]) +
    geom_pointrange(data = model_est, aes(x = year, y = est,
                                          ymin = lwr, ymax = upr),
                    color = colvals[3], size = 0.25) +
    scale_x_continuous(breaks = seq(from = 2014, to = 2022, by = 4)) +
    labs(x = 'Year', y = 'Estimated bycatch (tons)') +
    geom_text(data=annotations[1,],aes(x=xpos,y=ypos,hjust=hjustvar, vjust=vjustvar,label=annotateText), color = colvals[1]) +
    geom_text(data=annotations[2,],aes(x=xpos,y=ypos,hjust=hjustvar, vjust=vjustvar,label=annotateText), color = colvals[2]) +
    geom_text(data=annotations[3,],aes(x=xpos,y=ypos,hjust=hjustvar, vjust=vjustvar,label=annotateText), color = colvals[3]) +
    theme(legend.position = 'none') + ggtitle(this_sp) 
  
  return(p1)
  
}


# -------------------------------------------------------------------------
# Plot darhma residuals:
plot_residuals = function(this_model){
  
  sim_res = simulate(this_model, nsim = 500, type = "mle-mvn")
  check_res = sdmTMB::dharma_residuals(sim_res, this_model, return_DHARMa = TRUE)
  png(filename = file.path(this_plot_folder, paste0('res_dharma', img_type)), 
      width = img_width, height = 100, units = 'mm', res = img_res)
  par(mar = c(4, 4, 1, 0.5))
  plot(check_res, title = NULL)
  dev.off()
  
}



# -------------------------------------------------------------------------
# Plot omega:

plot_omega = function(this_model){
  
  tmp_df = list()
  for(j in 1:n_comps) {
    tmp_df[[j]] = data.frame(Lon = this_model$spde$mesh$loc[,1], 
                             Lat = this_model$spde$mesh$loc[,2], 
                             omega = this_model$parlist$omega_s[,j],
                             comp = paste0('component_', j))
  }
  plot_dat = bind_rows(tmp_df)
  plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  p1 = ggplot(plot_dat) + geom_sf(aes(color = omega), size = 1.5) + scale_colour_gradient2() + labs(color = 'Omega')
  p1 = add_sf_map(p1)
  p1 = p1 + facet_wrap(~comp)
  return(p1)
  
}


# -------------------------------------------------------------------------
# Plot epsilon:

plot_epsilon = function(this_model){
  
  tmp_df = list()
  save_plots = list()
  all_years = this_model$time_lu$time_from_data
  for(k in 1:n_comps) {
    for(i in seq_along(all_years)) {
      tmp_df[[i]] = data.frame(year = all_years[i],
                               Lon = this_model$spde$mesh$loc[,1], 
                               Lat = this_model$spde$mesh$loc[,2], 
                               epsilon = this_model$parlist$epsilon_st[,i,k])
    }
    plot_dat = bind_rows(tmp_df)
    plot_dat = plot_dat %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
    p1 = ggplot(plot_dat) + geom_sf(aes(color = epsilon), size = 1) + scale_colour_gradient2() + labs(color = 'Epsilon')
    p1 = add_sf_map(p1)
    p1 = p1 + facet_wrap(~ year)
    save_plots[[k]] = p1
  }

  return(save_plots)
  
}


# -------------------------------------------------------------------------
# Function to make the time term tinyVAST:
make_time_term = function(n_sp, n_fac, par_lab = 'a', save_folder = NULL) {
  
  tm_p1 = NULL
  for(k in 1:n_fac){
    tm_p1 = c(tm_p1, 
              paste(
                rep(paste0('f', k), times = n_sp-k+1),
                (1+k-1):n_sp,
                sep = ' -> '))
  }
  tm_p1 = paste(tm_p1, paste0(par_lab, 1:(n_sp*n_fac - sum(0:(n_fac-1)))), sep = ', 0, ')
  tm_p2 = paste0(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' -> '),
    ', 1, NA, 1'
  )
  tm_p3 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' <-> '),
    ', 0, NA, 0'
  )
  tm_p4 = paste(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' <-> '),
    paste0('sd_', par_lab, '_f', 1:n_fac),
    sep = ', 0, '
  )
  tm_file = c(tm_p1, tm_p2, tm_p3, tm_p4)
  if(!is.null(save_folder)) writeLines(tm_file, file.path(save_folder, paste0("time_input_sp", n_sp, ".txt")))
  tm_mod = paste0("\n  ", paste(tm_file, collapse = '\n  '), "\n")
  return(tm_mod)
  
}

# -------------------------------------------------------------------------
# Function to make the space term tinyVAST:
make_space_term = function(n_sp, n_fac, par_lab = 'o', save_folder = NULL) {
  
  sem_p1 = NULL
  for(k in 1:n_fac){
    sem_p1 = c(sem_p1, 
               paste(
                 rep(paste0('f', k), times = n_sp-k+1),
                 (1+k-1):n_sp,
                 sep = ' -> '))
  }
  sem_p1 = paste(sem_p1, paste0(par_lab, 1:(n_sp*n_fac - sum(0:(n_fac-1)))), sep = ', ')
  sem_p2 = paste(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' <-> '),
    paste0('sd_', par_lab, '_f', 1:n_fac),
    sep = ', '
  )
  sem_p3 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' <-> '),
    ', NA, 0'
  )
  sem_file = c(sem_p1, sem_p2, sem_p3)
  if(!is.null(save_folder)) writeLines(sem_file, file.path(save_folder, paste0("sem_input_sp", n_sp, ".txt")))
  sem_mod = paste0("\n  ", paste(sem_file, collapse = '\n  '), "\n")
  return(sem_mod)
  
}

# -------------------------------------------------------------------------
# Function to make the spacetime term tinyVAST:
make_spacetime_term = function(n_sp, n_fac, par_lab = 'e', save_folder = NULL) {
  
  dsem_p1 = NULL
  for(k in 1:n_fac){
    dsem_p1 = c(dsem_p1, 
                paste(
                  rep(paste0('f', k), times = n_sp-k+1),
                  (1+k-1):n_sp,
                  sep = ' -> '))
  }
  dsem_p1 = paste(dsem_p1, paste0(par_lab, 1:(n_sp*n_fac - sum(0:(n_fac-1)))), sep = ', 0, ')
  dsem_p2 = paste0(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' -> '),
    ', 1, NA, 0'
  )
  dsem_p3 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' -> '),
    ', 1, NA, 0'
  )
  dsem_p4 = paste(
    paste(paste0('f', 1:n_fac), paste0('f', 1:n_fac), sep = ' <-> '),
    paste0('sd_', par_lab, '_f', 1:n_fac),
    sep = ', 0, '
  )
  dsem_p5 = paste0(
    paste(1:n_sp, 1:n_sp, sep = ' <-> '),
    ', 0, NA, 0'
  )
  dsem_file = c(dsem_p1, dsem_p2, dsem_p3, dsem_p4, dsem_p5)
  if(!is.null(save_folder)) writeLines(dsem_file, file.path(save_folder, paste0("dsem_input_sp", n_sp, ".txt")))
  dsem_mod = paste0("\n  ", paste(dsem_file, collapse = '\n  '), "\n")
  return(dsem_mod)
  
}


# -------------------------------------------------------------------------
# Calculate Moran Index p-value:
get_moran = function(var_vec, lon, lat) {
  require(ape)
  require(fields)
  if(all(var_vec == 0)) {
    out = NA
  } else {
    coords = cbind(lon, lat)
    w = fields:::rdist(coords)
    moran_res = ape::Moran.I(x = var_vec, w = w)
    out = round(moran_res$p.value, digits = 3)
  }
  return(out)
}


# -------------------------------------------------------------------------
# Shorten species name:
short_sp_name = function(sp_name) {
  
  list_split = strsplit(sp_name," ")
  save_out = character(length(list_split))
  for(k in seq_along(list_split)) {
    n_words = length(list_split[[k]])
    if(n_words == 1) { 
      save_out[k] = list_split[[k]][1]
    } else if(n_words == 2) {
      save_out[k] = paste0(substring(text = list_split[[k]][1], first = 1, last = 1),
                        '. ', list_split[[k]][2])
    } else {
      save_out[k] = paste(list_split[[k]], collapse = " ")
    }
  }
  return(save_out)
  
}

# -------------------------------------------------------------------------
# Produce table with fixed effects sdmTMB:
get_summary_sdmTMB = function(model, n_comp, model_label = 'model') {
  save_df = list()
  for(i in 1:n_comp) {
    summdf = tidy(model, model = i) %>% mutate(z_score = estimate/std.error,
                                               p_value = 2*pnorm(-abs(z_score)))
    summdf = summdf %>% select(term, estimate, std.error, p_value) %>% 
      mutate(p_value = ifelse(p_value < 0.01, '<0.01', round(p_value, 2)),
             component = i, model = model_label)
    save_df[[i]] = summdf
  }
  out_df = bind_rows(save_df)
  return(out_df)
}

# -------------------------------------------------------------------------
# Remove terms from sdmTMB formula based on significance:
remove_terms_sdmTMB = function(model, n_comp, formula) {
  out_formula = formula
  for(k in 1:n_comp) {
    tmp_formula = as.character(formula[[k]])
    tmp_formula = paste(tmp_formula[2], tmp_formula[1], tmp_formula[3])
    tmp_formula = gsub(pattern = ' ', replacement = '', x = tmp_formula)
    summdf = tidy(model, model = k) %>% mutate(z_score = estimate/std.error,
                                               p_value = 2*pnorm(-abs(z_score)))
    nonsig_terms = summdf$term[which(summdf$p_value >= 0.05)]
    terms_vec = strsplit(x = tmp_formula, split = '[+]')[[1]]
    n_levels_qtr = length(levels(model$data$quarter))
    if(length(grep(pattern = 'quarter', x = nonsig_terms)) == (n_levels_qtr-1)) terms_vec = terms_vec[-which(terms_vec == 'quarter')]
    if(length(grep(pattern = 'sst', x = nonsig_terms)) == 1) terms_vec = terms_vec[-which(terms_vec == 'sst')]
    if(length(grep(pattern = 'trop_catch', x = nonsig_terms)) == 1) terms_vec = terms_vec[-which(terms_vec == 'trop_catch')]
    out_formula[[k]] = as.formula(paste(terms_vec, collapse = '+'))
  }
  return(out_formula)
}


# -------------------------------------------------------------------------
# Calculate ratio estimator:
calculate_ratio_bycatch = function(obs_df, eff_df, type = 'production') {
  
  # Define denominator here:
  if(type == 'production') {
    obs_df$denom = obs_df$trop_catch_nonstd
    eff_df$denom = eff_df$trop_catch_nonstd
  }
  if(type == 'n_sets') {
    obs_df$denom = 1
    eff_df$denom = 1
  }
  
  grData = obs_df %>% group_by(ID, year, quarter, sp_name) %>% 
    summarise(bycatch = sum(bycatch),
              denom = sum(denom), .groups = 'drop')
  raiseData = grData %>% mutate(ratio = bycatch/denom)
  # Replace Inf and NaN with zeros due to denom == 0:
  raiseData = raiseData %>% mutate(ratio = ifelse(is.infinite(ratio), 0, ratio))
  raiseData = raiseData %>% mutate(ratio = ifelse(is.nan(ratio), 0, ratio))
  
  # Calculate ratio per year/quarter/species for whole area:
  raiseData_tot = grData %>% group_by(year, quarter, sp_name) %>% 
    summarise(ratio_tot = sum(bycatch)/sum(denom), .groups = 'drop')
  
  # Summarise effort data:
  logbookData = eff_df %>% group_by(ID, year, quarter) %>% 
    summarise(denom = sum(denom), .groups = 'drop')
  
  # Merge dfs:
  raiseData = raiseData %>% select(ID, year, quarter, sp_name, ratio)
  mergedData = left_join(logbookData, raiseData, by = c('ID', 'year', 'quarter'))
  # I think this step is no needed...
  mergedData = left_join(mergedData, raiseData_tot, by = c('year', 'quarter', 'sp_name'))
  # Check if NA comes from grid with no bycatch:
  # identical( which(is.na(mergedData$sp_name)), which(is.na(mergedData$ratio)) )
  # Remove those NA safely:
  mergedData = mergedData %>% filter(!is.na(sp_name))
  # Check again if any NA and replace with tot ratio:
  if(any(is.na(mergedData$ratio))) {
    mergedData = mergedData %>% mutate(ratio = ifelse(is.na(ratio), ratio_tot, ratio))
  }
  
  # Estimates by year and quarter:
  mergedData = mergedData %>% mutate(est = ratio*denom)
  estimateData = mergedData %>% group_by(year, quarter, sp_name) %>% 
    summarise(est = sum(est, na.rm = TRUE), .groups = 'drop') 
  return(estimateData)
  
}



