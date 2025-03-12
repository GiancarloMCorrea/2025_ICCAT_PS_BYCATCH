
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

plot_predictions = function(plot_data, legend_position = c(0.65, 0.08), legend_title = 'Predicted CPUE', nCol = 5) {
  
  p1 = ggplot() +  
    geom_sf(data = plot_data, aes(fill = cpue_pred, color = cpue_pred)) + 
    scale_fill_viridis() + scale_color_viridis() +
    geom_sf(data = worldmap, fill = "gray60", color=NA) +
    coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
    xlab(NULL) + ylab(NULL) +
    theme(legend.position = legend_position, legend.direction="horizontal") +
    scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
    facet_wrap(~ year, ncol = nCol) +
    labs(fill = legend_title) + guides(color = 'none')
  
  return(p1)
  
}

plot_time_predictions = function(plot_data, obs_data = NULL, var_x, var_y, 
                                 var_lwr = NULL, var_upr = NULL,
                                 var_type,
                                 yLab = 'CPUE', add_legend = TRUE) {
  
  require(scales)
  CIrib = FALSE
  if((deparse(substitute(var_lwr)) %in% colnames(plot_data)) & (deparse(substitute(var_upr)) %in% colnames(plot_data))) CIrib = TRUE
  add_nominal = TRUE
  if(is.null(obs_data)) add_nominal = FALSE
  
  p1 = ggplot(data = plot_data, aes(x = {{var_x}}, y = {{var_y}})) +
    geom_line(aes(color = {{var_type}})) +
    ylab(yLab) + xlab(NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
	scale_y_continuous(expand = c(0, 0.15)) +
    coord_cartesian(ylim = c(0, NA))
  if(CIrib) {
    p1 = p1 + geom_ribbon(aes(ymin = {{var_lwr}}, ymax = {{var_upr}}, fill = {{var_type}}), alpha = 0.3) + 
    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
  }
  if(add_nominal) p1 = p1 + geom_point(data = obs_data, aes(x = {{var_x}}, y = {{var_y}}), color = 'red') 
  if(!add_legend) p1 = p1 + theme(legend.position = 'none')
  
  return(p1)
  
}


# -------------------------------------------------------------------------
# Function to make the time term:
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
# Function to make the space term:
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
# Function to make the spacetime term:
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
