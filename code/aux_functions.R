
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
  
  CIrib = FALSE
  if((deparse(substitute(var_lwr)) %in% colnames(plot_data)) & (deparse(substitute(var_upr)) %in% colnames(plot_data))) CIrib = TRUE
  add_nominal = TRUE
  if(is.null(obs_data)) add_nominal = FALSE
  
  p1 = ggplot(data = plot_data, aes(x = {{var_x}}, y = {{var_y}})) +
    geom_line(aes(color = {{var_type}})) +
    ylab(yLab) + xlab(NULL) 
  if(CIrib) {
    p1 = p1 + geom_ribbon(aes(ymin = {{var_lwr}}, ymax = {{var_upr}}, fill = {{var_type}}), alpha = 0.3) + 
    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))
  }
  if(add_nominal) p1 = p1 + geom_point(data = obs_data, aes(x = {{var_x}}, y = {{var_y}}), color = 'red') 
  if(!add_legend) p1 = p1 + theme(legend.position = 'none')
  
  return(p1)
  
}

