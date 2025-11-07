rm(list = ls())

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')

# -------------------------------------------------------------------------

# Read observed sim data (only for one simulation):
i_count = 1
save_sim = list()
these_frac = list.files(file.path(model_folder, "obs_sim_data"))
for(k in seq_along(these_frac)) {
  these_sims = list.files(file.path(model_folder, "obs_sim_data", these_frac[k]))
  # Read all files:
  for(j in seq_along(these_sims)) { 
    save_sim[[i_count]] = readRDS(file = file.path(model_folder, "obs_sim_data", these_frac[k], these_sims[j]))
    i_count = i_count + 1
  }
}
# Merge data:
obs_sim = bind_rows(save_sim)
# Define factors for sample coverage:
obs_sim = obs_sim %>% mutate(samp_frac = factor(samp_frac, levels = frac_vector,
                                                labels = paste0(frac_vector*100, "%")))


# Plot:
dir.create(file.path(plot_folder, 'sim_spatial'), showWarnings = FALSE, recursive = TRUE)
these_sp = sort(unique(obs_sim$sp_name))
for(isp in seq_along(these_sp)) {
  plot_data = obs_sim %>% filter(sp_name == these_sp[isp])
  plot_data = plot_data %>% filter(bycatch > 0)
  plot_data = plot_data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  p2 = ggplot(plot_data) + geom_sf(aes(size = bycatch, color = bycatch), alpha = 0.5) +
            scale_color_viridis_c(guide = "legend") + 
            scale_size(range = c(0, 5)) +
            guides(colour = guide_legend(title = 'Bycatch (t)'), size = guide_legend(title = 'Bycatch (t)')) +
            theme(legend.position = c(0.85, 0.15),
                  strip.text = element_text(size = 12),
                  strip.background = element_blank()) + labs(title = these_sp[isp])
  p2 = add_sf_map(p2)
  p2 = p2 + facet_wrap(~ samp_frac, drop = FALSE)
  ggsave(paste0('map_', isp, img_type), plot = p2, path = file.path(plot_folder, 'sim_spatial'),
         width = img_width, height = 150, units = 'mm', dpi = img_res)
}
