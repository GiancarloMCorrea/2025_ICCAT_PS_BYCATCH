rm(list = ls())

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')

# -------------------------------------------------------------------------
# Read observed sim data (only for one simulation):
i_count = 1
save_sim = list()
these_frac = list.files(file.path(model_folder, "sim_est"))
for(k in seq_along(these_frac)) {
  these_sims = list.files(file.path(model_folder, "sim_est", these_frac[k]))
  # Read all files:
  for(j in seq_along(these_sims)) { 
    save_sim[[i_count]] = readRDS(file = file.path(model_folder, "sim_est", these_frac[k], these_sims[j]))
    i_count = i_count + 1
  }
}
# Merge data:
est_sim = bind_rows(save_sim)
# Aggregate over years:
est_sim = est_sim %>% group_by(sp_name, sim, samp_frac) %>% summarise(est = sum(est), .groups = 'drop')
est_sim = est_sim %>% mutate(samp_frac = factor(samp_frac, levels = sort(unique(samp_frac)))) 


# -------------------------------------------------------------------------
# Now read true values:
these_sp = sort(unique(est_sim$sp_name))
true_values = list()
for(isp in seq_along(these_sp)) {
  # Read true values:
  true_values[[isp]] = read.csv(file = file.path("model/1_single", this_type, these_sp[isp], 'pred_est_time.csv'))
  true_values[[isp]]$sp_name = these_sp[isp]
}
# Merge:
true_values = bind_rows(true_values)
true_values = true_values %>% group_by(sp_name) %>% summarise(est = sum(est), .groups = 'drop')

# -------------------------------------------------------------------------

# Plot:
p1 = ggplot(data = est_sim, aes(x = samp_frac, y = est)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = est), color = 'red', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  labs(x = "Sampling fraction", y = "Estimated bycatch") +
  facet_wrap(~ sp_name, scales = 'free_y') 
ggsave(paste0('estimates', img_type), plot = p1, path = plot_folder,
       width = img_width*1.5, height = 160, units = 'mm', dpi = img_res)
