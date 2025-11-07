rm(list = ls())

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')
nSims = 100 # number of sims run
colorPal = c("blue", "red") # for estimator

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
# Define factors:
est_sim = est_sim %>% mutate(samp_frac = factor(samp_frac, levels = frac_vector,
                                      labels = paste0(frac_vector*100, "%")),
                   sp_name = factor(sp_name, levels = sp_df$sp_levels) )

# -------------------------------------------------------------------------
# Report convergence rate of model-based estimator:
# Just select one year:
plot_data = est_sim %>% group_by(sim, samp_frac, sp_name) %>% summarise(category = mean(category, na.rm=TRUE))
plot_data = plot_data %>% mutate(category = factor(category, levels = 1:4, 
                                                   labels = c("Omega+Epsilon", "Omega", "Model failed", "Model not run")))

c1 = ggplot(data = plot_data, aes(x = samp_frac, fill = category)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d() + 
  scale_y_continuous(labels = scales::percent) +
  xlab("Sampling coverage") + ylab("Percentage of replicates") +
  guides(fill = guide_legend(title = NULL)) +
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9),
        axis.text.y = element_text(size = 10)) +
  facet_wrap(~ sp_name, ncol = 3)
ggsave(paste0('conv_model', img_type), plot = c1, path = plot_folder, 
       width = img_width , height = 200, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Calculate relative error:
re_data = est_sim %>% filter(!(category == 3)) %>%
            mutate(re_ratio = ((est_ratio - true)/true)*100,
                   re_model = ((est_model - true)/true)*100)
# when true = 0, NaN are produced. In this case, ratio and model based will be zero as well:
# fill them with 0
re_data = re_data %>% mutate(re_ratio = if_else(is.nan(re_ratio), 0, re_ratio),
                             re_model = if_else(is.nan(re_model), 0, re_model))

# New format for plotting:
re_data = re_data %>% pivot_longer(cols = c("re_ratio", "re_model"),
                                   values_to = "re", names_to = c("est_type"))

# Summarise global:
agg_data = re_data %>% group_by(samp_frac, sp_name, est_type, year) %>%
  dplyr::summarise(q025 = quantile(re, probs = 0.025), 
                   q50 = quantile(re, probs = 0.5),
                   q975 = quantile(re, probs = 0.975)) 
agg_data = agg_data %>% mutate(est_type = factor(est_type, levels = c("re_ratio", "re_model"),
                                                 labels = c("Ratio", "Model-based")))
# Plot aggregated
plot_dat = agg_data %>%
  group_by(samp_frac, sp_name, est_type) %>%
  dplyr::summarise(q025 = median(q025), q50 = median(q50), q975 = median(q975))


p2 = ggplot(plot_dat, aes(x=samp_frac, y=q50, colour=est_type)) +
  geom_linerange(aes(ymin = q025, ymax = q975), position=position_dodge(0.4)) +
  geom_pointrange(aes(ymin = q025, ymax = q975),  
                  position=position_dodge(0.4), size = 0.2) +
  scale_color_manual(values = colorPal) +
  geom_hline(yintercept=0, color=1, linetype='dashed') +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(colour = "Estimator") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.text=element_text(size=10)) +
  xlab("Sampling coverage") + ylab("Relative error (%)") +
  facet_wrap(~sp_name, ncol = 3)
ggsave(paste0('agg_re', img_type), plot = p2, path = plot_folder, 
       width = img_width , height = 200, units = 'mm', dpi = img_res)


# Plot by year
# Make this plot by species group:
all_sp_types = unique(sp_df$sp_type)
for(k in seq_along(all_sp_types)) {
  plot_dat = agg_data %>% filter(sp_name %in% (sp_df %>% filter(sp_type == all_sp_types[k]) %>% pull(sp_levels)))
  
  p3 = ggplot(plot_dat, aes(x=year, y=q50)) +
    geom_line(aes(color = est_type)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = est_type), alpha = 0.3) +
    geom_hline(yintercept=0, color=1, linetype='dashed') +
    scale_color_manual(values = colorPal) +
    scale_fill_manual(values = colorPal) +
    coord_cartesian(ylim = c(-100, 100)) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill="white"),
          legend.text=element_text(size=10)) +
    xlab(NULL) + ylab("Relative error (%)") +
    facet_grid(samp_frac~sp_name) +
    guides(colour=guide_legend(title="Estimator"), fill=guide_legend(title="Estimator"))
  ggsave(paste0('yr_re_', gsub(pattern = " ", replacement = "", x = all_sp_types[k]), img_type), 
         plot = p3, path = plot_folder, width = img_width , height = 200, units = 'mm', dpi = img_res)
}


# # -------------------------------------------------------------------------
# # Plot number of years with missing information:
# 
# ny_true = est_sim %>% filter(samp_frac == 1) %>% 
#   group_by(sp_name, sim) %>% 
#   summarise(n_years_true = sum(true > 0), .groups = 'drop')
# ny_est = est_sim %>% group_by(sp_name, samp_frac, sim) %>% 
#   summarise(n_years_est = sum(est > 0), .groups = 'drop')
# ny_merged = left_join(ny_est, ny_true, by = c("sp_name", "sim"))
# ny_merged = ny_merged %>% mutate(ny_diff = n_years_true - n_years_est)
# ny_merged = ny_merged %>% mutate(samp_frac = factor(samp_frac, levels = sort(unique(samp_frac))),
#                                  ny_diff = factor(ny_diff)) 
# 
# # Plot:
# mycolpal = rainbow(n = length(unique(ny_merged$ny_diff)))
# 
# p1 = ggplot(data = ny_merged, aes(x = samp_frac)) +
#   geom_bar(aes(color = ny_diff, fill = ny_diff), width = 0.8) +
#   scale_color_manual(values = mycolpal) +
#   scale_fill_manual(values = mycolpal) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
#   labs(x = "Sampling fraction", y = "% simulations", fill = NULL, color = NULL) +
#   facet_wrap(~ sp_name, scales = 'free_y') 
# ggsave(paste0('nyr_missing', img_type), plot = p1, path = plot_folder,
#        width = img_width*1.5, height = 160, units = 'mm', dpi = img_res)
# 
# 
# # -------------------------------------------------------------------------
# # Plot number of species missing:
# 
# nsp_true = est_sim %>% filter(samp_frac == 1, true > 0) %>% 
#   group_by(year, sim) %>% 
#   summarise(n_sp_true = n_distinct(sp_name), .groups = 'drop')
# nsp_est = est_sim %>% filter(est > 0) %>% group_by(year, samp_frac, sim) %>% 
#   summarise(n_sp_est = n_distinct(sp_name), .groups = 'drop')
# nsp_merged = left_join(nsp_true, nsp_est, by = c("year", "sim"))
# nsp_merged = nsp_merged %>% mutate(nsp_diff = n_sp_true - n_sp_est)
# nsp_merged = nsp_merged %>% mutate(samp_frac = factor(samp_frac, levels = sort(unique(samp_frac))),
#                                    nsp_diff = factor(nsp_diff)) 
# 
# # Plot:
# mycolpal = rainbow(n = length(unique(nsp_merged$nsp_diff)))
# 
# p1 = ggplot(data = nsp_merged, aes(x = samp_frac)) +
#   geom_bar(aes(color = nsp_diff, fill = nsp_diff), width = 0.8) +
#   scale_color_manual(values = mycolpal) +
#   scale_fill_manual(values = mycolpal) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
#   labs(x = "Sampling fraction", y = "% simulations", fill = NULL, color = NULL) +
#   guides(color = guide_legend(ncol = 2)) +
#   facet_wrap(~ year, scales = 'free_y') 
# ggsave(paste0('nsp_missing', img_type), plot = p1, path = plot_folder,
#        width = img_width, height = 110, units = 'mm', dpi = img_res)
