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
# Calculate relative error:
est_sim = est_sim %>% mutate(re = ((est - true)/true)*100)
# Remove NaN because true = 0:
est_sim = est_sim %>% filter(!is.nan(re)) %>% mutate(re = ((est - true)/true)*100)

# -------------------------------------------------------------------------
# Bias by year and species:
bias_df = est_sim %>% group_by(sp_name, year, samp_frac) %>%
          summarise(bias = median(re), .groups = 'drop')
plot_dat = expand.grid(sp_name = sort(unique(est_sim$sp_name)),
                       year = sort(unique(est_sim$year)),
                       samp_frac = sort(unique(est_sim$samp_frac)))
plot_dat = left_join(plot_dat, bias_df, by = c("sp_name", "year", "samp_frac"))

# Plot:
mycolpal = rainbow(n = length(unique(bias_df$year)))

p1 = ggplot(data = plot_dat, aes(x = year, y = bias)) +
  geom_line(aes(color = factor(samp_frac))) +
  geom_point(aes(color = factor(samp_frac)), size = 0.8) +
  scale_color_manual(values = mycolpal) +
  scale_x_continuous(breaks = seq(from = min(bias_df$year), to = max(bias_df$year), by = 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        legend.position = c(0.7, 0.03), legend.direction = "horizontal") +
  labs(x = NULL, y = "Bias (%)", color = "Sampling fraction") +
  facet_wrap(~ sp_name) 
ggsave(paste0('bias', img_type), plot = p1, path = plot_folder,
       width = img_width*1.25, height = 160, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Precision by year and species:
precision_df = est_sim %>% group_by(sp_name, year, samp_frac) %>%
  summarise(precision = sd(re), .groups = 'drop')
plot_dat = expand.grid(sp_name = sort(unique(est_sim$sp_name)),
                       year = sort(unique(est_sim$year)),
                       samp_frac = sort(unique(est_sim$samp_frac)))
plot_dat = left_join(plot_dat, precision_df, by = c("sp_name", "year", "samp_frac"))

# Plot:
mycolpal = rainbow(n = length(unique(bias_df$year)))

p1 = ggplot(data = plot_dat, aes(x = year, y = precision)) +
  geom_line(aes(color = factor(samp_frac))) +
  geom_point(aes(color = factor(samp_frac)), size = 0.8) +
  scale_color_manual(values = mycolpal) +
  scale_x_continuous(breaks = seq(from = min(bias_df$year), to = max(bias_df$year), by = 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        legend.position = c(0.7, 0.03), legend.direction = "horizontal") +
  labs(x = NULL, y = "Precision", color = "Sampling fraction") +
  facet_wrap(~ sp_name, scales = 'free_y') 
ggsave(paste0('precision', img_type), plot = p1, path = plot_folder,
       width = img_width*1.25, height = 160, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Plot boxplot aggregated over years

# Aggregate over years:
agg_est_sim = est_sim %>% group_by(sp_name, sim, samp_frac) %>% 
  summarise(re = median(re, na.rm = TRUE), .groups = 'drop')
agg_est_sim = agg_est_sim %>% filter(!is.nan(re))
agg_est_sim = agg_est_sim %>% mutate(samp_frac = factor(samp_frac, levels = sort(unique(samp_frac)))) 

# Plot:
p1 = ggplot(data = agg_est_sim, aes(x = samp_frac, y = re)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  labs(x = "Sampling fraction", y = "Relative error (%)") +
  facet_wrap(~ sp_name, scales = 'free_y') 
ggsave(paste0('rel_error', img_type), plot = p1, path = plot_folder,
       width = img_width*1.5, height = 160, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Plot bias and precision by sampling fraction:
summ_est = agg_est_sim %>% group_by(sp_name, samp_frac) %>% 
  summarise(Bias = mean(re), 
            # Precision = quantile(re, probs = 0.95) - quantile(re, probs = 0.05),
            Precision = sd(re),
            .groups = 'drop')
summ_est = summ_est %>% mutate(samp_frac = as.numeric(as.character(samp_frac)))
summ_est = summ_est %>% pivot_longer(cols = c("Bias", "Precision"), 
                                 names_to = "metric", values_to = "value")

# Plot:
p1 = ggplot(data = summ_est, aes(x = samp_frac, y = value)) +
  geom_line(aes(group = sp_name), color = 'gray50') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  labs(x = "Sampling fraction", y = "Value") +
  scale_x_continuous(breaks = sort(unique(summ_est$samp_frac))) +
  facet_wrap(~ metric, scales = 'free_y') 
ggsave(paste0('summ_rel', img_type), plot = p1, path = plot_folder,
       width = img_width, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Plot number of years with missing information:

ny_true = est_sim %>% filter(samp_frac == 1) %>% 
  group_by(sp_name, sim) %>% 
  summarise(n_years_true = sum(true > 0), .groups = 'drop')
ny_est = est_sim %>% group_by(sp_name, samp_frac, sim) %>% 
  summarise(n_years_est = sum(est > 0), .groups = 'drop')
ny_merged = left_join(ny_est, ny_true, by = c("sp_name", "sim"))
ny_merged = ny_merged %>% mutate(ny_diff = n_years_true - n_years_est)
ny_merged = ny_merged %>% mutate(samp_frac = factor(samp_frac, levels = sort(unique(samp_frac))),
                                 ny_diff = factor(ny_diff)) 

# Plot:
mycolpal = rainbow(n = length(unique(ny_merged$ny_diff)))

p1 = ggplot(data = ny_merged, aes(x = samp_frac)) +
  geom_bar(aes(color = ny_diff, fill = ny_diff), width = 0.8) +
  scale_color_manual(values = mycolpal) +
  scale_fill_manual(values = mycolpal) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  labs(x = "Sampling fraction", y = "% simulations", fill = NULL, color = NULL) +
  facet_wrap(~ sp_name, scales = 'free_y') 
ggsave(paste0('nyr_missing', img_type), plot = p1, path = plot_folder,
       width = img_width*1.5, height = 160, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Plot number of species missing:

nsp_true = est_sim %>% filter(samp_frac == 1, true > 0) %>% 
  group_by(year, sim) %>% 
  summarise(n_sp_true = n_distinct(sp_name), .groups = 'drop')
nsp_est = est_sim %>% filter(est > 0) %>% group_by(year, samp_frac, sim) %>% 
  summarise(n_sp_est = n_distinct(sp_name), .groups = 'drop')
nsp_merged = left_join(nsp_true, nsp_est, by = c("year", "sim"))
nsp_merged = nsp_merged %>% mutate(nsp_diff = n_sp_true - n_sp_est)
nsp_merged = nsp_merged %>% mutate(samp_frac = factor(samp_frac, levels = sort(unique(samp_frac))),
                                   nsp_diff = factor(nsp_diff)) 

# Plot:
mycolpal = rainbow(n = length(unique(nsp_merged$nsp_diff)))

p1 = ggplot(data = nsp_merged, aes(x = samp_frac)) +
  geom_bar(aes(color = nsp_diff, fill = nsp_diff), width = 0.8) +
  scale_color_manual(values = mycolpal) +
  scale_fill_manual(values = mycolpal) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  labs(x = "Sampling fraction", y = "% simulations", fill = NULL, color = NULL) +
  guides(color = guide_legend(ncol = 2)) +
  facet_wrap(~ year, scales = 'free_y') 
ggsave(paste0('nsp_missing', img_type), plot = p1, path = plot_folder,
       width = img_width, height = 110, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Now read sdmTMB values:
# these_sp = sort(unique(est_sim$sp_name))
# true_values = list()
# for(isp in seq_along(these_sp)) {
#   # Read sdmTMB values:
#   true_values[[isp]] = read.csv(file = file.path("model/1_single", this_type, these_sp[isp], 'pred_est_time.csv'))
#   true_values[[isp]]$sp_name = these_sp[isp]
# }
# # Merge:
# true_values = bind_rows(true_values)
# true_values = true_values %>% group_by(sp_name) %>% summarise(est = sum(est), .groups = 'drop')

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

# Plot:
# p1 = ggplot(data = est_sim, aes(x = samp_frac, y = est)) +
#   geom_boxplot(outlier.size = 0.5) +
#   geom_hline(data = true_values, aes(yintercept = est), color = 'red', linetype = 'dashed') +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
#   labs(x = "Sampling fraction", y = "Estimated bycatch") +
#   facet_wrap(~ sp_name, scales = 'free_y') 
# ggsave(paste0('estimates', img_type), plot = p1, path = plot_folder,
#        width = img_width*1.5, height = 160, units = 'mm', dpi = img_res)
