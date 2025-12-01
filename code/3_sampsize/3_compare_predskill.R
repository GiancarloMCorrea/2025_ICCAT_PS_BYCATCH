rm(list = ls())

# Define type of school to be analyzed:
source('code/3_sampsize/load_libs.R')
nSims = 100 # number of sims run
colorPal = wes_palette('AsteroidCity1', n = 5, type = c("discrete"))[c(3,4)] # for estimator

# -------------------------------------------------------------------------
# Read observed sim data (only for one simulation):
i_count = 1
save_sim = list()
these_frac = list.files(file.path(model_folder, "crossval"))
for(k in seq_along(these_frac)) {
  these_sims = list.files(file.path(model_folder, "crossval", these_frac[k]))
  # Read all files:
  for(j in seq_along(these_sims)) { 
    save_sim[[i_count]] = readRDS(file = file.path(model_folder, "crossval", these_frac[k], these_sims[j]))
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
# Plot mae:
plot_dat = est_sim %>% filter(!(category == 3)) %>% mutate(rmse = round(rmse, digits = 3),
                                                           mae = round(mae, digits = 3))
if(this_type == 'FOB') maxY = 1
if(this_type == 'FSC') maxY = 3

p1 = ggplot(data = plot_dat, aes(x = samp_frac, y = rmse)) +
  geom_boxplot(outlier.size = 1) +
  xlab("Sampling coverage") + ylab("RMSE") +
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9),
        axis.text.y = element_text(size = 10)) +
  coord_cartesian(ylim = c(0, maxY)) +
  facet_wrap(~ sp_name, ncol = 3, scales = 'free_y')
ggsave(paste0('crossval', img_type), plot = p1, path = plot_folder, 
       width = img_width , height = 200, units = 'mm', dpi = img_res)

