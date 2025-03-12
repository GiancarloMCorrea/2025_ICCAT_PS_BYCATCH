rm(list = ls())
# -------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Define folder to save results
plot_folder = file.path('figures')
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
set_types = c('FOB', 'FSC')
est_df = list()
obs_df = list()
for(i in seq_along(set_types)) {
  # Read estimates:
  all_files = list.files(file.path('figures', set_types[i], 'tinyVAST'))
  all_files = all_files[-grep(pattern = '.png', x = all_files)] # remove png files
  # Start extracting bycatch estimates by species:
  save_df = list()
  for(k in seq_along(all_files)) {
    tmp_file = file.path('figures', set_types[i], 'tinyVAST', all_files[k], 'Bycatch_est_ts.csv')
    cond_file = file.exists(tmp_file)
    if(cond_file) {
      tmp_df = read.csv(tmp_file)
      tmp_df$set = set_types[i]
      save_df[[k]] = tmp_df
    }
  }
  # Extract observations:
  obs_save_df = list()
  for(k in seq_along(all_files)) {
    tmp_df = read.csv(file.path('data', set_types[i], 'bycatch_est_obs', paste0(all_files[k], '.csv')))
    tmp_df$set = set_types[i]
    tmp_df$species = all_files[k]
    obs_save_df[[k]] = tmp_df
  }
  # Merge all estimates and obs:
  est_df[[i]] = bind_rows(save_df)
  obs_df[[i]] = bind_rows(obs_save_df)
}

# Names:
names(est_df) = set_types
names(obs_df) = set_types

# -------------------------------------------------------------------------
# Make plots FOB:
this_set = 'FSC'

# These datasets:
df_est = est_df[[this_set]]
df_obs = obs_df[[this_set]]

# Make plot:
p1 = ggplot(data = df_est, aes(x = Year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) + 
  geom_point(data = df_obs, aes(x = Year, y = est), color = 'red', size = 0.75) +
  ylab('Bycatch estimates (t)') + xlab(NULL) +
  scale_y_continuous(expand = c(0, 0.15)) +
  scale_x_continuous(breaks = c(2014, 2018, 2022)) +
  coord_cartesian(ylim = c(0, NA)) + 
  theme(strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7)) +
  facet_wrap(~species, scales = 'free_y')
ggsave(filename = paste0('estimates_all_taxa', img_type), path = file.path(plot_folder, this_set), 
       plot = p1, width = img_width, height = 140, units = 'mm', dpi = img_res)
