rm(list = ls())

# Define type of school to be analyzed:
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
# Load data
wtData = readRDS(file.path(data_folder, 'weight_data.rds'))
effPoints = readRDS(file.path(data_folder, 'effPoints.rds'))

# -------------------------------------------------------------------------
# Method 1: raised by production (see Amande et al 2010): 10.1051/alr/2011003
# Strata are grid/year/quarter. 

estimateProd = calculate_ratio_bycatch(obs_df = wtData, eff_df = effPoints, type = 'production')
estimateProd = estimateProd %>% rename(est_prod = est)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Method 2: raised by effort (number of sets)
# Strata are grid/year/quarter. 

estimateSets = calculate_ratio_bycatch(obs_df = wtData, eff_df = effPoints, type = 'n_sets')
estimateSets = estimateSets %>% rename(est_sets = est)

# -------------------------------------------------------------------------
# Save ratio estimates:
outData = left_join(estimateProd, estimateSets, by = c('year', 'quarter', 'sp_name'))
saveRDS(outData, file = file.path(data_folder, 'ratio_estimates.rds'))

# -------------------------------------------------------------------------
# Compare plot by year:
plot_data = outData %>% 
  pivot_longer(cols = c(est_prod, est_sets), names_to = 'est_type', values_to = 'est_value') %>%
  mutate(est_type = factor(est_type, levels = c('est_prod', 'est_sets'), 
                           labels = c('Production', 'Effort (sets)')))
plot_data = plot_data %>% group_by(year, sp_name, est_type) %>%
  summarise(est_value = sum(est_value, na.rm = TRUE))

ggplot(plot_data, aes(x = year, y = est_value, color = est_type)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(from = 2014, to = 2022, by = 4)) +
  facet_wrap(~ sp_name, scales = 'free_y') +
  labs(x = 'Year', y = 'Estimated bycatch (tons)', color = 'Estimation type') +
  theme_classic() +
  theme(legend.position = c(0.8, 0.03))
ggsave(file.path(plot_folder, 'ratio_estimates_by_year.png'), width = img_width*1.5, 
       height = 180, units = 'mm', dpi = img_res)
