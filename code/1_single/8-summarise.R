rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FSC' # FOB or FSC
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
# Read data:
ratio_df = readRDS(file = file.path(data_folder, 'ratio_estimates.rds'))
ratio_df = ratio_df %>% group_by(year, sp_name) %>%
  summarise(est = sum(est_prod, na.rm = TRUE)) # select production method
# Read model estimates:
all_sp = list.files(file.path(model_folder))
model_df = list()
for(k in seq_along(all_sp)) {
  this_file = file.path(model_folder, all_sp[k], 'pred_est_time.csv')
  if(file.exists(this_file)) model_df[[k]] = read.csv(this_file)
}
# Merge:
model_df = bind_rows(model_df)

# -------------------------------------------------------------------------
# Make plot comparing ratio and model estimates for all species:
p1 = ggplot(ratio_df, aes(x = year, y = est)) +
  geom_line() +
  geom_point() +
  geom_pointrange(data = model_df, aes(x = year, y = est,
                                       ymin = lwr, ymax = upr),
                  color = 'red', size = 0.25) +
  scale_x_continuous(breaks = seq(from = 2014, to = 2022, by = 4)) +
  facet_wrap(~ sp_name, scales = 'free_y') +
  labs(x = 'Year', y = 'Estimated bycatch (tons)')
ggsave(file.path(plot_folder, 'compare_estimates.png'), width = img_width*1.5, plot = p1,
       height = 180, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------

# Read model estimates:
all_sp = list.files(file.path(model_folder))
coef_df = list()
for(k in seq_along(all_sp)) {
  this_file = file.path(model_folder, all_sp[k], 'mod_summ.rds')
  if(file.exists(this_file)) coef_df[[k]] = readRDS(this_file)
}
# Merge:
coef_df = bind_rows(coef_df)
coef_df = coef_df %>% filter(!grepl(pattern = "fyear", x = coef_df$term))
coef_df = coef_df %>% filter(!grepl(pattern = "quarter", x = coef_df$term))
coef_df = coef_df %>% mutate(term = factor(term, levels = c('sst', 'trop_catch'), 
                                           labels = c('SST', 'Target catch')),
                             component = factor(component, levels = c('1', '2'), 
                                           labels = c('Component 1', 'Component 2')))
# Add asterisk after sp name for tweedie models:
coef_df = coef_df %>% mutate(model = if_else(category == 3, paste0(model, '*'), model))


p2 = ggplot(data = coef_df, aes(x = term, y = model, fill = estimate)) +
  geom_tile(color = NA) +
  xlab(NULL ) + ylab(NULL) +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(fill = expression(beta)) +
  theme_classic() +
  facet_grid(~ component)
ggsave(file.path(plot_folder, 'table_effect.png'), width = img_width, 
       plot = p2, height = 180, units = 'mm', dpi = img_res)
