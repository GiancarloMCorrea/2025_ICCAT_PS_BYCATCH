rm(list = ls())

# Define type of school to be analyzed:
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

# Make plot:
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
# Make plot only for billfishes:

ratio1_df = ratio_df %>% filter(sp_name %in% c('Istiophoridae', 'I. albicans', 'X. gladius', 'M. nigricans'))
model1_df = model_df %>% filter(sp_name %in% c('Istiophoridae', 'I. albicans', 'X. gladius', 'M. nigricans'))

# Make plot:
p1 = ggplot(ratio1_df, aes(x = year, y = est)) +
  geom_line() +
  geom_point() +
  geom_pointrange(data = model1_df, aes(x = year, y = est,
                                       ymin = lwr, ymax = upr),
                  color = 'red', size = 0.25) +
  scale_x_continuous(breaks = seq(from = 2014, to = 2022, by = 4)) +
  facet_wrap(~ sp_name, scales = 'free_y') +
  labs(x = 'Year', y = 'Estimated bycatch (tons)')
ggsave(file.path(plot_folder, 'compare_estimates_billfishes.png'), width = img_width, plot = p1,
       height = 130, units = 'mm', dpi = img_res)


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
coef_df = coef_df %>% mutate(family = if_else(category == 3, 'Tweedie', 'Hurdle approach'))

# Save data for plotting later
sp_data = coef_df %>% group_by(model) %>% summarise(family = unique(family))
std_data = rbind(expand.grid(model = unique(coef_df$model), family = 'Hurdle approach',
                             component = c('Component 1', 'Component 2')),
                 expand.grid(model = unique(coef_df$model), family = 'Tweedie',
                             component = c('Component 1')))

# Continue processing:
coef_df = coef_df %>% filter(!grepl(pattern = "fyear", x = coef_df$term))
coef_df = coef_df %>% filter(!grepl(pattern = "quarter", x = coef_df$term))
coef_df = coef_df %>% mutate(term = factor(term, levels = c('sst', 'trop_catch'), 
                                           labels = c('SST', 'Target catch')),
                             component = factor(component, levels = c('1', '2'), 
                                           labels = c('Component 1', 'Component 2')))
# Create base data frame with all sp
coef_df = coef_df %>% select(term, estimate, component, model, family)
std_data = replicate_df(std_data, time_name = 'term', time_values = c('Target catch', 'SST'))
plot_data = left_join(std_data, coef_df)
# Put zero for species with coeff=0
this_family = 'Hurdle approach'
these_sp = sp_data %>% filter(family == this_family) %>% pull(model)
plot_data = plot_data %>% mutate(estimate = if_else(family == this_family & model %in% these_sp & is.na(estimate), 0, estimate))
this_family = 'Tweedie'
these_sp = sp_data %>% filter(family == this_family) %>% pull(model)
plot_data = plot_data %>% mutate(estimate = if_else(family == this_family & model %in% these_sp & is.na(estimate), 0, estimate))

# Replace Xiphioidea by Istiophoridae (delete this later)
plot_data = plot_data %>% mutate(model = if_else(model == 'Xiphioidea', 'Istiophoridae', model))

p2 = ggplot(data = plot_data, aes(x = term, y = model, fill = estimate)) +
  geom_tile(color = NA, na.rm=TRUE) +
  xlab(NULL ) + ylab(NULL) +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(fill = expression(beta)) +
  theme_classic() +
  facet_nested(~ family + component)
ggsave(file.path(plot_folder, 'table_effect.png'), width = img_width, 
       plot = p2, height = 160, units = 'mm', dpi = img_res)
