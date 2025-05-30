rm(list = ls())

# Define type of school to be analyzed:
this_type = 'FOB' # or FSC
source('code/1_single/load_libs.R')

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path(data_folder, 'weight_data.rds'))
numbers_data = readRDS(file = file.path(data_folder, 'numbers_data.rds'))
effPoints = readRDS(file.path(data_folder, 'effPoints.rds'))
obsPoints = readRDS(file.path(data_folder, 'obsPoints.rds'))

# -------------------------------------------------------------------------
# Make figures dominance sp groups:
n_sp = 25 # first N species to plot

# Weight:
plot_data = weight_data %>% group_by(year, sp_name) %>% summarise(bycatch = sum(bycatch))
cumsp_data = plot_data %>% group_by(sp_name) %>% summarise(bycatch = sum(bycatch))
cumsp_data = arrange(cumsp_data, desc(bycatch))
cumsp_data = cumsp_data %>% dplyr::filter(bycatch > 0)
write.csv(cumsp_data, file = file.path(data_folder, 'cumsp_data_weight.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>% mutate(sp_name = factor(sp_name, levels = cumsp_data$sp_name))

p1 = ggplot(data = plot_data, aes(x = sp_name, y = bycatch)) +
  geom_col() + 
  xlab(NULL) + ylab('Weight (kg)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ year, scales = 'free_y')
ggsave(paste0('sp_dom_weight', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# Numbers:
plot_data = numbers_data %>% group_by(year, sp_name) %>% summarise(bycatch = sum(bycatch))
cumsp_data = plot_data %>% group_by(sp_name) %>% summarise(bycatch = sum(bycatch))
cumsp_data = arrange(cumsp_data, desc(bycatch))
cumsp_data = cumsp_data %>% dplyr::filter(bycatch > 0)
write.csv(cumsp_data, file = file.path(data_folder, 'cumsp_data_numbers.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>% mutate(sp_name = factor(sp_name, levels = cumsp_data$sp_name))

p1 = ggplot(data = plot_data, aes(x = sp_name, y = log(bycatch+1))) +
  geom_col() + 
  xlab(NULL) + ylab('log(Numbers)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ year)
ggsave(paste0('sp_dom_numbers', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Make plot year and quarter:
p1 = ggplot(obsPoints, aes(x = factor(year), fill = quarter)) + 
  geom_bar() + 
  scale_fill_viridis_d() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = 'bottom') +
  xlab(NULL) + ylab('Number of sets (observations)') +
  guides(fill = guide_legend(title = 'quarter'))
p2 = ggplot(effPoints, aes(x = factor(year), fill = quarter)) + 
  geom_bar() + 
  scale_fill_viridis_d() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = 'bottom') +
  xlab(NULL) + ylab('Number of sets (effort)') +
  guides(fill = guide_legend(title = 'quarter'))
my_plot = grid.arrange(p1, p2, ncol = 2)
ggsave(filename = 'n_obs.png', path = plot_folder, plot = my_plot, 
       width = img_width*1.5, height = 100, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Map points:
obsSF = obsPoints %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
effSF = effPoints %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(obsSF) + geom_sf(size = 1, alpha = 0.5) + ggtitle('Observations')
p1 = add_sf_map(p1)
p2 = ggplot(effSF) + geom_sf(size = 1, alpha = 0.5) + ggtitle('Effort')
p2 = add_sf_map(p2)
my_plot = grid.arrange(p1, p2, ncol = 2)
ggsave(paste0('map_sets', img_type), plot = my_plot, path = plot_folder,
       width = img_width*1.2, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Fishing locations OBS by season/month:
plot_data = obsPoints %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(plot_data) + geom_sf(aes(color = quarter), size = 0.5, alpha = 0.35)
p1 = add_sf_map(p1)
p1 = p1 + scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom') + facet_wrap(~ year)
ggsave(paste0('obs_map_sets_qtr', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 150, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Number bycatch sp per set
plot_data = weight_data %>% mutate(pos_by = if_else(bycatch > 0, 1, 0)) %>% 
  group_by(year, quarter, id_set) %>% 
  summarise(bycatch = sum(bycatch), n_sp_catch = sum(pos_by))

p1 = ggplot(plot_data, aes(x = factor(year), y = n_sp_catch)) + 
  geom_boxplot(aes(fill = quarter), width = 0.8, outlier.size = 0.7) +
  scale_fill_viridis_d() +
  xlab(NULL) + ylab('N-sp/set') +
  theme(legend.position = 'bottom')
ggsave(paste0('obs_n-sp_set', img_type), path = plot_folder, plot = p1,
       width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Proportion of sets with positive bycatch:
plot_data = weight_data %>%   
  group_by(year, quarter, id_set) %>% 
  summarise(n_sp_catch = length(which(bycatch > 0)))
plot_data = plot_data %>% group_by(year, quarter) %>% 
  summarise(prop_bycatch = (length(which(n_sp_catch > 0))/n())*100)

p1 = ggplot(data = plot_data, aes(x = factor(year), y = prop_bycatch)) +
  geom_col(aes(fill = quarter), position = 'dodge') +
  scale_fill_viridis_d() +
  xlab(NULL) + ylab('Sets with bycatch (%)') +
  theme(legend.position = 'bottom')
ggsave(paste0('obs_prop-pos_set', img_type), path = plot_folder, plot = p1,
       width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Plot env variable on map by year:

# Fishing locations OBS by season/month:
plot_data = obsPoints %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
this_var = "trop_catch"

# Make plot:
p1 = ggplot(plot_data) + geom_sf(aes(color = get(this_var)), size = 0.5, alpha = 0.5)
p1 = add_sf_map(p1)
p1 = p1 + scale_color_viridis() + 
  theme(legend.position = 'bottom') + 
  guides(color=guide_legend(title=this_var)) +
  facet_wrap(~ year)
ggsave(paste0('obs_map_', this_var, img_type), path = plot_folder, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Plot relationship between covariates with bycatch log transformed (only positives):

plot_data = weight_data %>% dplyr::filter(bycatch > 0) %>% 
              mutate(logBycatch = log(bycatch))

p2 = ggplot(data = plot_data, aes(x = sst, y = logBycatch)) +
    geom_point() +
    xlab("sst") + ylab('log(bycatch)') +
    geom_smooth(se = FALSE) +
    facet_wrap(~ sp_name)
ggsave(filename = paste0('sst_rel', img_type), path = plot_folder, 
         plot = p2, width = img_width, height = 160, units = 'mm', dpi = img_res)
