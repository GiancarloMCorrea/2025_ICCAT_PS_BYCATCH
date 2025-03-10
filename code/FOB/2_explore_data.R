rm(list = ls())
# -------------------------------------------------------------------------
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
source('code/parameters_for_plots.R')
source('code/aux_functions.R')

# Select species and school type:
this_type = 'FOB'

# Define folder to save results
plot_folder = file.path('figures', this_type)
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Data folder:
data_folder = file.path('data', this_type)
plot_dir = file.path('figures', this_type)

# -------------------------------------------------------------------------
# Read data in:
weight_data = readRDS(file = file.path(data_folder, 'weight_data.rds'))
numbers_data = readRDS(file = file.path(data_folder, 'numbers_data.rds'))
MyGridSets = readRDS(file = file.path(data_folder, 'extraRegion_tinyVAST.rds'))
load(file.path(data_folder, 'effPoints.RData'))
load(file.path(data_folder, 'obsPoints.RData'))
obsPoints = obsPoints %>% dplyr::filter(year >= 2013)

# -------------------------------------------------------------------------
# Make figures dominance sp groups:
n_sp = 25 # first N species

# Weight:
plot_data = weight_data %>% group_by(Year, sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = plot_data %>% group_by(sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = arrange(cumsp_data, desc(Catch))
cumsp_data = cumsp_data %>% dplyr::filter(Catch > 0)
write.csv(cumsp_data, file = file.path(data_folder, 'cumsp_data_weight.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>% mutate(sp_name = factor(sp_name, levels = cumsp_data$sp_name))

p1 = ggplot(data = plot_data, aes(x = sp_name, y = Catch)) +
  geom_col() + 
  xlab(NULL) + ylab('Weight (kg)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ Year)
ggsave(paste0('sp_dom_weight', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# Numbers:
plot_data = numbers_data %>% group_by(Year, sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = plot_data %>% group_by(sp_name) %>% summarise(Catch = sum(Catch))
cumsp_data = arrange(cumsp_data, desc(Catch))
cumsp_data = cumsp_data %>% dplyr::filter(Catch > 0)
write.csv(cumsp_data, file = file.path(data_folder, 'cumsp_data_numbers.csv'), row.names = FALSE)
cumsp_data = cumsp_data[1:n_sp, ]
plot_data = plot_data %>% dplyr::filter(sp_name %in% cumsp_data$sp_name) %>% mutate(sp_name = factor(sp_name, levels = cumsp_data$sp_name))

p1 = ggplot(data = plot_data, aes(x = sp_name, y = log(Catch+1))) +
  geom_col() + 
  xlab(NULL) + ylab('log(Numbers)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  facet_wrap(~ Year)
ggsave(paste0('sp_dom_numbers', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 140, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Fishing locations OBS by season/month:
plot_data = obsPoints %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>% 
              mutate(month = as.numeric(month),
                     season = ifelse(month %in% 1:3, '1',
                                    ifelse(month %in% 4:6, '2',
                                          ifelse(month %in% 7:9, '3', '4'))))

# Make plot:
p1 = ggplot(plot_data) + geom_sf(aes(color = season), size = 0.5, alpha = 0.35)
p1 = add_sf_map(p1)
p1 = p1 + scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom') + facet_wrap(~ year)
ggsave(paste0('obs_map_sets_seas', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 150, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Number bycatch sp per set
plot_data = weight_data %>% group_by(Year, month, id_set) %>% 
  summarise(bycatch_catch = sum(Catch), n_sp_catch = length(which(Catch > 0)))
plot_data = plot_data %>% mutate(month = as.numeric(month),
                                 season = ifelse(month %in% 1:3, '1',
                                                 ifelse(month %in% 4:6, '2',
                                                        ifelse(month %in% 7:9, '3', '4'))))

p1 = ggplot(plot_data, aes(x = factor(Year), y = n_sp_catch)) + 
  geom_boxplot(aes(fill = season), width = 0.8) +
  scale_fill_brewer(palette = 'Set1') +
  xlab(NULL) + ylab('N-sp/set') +
  theme(legend.position = 'bottom')
ggsave(paste0('obs_n-sp_set', img_type), path = plot_folder, plot = p1,
       width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Proportion of sets with positive bycatch:
plot_data = weight_data %>% group_by(Year, month, id_set) %>% 
  summarise(n_sp_catch = length(which(Catch > 0)))
plot_data = plot_data %>% mutate(month = as.numeric(month),
                                 season = ifelse(month %in% 1:3, '1',
                                                 ifelse(month %in% 4:6, '2',
                                                        ifelse(month %in% 7:9, '3', '4'))))
plot_data = plot_data %>% group_by(Year, season) %>% 
  summarise(prop_bycatch = (length(which(n_sp_catch > 0))/n())*100)

p1 = ggplot(data = plot_data, aes(x = factor(Year), y = prop_bycatch)) +
  geom_col(aes(fill = season), position = 'dodge') +
  scale_fill_brewer(palette = 'Set1') +
  xlab(NULL) + ylab('Sets with bycatch (%)') +
  theme(legend.position = 'bottom')
ggsave(paste0('obs_prop-pos_set', img_type), path = plot_folder, plot = p1,
       width = img_width*0.75, height = 90, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Plot N sets per grid:
plot_data = MyGridSets %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)

# Make plot:
p1 = ggplot(plot_data) + geom_sf(aes(color = n_sets), size = 0.5)
p1 = add_sf_map(p1)
p1 = p1 + scale_color_viridis() + 
  theme(legend.position = 'bottom') + facet_wrap(~ Year)
ggsave(paste0('eff_map_grid-sets', img_type), path = plot_folder, plot = p1,
       width = img_width, height = 150, units = 'mm', dpi = img_res)
